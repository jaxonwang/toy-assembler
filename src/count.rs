use fastq::{Parser, Record, RefRecord};
use serde::{Deserialize, Serialize};
use std::boxed::Box;
use std::cell::RefCell;
use std::collections::{HashMap, HashSet};
use std::hash::Hash;
use std::hash::Hasher;
use std::io;
use std::rc::Rc;
use std::sync::{Arc, Mutex};
use std::time::Instant;

extern crate chrono;
extern crate crossbeam;
extern crate fastq;
extern crate fern;
extern crate log;
extern crate mpi;
extern crate serde;

use log::{debug, info};
use mpi::topology::SystemCommunicator;
use mpi::traits::*;
use mpi::Threading;

mod msg;
use self::msg::Packed;

#[derive(Eq, PartialEq, Clone, Copy)]
enum Direction {
    Left,
    Right,
}

use self::Direction::*;

impl Direction {
    fn reverse(&self) -> Direction {
        match self {
            Left => Right,
            Right => Left,
        }
    }
}

type CountNumber = u32;
type KMer = Vec<u8>;

fn kmer_str(k: &KMer) -> String {
    String::from_utf8_lossy(&k[..]).to_string()
}

#[derive(Debug, Copy, Clone, Eq, Serialize, Deserialize)]
struct KMerData {
    adj: [u32; 8], // number of edges adjcent to a kmer with base [leftA, leftT, leftG, leftC, rightG, rC, rA, rT];
    count: CountNumber,
}
impl PartialEq for KMerData {
    fn eq(&self, other: &Self) -> bool {
        self.count == other.count && self.adj == other.adj
    }
}

type CountTable = HashMap<KMer, KMerData>;

const BASE_A: u8 = b'A';
const BASE_T: u8 = b'T';
const BASE_C: u8 = b'C';
const BASE_G: u8 = b'G';

#[allow(dead_code)]
fn partition_by_borders<O: Ord>(kmer: &O, borders: &[O]) -> usize {
    // return i if kmer lies in (borders[i-1], borders[i]]
    use std::cmp::Ordering::*;
    let mut i = 0usize;
    let mut j = borders.len();
    let mut middle = j;
    while i < j {
        middle = (i + j) / 2;
        match kmer.cmp(&borders[middle]) {
            Greater => i = middle + 1,
            Equal => break,
            Less => j = middle,
        };
    }
    if i == middle + 1 && i >= j {
        return i;
    }
    middle
}

fn _partition_by_hash(kmer: &KMer, partition_num: usize, window_size: usize) -> usize {
    let dohash = |window: &[u8]| {
        let mut hasher = std::collections::hash_map::DefaultHasher::new();
        window.hash(&mut hasher);
        hasher.finish()
    };

    let mut hashes = vec![];
    for start in 0..kmer.len() - window_size {
        hashes.push(dohash(&kmer[start..start + window_size]));
    }
    (*hashes.iter().min().unwrap() % partition_num as u64) as usize
}

fn partition_by_hash(kmer: &KMer, partition_num: usize) -> usize {
    let window_size = (kmer.len() as f64).sqrt() as usize; //NOTE:heuristicly choose sqrt
    if window_size < 1 {
        panic!("kmer too short for locality sensitive hashing")
    }
    _partition_by_hash(kmer, partition_num, window_size)
}

mod mpiconst {
    pub const ROOT_RANK: mpi::topology::Rank = 0;
    pub const DONE_MARK: &[u8] = "Done".as_bytes();
}

fn update_count_table(count_table: &mut CountTable, kmer: KMer, kmerdata: KMerData) {
    // merge data of the same kmer
    match count_table.get_mut(&kmer) {
        Some(kd) => {
            kd.count += kmerdata.count;
            for i in 0..kd.adj.len() {
                kd.adj[i] += kmerdata.adj[i];
            }
        }
        None => {
            count_table.insert(kmer, kmerdata);
        }
    };
}

fn get_complement_base(base: u8) -> u8 {
    match base {
        BASE_A => BASE_T,
        BASE_T => BASE_A,
        BASE_C => BASE_G,
        BASE_G => BASE_C,
        other => panic!(format!("unknown base: {}", other)),
    }
}

fn get_reverse_complement(kmer: &KMer) -> KMer {
    kmer.iter()
        .rev()
        .map(|x| get_complement_base(*x))
        .clone()
        .collect()
}

fn get_reverse_complement_adj(adj: &[u32; 8]) -> [u32; 8] {
    // number of edges adjcent to a kmer with base [leftA, leftT, leftG, leftC, rightG, rC, rA, rT];
    let mut a = [0u32; 8];
    for i in 0..adj.len() {
        a[7 - i] = adj[i];
    }
    a
}

fn adj_position(direction: Direction, base: u8) -> usize {
    match direction {
        Left => match base {
            BASE_A => 0,
            BASE_T => 1,
            BASE_G => 2,
            BASE_C => 3,
            _ => panic!("bad base"),
        },
        Right => match base {
            BASE_G => 4,
            BASE_C => 5,
            BASE_A => 6,
            BASE_T => 7,
            _ => panic!("bad base"),
        },
    }
}

fn adj_base(pos: usize) -> u8 {
    match pos {
        0 => BASE_A,
        1 => BASE_T,
        2 => BASE_G,
        3 => BASE_C,
        4 => BASE_G,
        5 => BASE_C,
        6 => BASE_A,
        7 => BASE_T,
        _ => panic!("bad pos"),
    }
}

fn adj_new_edge_with_base(adj: &mut [u32; 8], direction: Direction, base: u8) {
    adj[adj_position(direction, base)] += 1;
}

fn kmer_side_neighbor(kmer: &KMer, data: &KMerData, direction: Direction) -> Vec<(KMer, u8)> {
    // return the out edge bases
    let range = match direction {
        Left => 0..4,
        Right => 4..8,
    };
    let mut neighbors: Vec<(KMer, u8)> = vec![];
    for i in range {
        if data.adj[i] > 0 {
            let mut neighbor_kmer: KMer = vec![0u8; kmer.len()];
            if let Left = direction {
                neighbor_kmer[1..].copy_from_slice(&kmer[..kmer.len() - 1]);
                neighbor_kmer[0] = adj_base(i);
            } else {
                neighbor_kmer[..kmer.len() - 1].copy_from_slice(&kmer[1..]);
                neighbor_kmer[kmer.len() - 1] = adj_base(i);
            }
            neighbors.push((neighbor_kmer, adj_base(i)));
        }
    }
    neighbors
}

trait KMerIndexing<'a, V>: std::ops::Index<&'a KMer, Output = V> {
    fn get_adj_mut(&mut self, kmer: &KMer) -> Option<&mut [u32; 8]>;
}

impl KMerIndexing<'_, KMerData> for CountTable {
    fn get_adj_mut(&mut self, kmer: &KMer) -> Option<&mut [u32; 8]> {
        self.get_mut(kmer).map(|x| &mut x.adj)
    }
}
impl KMerIndexing<'_, KMerExtendData> for KMerExtendTable {
    fn get_adj_mut(&mut self, kmer: &KMer) -> Option<&mut [u32; 8]> {
        self.get_mut(kmer).map(|x| &mut x.kmer_data.adj)
    }
}

fn remove_edges_connect_to_kmer<'a, T, V>(
    table: &mut T,
    kmer: &KMer,
    kmer_data: &KMerData,
    direction: Direction,
) where
    T: KMerIndexing<'a, V>,
{
    let neighbors = kmer_side_neighbor(kmer, kmer_data, direction);
    for (nk, _) in neighbors {
        // remove neighbors edges
        let mut edge = match direction {
            Left => kmer[kmer.len() - 1], // edge is at in k's sides
            Right => kmer[0],
        };
        let (nk, changed) = get_canonical(nk);
        let neightbor_end;
        if changed {
            edge = get_complement_base(edge);
            neightbor_end = direction;
        } else {
            neightbor_end = direction.reverse();
        }
        match table.get_adj_mut(&nk) {
            Some(adj) => adj[adj_position(neightbor_end, edge)] = 0,
            None => (),
        }
    }
}

fn adj_one_side_degree(adj: &[u32; 8], direction: Direction) -> usize {
    let range = match direction {
        Left => 0..4,
        Right => 4..8,
    };
    adj[range]
        .iter()
        .filter_map(|c| if *c > 0 { Some(1) } else { None })
        .sum()
}

fn adj_degree(adj: &[u32; 8]) -> usize {
    adj.iter()
        .filter_map(|c| if *c > 0 { Some(1) } else { None })
        .sum()
}

fn get_canonical(kmer: KMer) -> (KMer, bool) {
    // bool is true if there is a change
    use std::cmp::Ordering::*;
    for i in 0..kmer.len() / 2 + 1 {
        match kmer[i].cmp(&get_complement_base(kmer[kmer.len() - 1 - i])) {
            Less => return (get_reverse_complement(&kmer), true),
            Greater => return (kmer, false),
            Equal => continue,
        }
    }
    panic!("KMer is plalindrome")
}

struct KMerCounter {
    kmer_len: usize,
    min_cov: usize,
}

impl KMerCounter {
    pub fn count_kmers_from_fastq(
        &self,
        count_table: &mut CountTable,
        r: impl io::Read,
    ) -> io::Result<()> {
        let parser = Parser::new(r);

        info!("reading file and generating kmers");
        parser
            .each(|record| self.all_reads(count_table, record))
            .map(|_| ()) // I dont' care the bool in the result
    }

    fn all_reads(&self, count_table: &mut CountTable, record: RefRecord) -> bool {
        let read: &[u8] = Record::seq(&record);

        // drop too short read
        if read.len() < self.kmer_len {
            return true;
        }
        // drop read with unknown base
        for b in read.iter() {
            // illumina unknown base
            if *b == b'N' || *b == b'.' {
                return true;
            }
        }

        // TODO: further optimization
        for start_pos in 0..(read.len() - self.kmer_len + 1) {
            let kmer = Vec::from(&read[start_pos..start_pos + self.kmer_len]);
            let mut adj = [0u32; 8];
            if start_pos > 0 {
                // edge to left
                adj_new_edge_with_base(&mut adj, Left, read[start_pos - 1]);
            }
            if start_pos + self.kmer_len < read.len() {
                //edge to right
                adj_new_edge_with_base(&mut adj, Right, read[start_pos + self.kmer_len]);
            }

            let (canonical, changed) = get_canonical(kmer);
            if changed {
                adj = get_reverse_complement_adj(&adj);
            }
            let kmer_data = KMerData { adj, count: 1 };

            update_count_table(count_table, canonical, kmer_data);
        }
        true
    }

    pub fn distributed_count(&self, count_table: &CountTable, world: &SystemCommunicator) {
        let p = world.size();

        let mut new_table: CountTable = HashMap::new();
        let local_table = Arc::new(Mutex::new(HashMap::new())); // stores entries partitioned to local

        let mut mpi_time = 0f64;

        info!("start counting");

        crossbeam::scope(|scope| {
            // dispatch kemrs by borders
            scope.spawn(|_| self.dispatcher(count_table, world, local_table.clone()));

            let mut done = 0;
            while done != p - 1 {
                debug!("try_to receive");
                let start = Instant::now();
                let (msg, s) = world.any_process().receive_vec::<u8>();
                mpi_time += (Instant::now() - start).as_secs_f64();
                if msg == mpiconst::DONE_MARK {
                    done += 1;
                    debug!("{} Done!", s.source_rank());
                    continue;
                }
                let mut pos = 0usize;
                // unpack the msg one by one
                while pos != msg.len() {
                    let (ret, size) = KMer::unpack(&msg[pos..]);
                    let kmer = ret.unwrap();
                    pos += size;
                    let (ret, size) = KMerData::unpack(&msg[pos..]);
                    let kmerdata = ret.unwrap();
                    pos += size;
                    update_count_table(&mut new_table, kmer, kmerdata);
                }
            }
        })
        .unwrap();
        // update from local
        let guard = local_table.lock().unwrap();
        for (k, v) in guard.iter() {
            update_count_table(&mut new_table, k.clone(), *v);
        }
        info!(
            "Count Done: kmers number {}, mpi recv time {}",
            new_table.len(),
            mpi_time
        );
        // for (k, v) in new_table.iter() {
        //     println!("{}: {:?}", String::from_utf8_lossy(k), v);
        // }
    }

    fn dispatcher(
        &self,
        count_table: &CountTable,
        world: &SystemCommunicator,
        local_table: Arc<Mutex<CountTable>>,
    ) {
        let myrank = world.rank();
        let buffer_size = 128usize;
        let max_num_in_buf = 1usize;
        let mut num_in_buf = 0usize;
        let mut mpi_time = 0f64;
        let mut buffers: Vec<Vec<u8>> = vec![];
        for _ in 0..world.size() {
            buffers.push(Vec::<u8>::with_capacity(buffer_size));
        }

        for (k, v) in count_table.iter() {
            let dst = partition_by_hash(k, world.size() as usize);
            if dst == myrank as usize {
                // if local set localtable
                let mut guard = local_table.lock().unwrap();
                update_count_table(&mut *guard, k.clone(), *v);
                continue;
            }
            if num_in_buf == max_num_in_buf {
                let dst_process = world.process_at_rank(dst as mpi::topology::Rank);
                let start = Instant::now();
                dst_process.send(&buffers[dst][..]); // send full buf
                mpi_time += (Instant::now() - start).as_secs_f64();
                // empty buffer
                buffers[dst].truncate(0);
            } else {
                k.pack(&mut buffers[dst]);
                v.pack(&mut buffers[dst]);
                num_in_buf += 1;
            }

            // debug!(
            //     "send {} to {}",
            //     String::from_utf8_lossy(&msg),
            //     dst
            // );
        }
        for dst in 0..world.size() {
            if dst != myrank {
                let start = Instant::now();
                world.process_at_rank(dst).send(&buffers[dst as usize][..]);
                world.process_at_rank(dst).send(mpiconst::DONE_MARK);
                mpi_time += (Instant::now() - start).as_secs_f64();
            }
        }
        info!("dispatch done, mpi send time {}", mpi_time);
    }

    fn kmer_coverage_filter(&self, count_table: &mut CountTable) {
        let mut low_cov_kmers: Vec<(KMer, KMerData)> = vec![];
        for (k, data) in count_table.iter() {
            if (data.count as usize) < self.min_cov {
                low_cov_kmers.push((k.clone(), data.clone()));
            }
        }
        for (k, data) in low_cov_kmers {
            for direction in [Left, Right].iter() {
                remove_edges_connect_to_kmer(count_table, &k, &data, *direction);
            }
            count_table.remove(&k);
        }
    }
}

struct LinearPath {
    // two sides
    path: Vec<u8>,
    coverage: u32,
}

impl LinearPath {
    fn get_ended_kmer(&self, direction: Direction, kmer_len: usize) -> KMer {
        let path_len = self.path.len();
        if path_len < kmer_len {
            panic!(format!(
                "path len: {} shorter than kmer len {}",
                path_len, kmer_len
            ));
        }
        match direction {
            Left => self.path[..kmer_len].to_vec(),
            Right => self.path[path_len - kmer_len..path_len].to_vec(),
        }
    }
}

fn adj_str(adj: &[u32; 8]) -> String {
    let sj = adj
        .iter()
        .enumerate()
        .map(|(i, v)| match *v > 0 {
            true => adj_base(i),
            false => b'-',
        })
        .collect();
    kmer_str(&sj)
}

#[derive(Serialize, Deserialize)]
struct CompactNode {
    node: Vec<u8>,
    coverage: u32,
}

impl CompactNode {
    fn new(p: &[u8], cov: u32) -> Self {
        CompactNode {
            node: p.to_vec(),
            coverage: cov,
        }
    }
}

#[derive(Serialize, Deserialize)]
struct SubGraph {
    //CSC format
    nodes: Vec<CompactNode>,    //sorted
    adjs: Vec<Vec<[usize; 8]>>, // adj[0] -> list of index node 0 connect to
}

impl SubGraph {
    fn new_single_kmer(k: &KMer, cov: u32) -> Self {
        SubGraph {
            nodes: vec![CompactNode::new(&k[..], cov)],
            adjs: vec![vec![]],
        }
    }
}

#[derive(Clone)]
struct GraphConnect {
    ptr: Arc<SubGraph>,
    node_idx: usize,      // which node connect to?
    direction: Direction, // is_left -> true, which side connect to?
}

#[derive(Clone)]
struct LinearConnect {
    ptr: Arc<LinearPath>,
    direction: Direction,
}

enum Connect {
    SomeGraph(GraphConnect),
    SomeLinear(LinearConnect),
}

impl Clone for Connect {
    fn clone(&self) -> Self {
        match self {
            Connect::SomeGraph(a) => Connect::SomeGraph(a.clone()),
            Connect::SomeLinear(a) => Connect::SomeLinear(a.clone()),
        }
    }
}

impl Connect {
    fn new_some_graph(ptr: Arc<SubGraph>, node_idx: usize, direction: Direction) -> Self {
        Self::SomeGraph(GraphConnect {
            ptr,
            node_idx,
            direction,
        })
    }
    fn new_some_linear(ptr: Arc<LinearPath>, direction: Direction) -> Self {
        Self::SomeLinear(LinearConnect { ptr, direction })
    }
    #[allow(dead_code)]
    fn unwrap_graph(self) -> Arc<SubGraph> {
        match self {
            Self::SomeGraph(g) => g.ptr,
            _ => panic!("Connect is not graph"),
        }
    }
    #[allow(dead_code)]
    fn unwrap_linear(self) -> Arc<LinearPath> {
        match self {
            Self::SomeLinear(l) => l.ptr,
            _ => panic!("Connect is not linear"),
        }
    }
    fn borrow_linear(&self) -> &LinearConnect {
        match self {
            Self::SomeLinear(ref l) => l,
            _ => panic!("Connect is not linear"),
        }
    }
}

struct KMerExtendData {
    kmer_data: KMerData,
    connect: Connect,
}

type KMerExtendTable = HashMap<KMer, KMerExtendData>;

fn local_linear_path(kmers_from_counter: CountTable, kmer_len: usize) -> KMerExtendTable {
    let mut visited = HashSet::new();
    let mut extend_table = HashMap::new();
    info!("start local linear path generation.");

    // try connect kmers
    for (k, v) in kmers_from_counter.iter() {
        // extend right
        let is_linear = |x| adj_degree(&kmers_from_counter[x].adj) <= 2;
        if visited.contains(k) {
            continue;
        }

        let mut travel = |is_left: Direction| -> Vec<(u8, u32)> {
            let mut side_path: Vec<(u8, u32)> = vec![]; // base, count
            let mut current: &KMer = k;
            let mut direction = is_left;
            // extend to one direction
            while !visited.contains(current) && is_linear(current) {
                visited.insert(current.clone()); // now k is linear && not visited. mark as visited.
                let base;
                if let Left = is_left {
                    if direction == is_left {
                        base = current[0];
                    } else {
                        base = get_complement_base(current[kmer_len - 1]);
                    }
                } else {
                    if direction == is_left {
                        base = current[kmer_len - 1];
                    } else {
                        base = get_complement_base(current[0]);
                    }
                }
                // TODO: clone here
                let current_data: &KMerData = kmers_from_counter.get(current).unwrap();
                side_path.push((base, current_data.count));
                let mut neighbor = kmer_side_neighbor(current, current_data, direction);
                match neighbor.len() {
                    0 => break, // this is an end
                    1 => (),
                    n => panic!(format!("impossible degree {}", n)),
                }
                let (neighbor, _) = neighbor.pop().unwrap();
                let (ref extend, changed) = get_canonical(neighbor);
                if changed {
                    // change to it reverse complement, to go the same direction, need change
                    direction = direction.reverse();
                }
                if kmers_from_counter.contains_key(extend) {
                    // is local
                    current = kmers_from_counter.get_key_value(extend).unwrap().0;
                }
            }
            if let Right = is_left {
                visited.remove(k); // remove for left travel enabled
            }
            side_path
        };

        let rpath = travel(Right);
        let lpath = travel(Left); //extend left

        if lpath.len() == 0 {
            // nothing to do for this kmer, store it as SubGraph
            let sub = Arc::new(SubGraph::new_single_kmer(k, v.count));
            extend_table.insert(
                k.clone(),
                KMerExtendData {
                    kmer_data: v.clone(),
                    connect: Connect::new_some_graph(sub, 0, Left),
                },
            );
            continue;
        }

        let mut path: Vec<u8> = vec![0; lpath.len() + rpath.len() - 2 + k.len()];
        // init as middle
        let mut coverage = lpath.last().unwrap().1;
        let lpath_len = lpath.len();
        // must use k(middle) as start, orthewise we don't know direction to continue
        path[lpath_len - 1..lpath_len + kmer_len - 1].copy_from_slice(&k[..]);

        // connect left
        for i in 1..lpath_len {
            let (base, count) = lpath[i];
            path[lpath_len - i - 1] = base;
            coverage += count;
        }
        // connect right
        for i in 1..rpath.len() {
            let (base, count) = rpath[i];
            path[lpath_len + kmer_len - 2 + i] = base;
            coverage += count;
        }
        path = get_canonical(path).0;

        let left_end_kmer = get_canonical(path[..kmer_len].to_vec()).0;
        let right_end_kmer = get_canonical(path[path.len() - kmer_len..path.len()].to_vec()).0;
        let lp_obj = Arc::new(LinearPath { path, coverage });
        let l_kmer_data = kmers_from_counter[&left_end_kmer].clone();
        let r_kmer_data = kmers_from_counter[&right_end_kmer].clone();
        extend_table.insert(
            // insert look up for left
            left_end_kmer,
            KMerExtendData {
                kmer_data: l_kmer_data,
                connect: Connect::new_some_linear(lp_obj.clone(), Left),
            },
        );
        extend_table.insert(
            // for right
            right_end_kmer,
            KMerExtendData {
                kmer_data: r_kmer_data,
                connect: Connect::new_some_linear(lp_obj.clone(), Right),
            },
        );
    }
    info!("linear path generation done.");
    extend_table
}

struct ContigGenerator {
    kmer_len: usize,
    world: Arc<SystemCommunicator>,
    tip_minimal_coverage: f64,
    tip_minimal_length: usize,
    tip_travel_stop_length: usize,
}

#[derive(Serialize, Deserialize, Clone)]
enum TipRemovalMessageType {
    CallingFrame,
    Return,
    Done,
}

#[derive(Serialize, Deserialize, Clone)]
struct TipRemovalCallingFrame {
    callers: Vec<mpi::topology::Rank>,
    from_base: u8, // if 0 means it is the start of travel, that's a local request
    kmer: KMer,    // the kmer to visit next
    path_len: usize,
    path_coverage: u32, // count for the path, not count/len
}

#[derive(Serialize, Deserialize, Clone)]
struct TipRemovalReturn {
    callers: Vec<mpi::topology::Rank>, // who to returen
    kmer: KMer,                        // this
    legal: bool,
}

enum TipRemovalMessageEnum {
    CallingFrame(Box<TipRemovalCallingFrame>),
    Return(Box<TipRemovalReturn>),
    Local,
}

impl ContigGenerator {
    fn is_legal_tip(&self, path_len: usize, path_coverage: u32) -> bool {
        // only judged at the end of a kmer
        path_len >= self.tip_minimal_length
            && (path_coverage as f64 / path_len as f64) >= self.tip_minimal_coverage
    }

    fn tips_travel_kmer(
        &self,
        kmer_table: Arc<KMerExtendTable>,
        frame: TipRemovalCallingFrame,
        mut illegals: Arc<Mutex<Vec<KMer>>>,
    ) {
        use TipRemovalMessageEnum::*;
        match self.tips_travel_kmer_logic(kmer_table, frame, illegals) {
            CallingFrame(cf) => {
                let mut msg = vec![];
                TipRemovalMessageType::CallingFrame.pack(&mut msg);
                cf.pack(&mut msg);
                let dst = partition_by_hash(&cf.kmer, self.world.size() as usize);
                self.world
                    .process_at_rank(dst as mpi::topology::Rank)
                    .send(&msg[..]);
            }
            Return(rf) => {
                let mut msg = vec![];
                TipRemovalMessageType::Return.pack(&mut msg);
                rf.pack(&mut msg);
                let dst = partition_by_hash(&rf.kmer, self.world.size() as usize);
                self.world
                    .process_at_rank(dst as mpi::topology::Rank)
                    .send(&msg[..]);
            }
            Local => (),
        };
    }

    fn tips_travel_kmer_logic(
        &self,
        kmer_table: Arc<KMerExtendTable>,
        frame: TipRemovalCallingFrame,
        mut illegals: Arc<Mutex<Vec<KMer>>>,
    ) -> TipRemovalMessageEnum {
        use TipRemovalMessageEnum::*;
        assert!(
            kmer_table.contains_key(&frame.kmer),
            "You asked the wrong peer!"
        );
        let ref extended_data = kmer_table[&frame.kmer];

        let assemble_ret = |this_path_len, this_path_count| -> TipRemovalMessageEnum {
            let mut ret_kmer = vec![frame.from_base];
            ret_kmer.extend_from_slice(&frame.kmer[..self.kmer_len - 1]);
            let mut callers = frame.callers.clone();
            Return(Box::new(TipRemovalReturn {
                callers,
                kmer: ret_kmer,
                legal: self.is_legal_tip(
                    frame.path_len + this_path_len,
                    frame.path_coverage + this_path_count,
                ),
            }))
        };

        let mark_remove_if_illegal = |legal, kmer, mut illegals: Arc<Mutex<Vec<KMer>>>| {
            if legal {
                let mut list = illegals.lock().unwrap();
                list.push(kmer);
            }
        };

        if adj_degree(&extended_data.kmer_data.adj) != 2 {
            // reach end
            assert!(frame.from_base != 0, "the from base can not be 0! should not get a local request when this kmer is not linear!!");
            let ret = assemble_ret(0, 0);
            if let Return(ref ret) = ret {
                // mark as to remove
                mark_remove_if_illegal(ret.legal, frame.kmer.clone(), illegals.clone());
            }
            return ret;
        };

        let linear_connect = extended_data.connect.borrow_linear();
        let ref linear_path = linear_connect.ptr;
        let linear_path_len = linear_path.path.len();
        if linear_path_len + frame.path_len >= self.tip_travel_stop_length {
            // if is long engouh, stop traveling, return here to reduce one cross node travel
            if frame.from_base != 0 {
                // is from remove, assemble a message
                let ret = assemble_ret(linear_path_len, linear_path.coverage);
                if let Return(ref ret) = ret {
                    // mark as to remove
                    mark_remove_if_illegal(ret.legal, frame.kmer.clone(), illegals.clone());
                }
                return ret;
            } else {
                let legal = self.is_legal_tip(linear_path_len, linear_path.coverage);
                mark_remove_if_illegal(legal, frame.kmer.clone(), illegals.clone());
                return Local;
            }
        }

        let other_end_kmer =
            ContigGenerator::tip_travel_to_end_kmer(kmer_table.clone(), &frame.kmer);
        let (next_kmer, changed) =
            ContigGenerator::tip_travel_end_to_next_hop(kmer_table.clone(), &other_end_kmer);

        let mut callers = frame.callers.clone();
        callers.push(self.world.rank());
        let mut from_base = match linear_connect.direction {
            Left => other_end_kmer[self.kmer_len - 1],
            Right => other_end_kmer[0],
        };
        if changed {
            from_base = get_complement_base(from_base);
        }
        let call_frame = TipRemovalCallingFrame {
            callers,
            from_base,
            kmer: next_kmer,
            path_len: linear_path_len + frame.path_len,
            path_coverage: linear_path.coverage + frame.path_coverage,
        };
        TipRemovalMessageEnum::CallingFrame(Box::new(call_frame))
    }

    fn tip_travel_to_end_kmer(kmer_table: Arc<KMerExtendTable>, kmer: &KMer) -> KMer {
        // in a local linear_path, travel from one end to the other end
        let linear_connect = kmer_table[kmer].connect.borrow_linear();
        let ref linear_path = linear_connect.ptr;
        linear_path.get_ended_kmer(linear_connect.direction.reverse(), kmer.len())
    }

    fn tip_travel_end_to_next_hop(kmer_table: Arc<KMerExtendTable>, kmer: &KMer) -> (KMer, bool) {
        // kmer is one end of a local linear_path, to get to it's contigous remote kmer
        let ref kmer_data = kmer_table[kmer].kmer_data;
        let linear_connect = kmer_table[kmer].connect.borrow_linear();
        let direction = linear_connect.direction;
        let next_hops = kmer_side_neighbor(kmer, kmer_data, direction);
        assert!(next_hops.len() == 1, "out edge is not 1");
        get_canonical(next_hops[0].0.clone())
    }

    fn remove_low_coverage_tips(&self, kmer_table: KMerExtendTable) {
        // must be call after all local linear path are generated
        let kmer_table = Arc::new(kmer_table);

        let mut tip_path: Vec<(&KMer, Direction)> = vec![];
        for (k, data) in kmer_table.iter() {
            let degree_left = adj_one_side_degree(&data.kmer_data.adj, Left);
            let degree_right = adj_one_side_degree(&data.kmer_data.adj, Right);
            let direction = match (degree_left, degree_right) {
                (0, 0) => panic!("imposible kmer without degree"),
                (0, 1) => Right,
                (1, 0) => Left,
                _ => continue, // not a tip
            };
            tip_path.push((k, direction));
        }
        let low_quality_tips: Vec<Vec<KMer>> = vec![];

        let mut illegals = Arc::new(Mutex::new(vec![]));

        crossbeam::scope(|scope| {
            scope.spawn(|_| self.tips_remove_msg_receiver(kmer_table.clone(), illegals.clone()));

            for (k, direction) in tip_path {
                // NOTE: this can be paralleled
                self.tips_travel_kmer(
                    kmer_table.clone(),
                    TipRemovalCallingFrame {
                        callers: vec![],
                        from_base: 0,
                        kmer: k.clone(),
                        path_len: 0,
                        path_coverage: 0,
                    },
                    illegals.clone(),
                );
            }
            for i in 0..self.world.size() {
                let mut msg = vec![];
                TipRemovalMessageType::Done.pack(&mut msg);
                self.world.process_at_rank(i).send(&msg[..]);
            }
        })
        .unwrap();
    }

    fn tips_remove_msg_receiver(
        &self,
        kmer_table: Arc<KMerExtendTable>,
        mut illegals: Arc<Mutex<Vec<KMer>>>,
    ) {
        crossbeam::scope(|scope| {
            use TipRemovalMessageType::*;
            let mut done = 0;
            let world_size = self.world.size() as usize;

            while done != world_size {
                let (msg, _) = self.world.any_process().receive_vec();
                let (msg_type, size) = TipRemovalMessageType::unpack(&msg[..]);
                let msg_type = msg_type.unwrap();
                match msg_type {
                    CallingFrame => {
                        let (frame, _) = TipRemovalCallingFrame::unpack(&msg[size..]);
                        let frame = frame.unwrap();
                        scope.spawn(|_| {
                            self.tips_travel_kmer(kmer_table.clone(), frame, illegals.clone())
                        });
                    }
                    Return => {
                        let (ret, _) = TipRemovalReturn::unpack(&msg[size..]);
                        let ret = ret.unwrap();
                        // mark local removal
                        if !ret.legal {
                            debug_assert!(kmer_table.contains_key(&ret.kmer));
                            let mut list = illegals.lock().unwrap();
                            list.push(ret.kmer.clone());
                        }
                        if ret.callers.len() > 0 {
                            // continue return
                            let ret = ret.clone();
                            let kmer_table = kmer_table.clone();
                            scope.spawn(move |_| {
                                let end = ContigGenerator::tip_travel_to_end_kmer(
                                    kmer_table.clone(),
                                    &ret.kmer,
                                );
                                let (next_hop, _) = ContigGenerator::tip_travel_end_to_next_hop(
                                    kmer_table.clone(),
                                    &end,
                                );
                                let dst = partition_by_hash(&next_hop, self.world.size() as usize);

                                let ret_callers = ret.callers[..ret.callers.len() - 1].to_vec();
                                let next_ret = TipRemovalReturn {
                                    callers: ret_callers,
                                    kmer: next_hop,
                                    legal: ret.legal,
                                };
                                let mut msg = vec![];
                                next_ret.pack(&mut msg);

                                self.world
                                    .process_at_rank(dst as mpi::topology::Rank)
                                    .send(&msg[..]);
                            });
                        }
                    }
                    Done => done += 1,

                    Local => panic!("should be no local request!"),
                }
            }
        })
        .unwrap();
    }
}
// pub fn main() {
//     let kmer_len = 30usize;
//
//     let fastq_path = match env::args().nth(1) {
//         None => {
//             println!("usage: a.out file");
//             std::process::exit(1)
//         }
//         Some(p) => p,
//     };
//     let f = fs::File::open(fastq_path).unwrap();
//
//     let mut count_table = HashMap::new();
//     count_kmers_from_fastq(&mut count_table, kmer_len, f).unwrap();
//
//     for (key, v) in count_table.iter() {
//         println!("{}: {}", String::from_utf8_lossy(key), v);
//     }
// }
//

fn setup_logger(process_id: Arc<String>) -> Result<(), fern::InitError> {
    fern::Dispatch::new()
        .format(move |out, message, record| {
            out.finish(format_args!(
                "{} {} {} {} {}",
                chrono::Local::now().format("%m-%d %H:%M:%S.%f"),
                record.target(),
                process_id,
                record.level(),
                message
            ))
        })
        .level(log::LevelFilter::Info)
        .chain(std::io::stdout())
        .apply()?;
    Ok(())
}

fn main() {
    let (universe, threading) = mpi::initialize_with_threading(Threading::Multiple).unwrap();
    assert_eq!(threading, mpi::environment::threading_support());

    let world = universe.world();

    let process_id = Arc::new(format!("P{}", world.rank()));
    setup_logger(process_id).unwrap();

    let counter = KMerCounter {
        kmer_len: 31,
        min_cov: 1,
    };
    //     let in1 = "@
    // AAATTTCC
    // +
    // AAAAAEEE
    // ";
    //     let f = in1.as_bytes();

    use std::env;
    use std::fs;
    let fastq_path = match env::args().nth(1) {
        None => {
            println!("usage: a.out file");
            std::process::exit(1)
        }
        Some(p) => p,
    };
    let f = fs::File::open(fastq_path).unwrap();

    let mut count_table = HashMap::new();
    counter.count_kmers_from_fastq(&mut count_table, f).unwrap();
    counter.distributed_count(&count_table, &world);
    counter.kmer_coverage_filter(&mut count_table);
    local_linear_path(count_table, counter.kmer_len);
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::Rng;

    fn new_kmer(k: &str) -> KMer {
        k.as_bytes().to_vec()
    }

    #[test]
    fn count_kmers_from_fastq_test() {
        let in1 = "@
AAATTGCC
+
AAAAAEEE
";

        let map_kmer = |x: &(&str, u32, [u32; 8])| {
            (
                new_kmer(x.0),
                KMerData {
                    count: x.1,
                    adj: x.2,
                },
            )
        };
        let counter = KMerCounter {
            kmer_len: 6,
            min_cov: 1,
        };
        let mut count_table = HashMap::new();
        counter
            .count_kmers_from_fastq(&mut count_table, in1.as_bytes())
            .unwrap();
        let out1: CountTable = [
            ("GGCAAT", 1, [0, 0, 0, 0, 0, 0, 0, 1]),
            ("CAATTT", 1, [0, 0, 1, 0, 0, 0, 0, 0]),
            ("GCAATT", 1, [0, 0, 1, 0, 0, 0, 0, 1]),
        ]
        // number of edges adjcent to a kmer with base [leftA, leftT, leftG, leftC, rightG, rC, rA, rT];
        .iter()
        .map(map_kmer)
        .collect();
        assert_eq!(count_table, out1);

        let mut count_table = HashMap::new();
        let counter = KMerCounter {
            kmer_len: 3,
            min_cov: 1,
        };
        counter
            .count_kmers_from_fastq(&mut count_table, in1.as_bytes())
            .unwrap();

        let out1: CountTable = [
            ("ATT", 2, [2, 0, 0, 0, 1, 0, 0, 1]),
            ("TGC", 1, [0, 1, 0, 0, 0, 1, 0, 0]),
            ("GGC", 1, [0, 0, 0, 0, 0, 0, 1, 0]),
            ("TTG", 1, [1, 0, 0, 0, 0, 1, 0, 0]),
            ("TTT", 1, [1, 0, 0, 0, 0, 0, 0, 0]),
        ]
        .iter()
        .map(map_kmer)
        .collect();
        assert_eq!(count_table, out1);

        let mut count_table = HashMap::new(); // test kmer cov filter
        let counter = KMerCounter {
            kmer_len: 3,
            min_cov: 2,
        };
        counter
            .count_kmers_from_fastq(&mut count_table, in1.as_bytes())
            .unwrap();
        counter.kmer_coverage_filter(&mut count_table);

        let out1: CountTable = [("ATT", 2, [2, 0, 0, 0, 0, 0, 0, 0])]
            .iter()
            .map(map_kmer)
            .collect();
        assert_eq!(count_table, out1);

        // for (key, v) in count_table.iter() {
        //     println!("{}: {}", String::from_utf8_lossy(key), v);
        // }
    }

    #[test]
    fn partition_by_borders_test() {
        let borders = vec![20, 30, 40, 50, 60];
        assert_eq!(partition_by_borders(&10, &borders[..]), 0);
        assert_eq!(partition_by_borders(&20, &borders[..]), 0);
        assert_eq!(partition_by_borders(&21, &borders[..]), 1);
        assert_eq!(partition_by_borders(&30, &borders[..]), 1);
        assert_eq!(partition_by_borders(&31, &borders[..]), 2);
        assert_eq!(partition_by_borders(&35, &borders[..]), 2);
        assert_eq!(partition_by_borders(&60, &borders[..]), 4);
        assert_eq!(partition_by_borders(&66, &borders[..]), 5);
    }

    fn test_locality<F>(hash_func: F) -> (i32, f64)
    where
        F: Fn(&KMer, usize) -> usize,
    {
        let mut rng = rand::thread_rng();
        let read: Vec<u8> = (0..10000).map(|_| rng.gen::<u8>() % 4).collect();
        let mut kmers = vec![];
        let kmer_len = 31;
        let partition_num = 10;
        for i in 0..read.len() - kmer_len {
            kmers.push(read[i..i + kmer_len].to_vec());
        }
        let hashes: Vec<usize> = kmers
            .iter()
            .map(|kmer| hash_func(kmer, partition_num))
            .collect();
        let mut locality_counts = 0;
        for i in 0..hashes.len() - 1 {
            if hashes[i] != hashes[i + 1] {
                locality_counts += 1;
            }
        }

        let mut balance_count = vec![0usize; partition_num];
        hashes.iter().for_each(|p| balance_count[*p] += 1);
        let mean = (balance_count.iter().sum::<usize>() / partition_num) as f64;
        let mut v1 = 0f64;
        for i in hashes.iter() {
            v1 += (mean - *i as f64) * (mean - *i as f64);
        }
        let variance = v1.sqrt() / (partition_num - 1) as f64;
        (locality_counts, variance)
    }

    #[test]
    fn partition_by_hash_test() {
        let simplehash = |kmer: &KMer, partition_num| {
            let mut hasher = std::collections::hash_map::DefaultHasher::new();
            kmer.hash(&mut hasher);
            (hasher.finish() % partition_num as u64) as usize
        };
        let (locality_counts1, _var1) = test_locality(simplehash);
        // println!("count by hash: {} {}", locality_counts1, var1);
        let (locality_counts2, _var2) = test_locality(partition_by_hash);
        // println!("count by locality hash: {} {}", locality_counts2, var2);
        assert!(locality_counts2 * 10 < locality_counts1);
        // benchmark code
        // for i in 5..30 {
        //     let (c, v) = test_locality(|a, b| _partition_by_hash(a, b, i));
        //     println!("{}, {}, {}", i, c, v);
        // }
    }

    #[test]
    fn get_reverse_complement_test() {
        let k = new_kmer("ATGC");
        assert_eq!(new_kmer("GCAT"), get_reverse_complement(&k));
    }

    #[test]
    fn adj_one_side_degree_test() {
        let adj = [4u32, 0, 0, 0, 0, 3, 0, 6];
        assert_eq!(adj_one_side_degree(&adj, Left), 1);
        assert_eq!(adj_one_side_degree(&adj, Right), 2);
        assert_eq!(adj_degree(&adj), 3);
        let adj = [0u32, 0, 0, 0, 0, 0, 0, 0];
        assert_eq!(adj_one_side_degree(&adj, Left), 0);
        assert_eq!(adj_one_side_degree(&adj, Right), 0);
    }

    #[test]
    fn kmer_side_neighbor_test() {
        let kmer = new_kmer("ATCGATCG");
        let data = KMerData {
            adj: [1, 0, 1, 0, 0, 1, 0, 1],
            count: 2,
        };
        let l = kmer_side_neighbor(&kmer, &data, Left);
        let l: Vec<KMer> = l.iter().map(|(k, _)| k.clone()).collect();
        let r = kmer_side_neighbor(&kmer, &data, Right);
        let r: Vec<KMer> = r.iter().map(|(k, _)| k.clone()).collect();
        assert_eq!(l, vec![new_kmer("AATCGATC"), new_kmer("GATCGATC")]);
        assert_eq!(r, vec![new_kmer("TCGATCGC"), new_kmer("TCGATCGT")]);
    }

    #[test]
    fn local_linear_path_test() {
        let test = |in_seq: &str, kmer_len| {
            let counter = KMerCounter {
                kmer_len,
                min_cov: 1,
            };
            let mut count_table: CountTable = HashMap::new();
            counter
                .count_kmers_from_fastq(&mut count_table, in_seq.as_bytes())
                .unwrap();
            for (_, v) in count_table.iter() {
                assert_eq!(v.count, 1);
            }
            let table = local_linear_path(count_table, counter.kmer_len);
            let seq = get_canonical(in_seq.split('\n').nth(1).unwrap().as_bytes().to_vec()).0;
            let left_most_kmer = get_canonical(seq[..counter.kmer_len].to_vec()).0;
            let right_most_kmer =
                get_canonical(seq[seq.len() - counter.kmer_len..seq.len()].to_vec()).0;
            let ref linear_path = table[&left_most_kmer].connect.borrow_linear().ptr;
            let ref linear_path_r = table[&right_most_kmer].connect.borrow_linear().ptr;
            assert_eq!(linear_path.path, linear_path_r.path);
            assert_eq!(linear_path.path, seq);
            assert_eq!(
                linear_path.coverage as usize,
                seq.len() - counter.kmer_len + 1
            );
        };
        let in_seq1 = "@
ATTAAACACGGTCAGATTCG
+
AAAAAEEEEEEEEEEEEEEE
";

        let in_seq2 = "@
ATATTGCACAGGTGGCAAACTCCAACTGTTTCTTACGGTTTTATCGCAACCGGCAAATCCAGATTTTTCAGTATTTACACAAAGGGAGAGGGATTCTCTTTGTTAAAAAGTGACCATACCCGATGTCGTGCCCGGGTCAGCGCCACGTACA
+
AAAAAEEE/AEAE66EEEE//E/E/EAAEEE<//AEE/EEA//<AE//E<E//EE/EEEEEE//</AA6EEEEEEA//E/AEEEA/<AE/EEEAE66</EE/E/E<E/EAEA//6/EE</E/<6////E<E66/EAE///<<AEE<A<A/A
";
        let in_seq3 = "@
ATTAAACACGGTCAGATTCGGATCGGTTTTATGCAGCGCATCCAGCTCCGGCGTGGTTTCACGCGGATAACCGTACAGACTCATGCGGCCACGCTGGGTCGACTCGCCAATCACCAGCACTAAGGTGCGCGGTTCGTTACCCGATTCAT
+
AAAAAEEEEEEEEEEEEEEEEEEEEEEAEEEEEEEEEEEEEEEEEEEEAEAEEEEEEEEEEEEAEEEEEEEEA<EEEEEEEEEEEEAEEEEAEEEAAEEEAEEEEAEEEEEAEE6<AAEEAAEA<EEEEEAAEEEE<AEA<<<<<<AEA
";
        test(in_seq1, 11);
        test(in_seq2, 31);
        test(in_seq3, 31);
    }
    #[test]
    fn get_ended_kmer_test() {
        let path = LinearPath {
            path: "ABCDEFGHIJKLMN".as_bytes().to_vec(),
            coverage: 1,
        };
        assert_eq!(kmer_str(&path.get_ended_kmer(Left, 5)), "ABCDE");
        assert_eq!(kmer_str(&path.get_ended_kmer(Right, 5)), "JKLMN");
    }

    #[test]
    fn tip_travel_end_to_next_hop() {
        let seq = "TTAAACACGGTCAGATTCGGATCGGTTTTATGCAGCGCATCCAGCTCCGG";
        let mut kmer_table: KMerExtendTable = HashMap::new();
        let linear_path = LinearPath {
            path: seq.as_bytes().to_vec(),
            coverage: 1,
        };
        let boxed_path = Arc::new(linear_path);
        let kmer_left = "TTAAACACGGT".as_bytes().to_vec();
        let kmer_right = "TCCAGCTCCGG".as_bytes().to_vec();
        let connect = Connect::SomeLinear(LinearConnect {
            ptr: boxed_path.clone(),
            direction: Left,
        });
        kmer_table.insert(
            kmer_left.clone(),
            KMerExtendData {
                kmer_data: KMerData {
                    adj: [0, 0, 1, 0, 0, 1, 0, 0], // left G, right C
                    count: 20,
                },
                connect,
            },
        );

        let connect = Connect::SomeLinear(LinearConnect {
            ptr: boxed_path.clone(),
            direction: Right,
        });
        kmer_table.insert(
            kmer_right.clone(),
            KMerExtendData {
                kmer_data: KMerData {
                    adj: [1, 0, 0, 0, 0, 0, 1, 0], // left A right A
                    count: 20,
                },
                connect,
            },
        );
        let kmer_table = Arc::new(kmer_table);
        // from left to right
        let end_kmer = ContigGenerator::tip_travel_to_end_kmer(kmer_table.clone(), &kmer_left);
        assert_eq!(end_kmer, kmer_right);
        let (r, changed) = ContigGenerator::tip_travel_end_to_next_hop(kmer_table.clone(), &end_kmer);
        assert_eq!(changed, true);
        assert_eq!(r, "TCCGGAGCTGG".as_bytes().to_vec());

        // from right to left
        let end_kmer = ContigGenerator::tip_travel_to_end_kmer(kmer_table.clone(), &kmer_right);
        assert_eq!(end_kmer, kmer_left);
        let (r, changed) = ContigGenerator::tip_travel_end_to_next_hop(kmer_table.clone(), &end_kmer);
        assert_eq!(changed, false);
        assert_eq!(r, "GTTAAACACGG".as_bytes().to_vec());
        
    }
}
