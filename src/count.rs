use fastq::{Parser, Record, RefRecord};
use serde::{Deserialize, Serialize};
use std::cell::{Ref, RefCell};
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

fn adj_position(is_left: bool, base: u8) -> usize {
    if is_left {
        match base {
            BASE_A => 0,
            BASE_T => 1,
            BASE_G => 2,
            BASE_C => 3,
            _ => panic!("bad base"),
        }
    } else {
        match base {
            BASE_G => 4,
            BASE_C => 5,
            BASE_A => 6,
            BASE_T => 7,
            _ => panic!("bad base"),
        }
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

fn adj_new_edge_with_base(adj: &mut [u32; 8], is_left: bool, base: u8) {
    adj[adj_position(is_left, base)] += 1;
}

fn kmer_side_neighbor(kmer: &KMer, data: &KMerData, is_left: bool) -> Vec<KMer> {
    // return the out edge bases
    let range = match is_left {
        true => 0..4,
        false => 4..8,
    };
    let mut neighbors: Vec<KMer> = vec![];
    for i in range {
        if data.adj[i] > 0 {
            let mut neighbor_kmer: KMer = vec![0u8; kmer.len()];
            if is_left {
                neighbor_kmer[1..].copy_from_slice(&kmer[..kmer.len() - 1]);
                neighbor_kmer[0] = adj_base(i);
            } else {
                neighbor_kmer[..kmer.len() - 1].copy_from_slice(&kmer[1..]);
                neighbor_kmer[kmer.len() - 1] = adj_base(i);
            }
            neighbors.push(neighbor_kmer);
        }
    }
    neighbors
}

fn adj_one_side_degree(adj: &[u32; 8], is_left: bool) -> usize {
    let range = match is_left {
        true => 0..4,
        false => 4..8,
    };
    adj[range]
        .iter()
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
                adj_new_edge_with_base(&mut adj, true, read[start_pos - 1]);
            }
            if start_pos + self.kmer_len < read.len() {
                //edge to right
                adj_new_edge_with_base(&mut adj, false, read[start_pos + self.kmer_len]);
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
        let mut low_cov_kmers: Vec<KMer> = vec![];
        for (k, data) in count_table.iter() {
            if (data.count as usize) < self.min_cov {
                low_cov_kmers.push(k.clone());
            }
        }
        for k in low_cov_kmers {
            count_table.remove(&k);
        }
    }
}

struct LinearPath {
    // two sides
    path: Vec<u8>,
    coverage: u32,
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
    start_kmer: KMer,
    end_kmer: KMer,
    coverage: u32,
}

impl CompactNode {
    fn new(p: &[u8], kmer_len: usize, cov: u32) -> Self {
        assert!(
            kmer_len >= p.len(),
            format!("vec too short: {} need: {}", p.len(), kmer_len)
        );
        CompactNode {
            start_kmer: p[..kmer_len].to_vec(),
            end_kmer: p[p.len() - kmer_len..p.len()].to_vec(),
            coverage: cov,
        }
    }
}

#[derive(Serialize, Deserialize)]
struct SubGraph {
    //CSC format
    nodes: Vec<CompactNode>, //sorted
    adjs: Vec<Vec<[usize;8]>>,   // adj[0] -> list of index node 0 connect to
}

impl SubGraph {
    fn new_single_kmer(k: &KMer, kmer_len: usize, cov: u32) -> Self {
        SubGraph {
            nodes: vec![CompactNode::new(&k[..], kmer_len, cov)],
            adjs: vec![vec![]],
        }
    }
}

enum Connect {
    SomeGraph(Rc<RefCell<SubGraph>>),
    SomeLinear(Rc<RefCell<LinearPath>>),
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
    fn new_some_graph(g: SubGraph) -> Self {
        Self::SomeGraph(Rc::new(RefCell::new(g)))
    }
    fn new_some_linear(l: LinearPath) -> Self {
        Self::SomeLinear(Rc::new(RefCell::new(l)))
    }
    #[allow(dead_code)]
    fn unwrap_graph(self) -> Rc<RefCell<SubGraph>> {
        match self {
            Self::SomeGraph(g) => g,
            _ => panic!("Connect is not graph"),
        }
    }
    #[allow(dead_code)]
    fn unwrap_linear(self) -> Rc<RefCell<LinearPath>> {
        match self {
            Self::SomeLinear(l) => l,
            _ => panic!("Connect is not linear"),
        }
    }
    fn borrow_linear(&self) -> Ref<LinearPath> {
        match self {
            Self::SomeLinear(l) => l.borrow(),
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
        let is_linear = |x| {
            adj_one_side_degree(&kmers_from_counter[x].adj, false)
                + adj_one_side_degree(&kmers_from_counter[x].adj, true)
                <= 2
        };
        if visited.contains(k) {
            continue;
        }

        let mut travel = |is_left| -> Vec<(u8, u32)> {
            let mut side_path: Vec<(u8, u32)> = vec![]; // base, count
            let mut current: &KMer = k;
            let mut direction = is_left;
            // extend to one direction
            while !visited.contains(current) && is_linear(current) {
                visited.insert(current.clone()); // now k is linear && not visited. mark as visited.
                let base;
                if is_left {
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
                let neighbor = neighbor.pop().unwrap();
                let (ref extend, changed) = get_canonical(neighbor);
                if changed {
                    // change to it reverse complement, to go the same direction, need change
                    direction = !direction;
                }
                if kmers_from_counter.contains_key(extend) {
                    // is local
                    current = kmers_from_counter.get_key_value(extend).unwrap().0;
                }
            }
            if !is_left {
                visited.remove(k); // remove for left travel enabled
            }
            side_path
        };

        let rpath = travel(false);
        let lpath = travel(true); //extend left

        if lpath.len() == 0 {
            // nothing to do for this kmer, store it as SubGraph
            let sub = SubGraph::new_single_kmer(k, kmer_len, v.count);
            extend_table.insert(
                k.clone(),
                KMerExtendData {
                    kmer_data: v.clone(),
                    connect: Connect::new_some_graph(sub),
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
        let cnc = Connect::new_some_linear(LinearPath { path, coverage });
        let l_kmer_data = kmers_from_counter[&left_end_kmer].clone();
        let r_kmer_data = kmers_from_counter[&right_end_kmer].clone();
        extend_table.insert(
            // insert look up for left
            left_end_kmer,
            KMerExtendData {
                kmer_data: l_kmer_data,
                connect: cnc.clone(),
            },
        );
        extend_table.insert(
            // for right
            right_end_kmer,
            KMerExtendData {
                kmer_data: r_kmer_data,
                connect: cnc.clone(),
            },
        );
    }
    info!("linear path generation done.");
    extend_table
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

        let out1: CountTable = [("ATT", 2, [2, 0, 0, 0, 1, 0, 0, 1])]
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
        assert_eq!(adj_one_side_degree(&adj, true), 1);
        assert_eq!(adj_one_side_degree(&adj, false), 2);
        let adj = [0u32, 0, 0, 0, 0, 0, 0, 0];
        assert_eq!(adj_one_side_degree(&adj, true), 0);
        assert_eq!(adj_one_side_degree(&adj, false), 0);
    }

    #[test]
    fn kmer_side_neighbor_test() {
        let kmer = new_kmer("ATCGATCG");
        let data = KMerData {
            adj: [1, 0, 1, 0, 0, 1, 0, 1],
            count: 2,
        };
        let l = kmer_side_neighbor(&kmer, &data, true);
        let r = kmer_side_neighbor(&kmer, &data, false);
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
            let linear_path = table[&left_most_kmer].connect.borrow_linear();
            let linear_path_r = table[&right_most_kmer].connect.borrow_linear();
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
}
