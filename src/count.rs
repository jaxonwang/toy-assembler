use fastq::{Parser, Record, RefRecord};
use std::collections::HashMap;
use std::convert::TryFrom;
use std::hash::Hash;
use std::hash::Hasher;
use std::io;
use std::sync::{Arc, Mutex};

extern crate chrono;
extern crate crossbeam;
extern crate fastq;
extern crate fern;
extern crate log;
extern crate mpi;

use log::{debug, info};
use mpi::topology::SystemCommunicator;
use mpi::traits::*;
use mpi::Threading;

type CountNumber = u32;
type KMer = Vec<u8>;

#[derive(Debug, Copy, Clone, Eq)]
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

const BASE_A: u8 = 'A' as u8;
const BASE_T: u8 = 'T' as u8;
const BASE_C: u8 = 'C' as u8;
const BASE_G: u8 = 'G' as u8;

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
    return middle;
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

const KMER_DATA_SIZE: usize = std::mem::size_of::<KMerData>();
type KDBytes = [u8; KMER_DATA_SIZE];
union KD {
    kmer_data: KMerData,
    bytes: KDBytes,
}

fn pack(buf: &mut [u8], kmer: &KMer, kmer_data: &KMerData) -> usize {
    buf[..kmer.len()].copy_from_slice(&kmer[..]);
    let kc = KD {
        kmer_data: kmer_data.clone(),
    };
    buf[kmer.len()..kmer.len()+KMER_DATA_SIZE].copy_from_slice(&unsafe { kc.bytes }[..]);
    return KMER_DATA_SIZE + kmer.len();
}

fn unpack(buf: &[u8]) -> (KMer, KMerData) {
    let kmer = buf[..buf.len() - KMER_DATA_SIZE].to_vec();
    let bytes: [u8; KMER_DATA_SIZE] =
        <[u8; KMER_DATA_SIZE]>::try_from(&buf[buf.len() - KMER_DATA_SIZE..]).unwrap();
    let kc = KD { bytes };
    return (kmer, unsafe { kc.kmer_data });
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

fn get_reverse_complement(kmer: &KMer) -> KMer {
    kmer.iter()
        .rev()
        .map(|base| match *base {
            BASE_A => BASE_T,
            BASE_T => BASE_A,
            BASE_C => BASE_G,
            BASE_G => BASE_C,
            other => panic!(format!("unknown base: {}", other)),
        })
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

fn adj_new_edge_with_base(adj: &mut [u32; 8], is_left: bool, base: u8) {
    adj[adj_position(is_left, base)] += 1;
}

fn get_canonical(kmer: KMer) -> (KMer, bool) {
    // bool is true if there is a change
    use std::cmp::Ordering::*;
    let rev_comp = get_reverse_complement(&kmer);
    match kmer.cmp(&rev_comp) {
        Less => (rev_comp, true),
        Greater => (kmer, false),
        Equal => panic!("KMer is plalindrome"),
    }
}

struct KMerCounter {
    kmer_len: usize,
}

impl KMerCounter {
    pub fn count_kmers_from_fastq(
        &self,
        count_table: &mut CountTable,
        r: impl io::Read,
    ) -> io::Result<()> {
        let parser = Parser::new(r);

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
            if *b == 'N' as u8 || *b == '.' as u8 {
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
        return true;
    }

    pub fn distributed_count(&self, count_table: CountTable, world: &SystemCommunicator) {
        let p = world.size();

        let mut new_table: CountTable = HashMap::new();
        let local_table = Arc::new(Mutex::new(HashMap::new())); // stores entries partitioned to local

        crossbeam::scope(|scope| {
            // dispatch kemrs by borders
            scope.spawn(|_| self.dispatcher(count_table, world, local_table.clone()));

            let mut done = 0;
            while done != p - 1{
                debug!("try_to receive");
                let (msg, s) = world.any_process().receive_vec::<u8>();
                if msg == mpiconst::DONE_MARK {
                    done += 1;
                    debug!("{} Done!", s.source_rank());
                    continue;
                }
                let item_size = self.kmer_len + KMER_DATA_SIZE;
                let mut pos = 0usize;
                // unpack the msg one by one
                while pos < msg.len(){
                    let (kmer, kmerdata) = unpack(&msg[pos..pos+item_size]);
                    update_count_table(&mut new_table, kmer, kmerdata);
                    pos += item_size;
                }
            }
        })
        .unwrap();
        // update from local
        let guard = local_table.lock().unwrap();
        for (k, v) in guard.iter() {
            update_count_table(&mut new_table, k.clone(), *v);
        }
        info!("Count Done: kmers number {}", new_table.len());
        // for (k, v) in new_table.iter() {
        //     println!("{}: {:?}", String::from_utf8_lossy(k), v);
        // }
    }

    fn dispatcher(
        &self,
        count_table: CountTable,
        world: &SystemCommunicator,
        local_table: Arc<Mutex<CountTable>>,
    ) {
        let myrank = world.rank();
        let buffer_size = 140;
        let item_size = self.kmer_len + KMER_DATA_SIZE;
        let item_num = buffer_size / item_size;
        let mut buffers:Vec<(Vec<u8>, usize)> = vec![];
        for _ in 0..world.size() {
           buffers.push((vec![0u8; item_size * item_num], 0usize));
        }
        for (k, v) in count_table.iter() {
            let dst = partition_by_hash(k, world.size() as usize);
            if dst == myrank as usize {
                // if local set localtable
                let mut guard = local_table.lock().unwrap();
                update_count_table(&mut *guard, k.clone(), *v);
                continue;
            }
            let pos = buffers[dst].1;
            if pos == buffers[dst].0.len(){
                let dst_process = world.process_at_rank(dst as mpi::topology::Rank);
                dst_process.send(&buffers[dst].0[..]); // send full buf
                // empty buffer
            }else{
                pack(&mut buffers[dst].0[pos..pos+item_size], k, v);
                buffers[dst].1 += item_size;
            }

            // debug!(
            //     "send {} to {}",
            //     String::from_utf8_lossy(&msg),
            //     dst
            // );
        }
        for dst in 0..world.size(){
            if dst != myrank {
                world.process_at_rank(dst).send(&buffers[dst as usize].0[0..buffers[dst as usize].1]);
                world.process_at_rank(dst).send(mpiconst::DONE_MARK);
            }
        } 
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
    setup_logger(process_id.clone()).unwrap();

    let counter = KMerCounter { kmer_len: 31 };
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
    counter.distributed_count(count_table, &world);
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

        let counter = KMerCounter { kmer_len: 6 };
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
        .map(|x| {
            (
                new_kmer(x.0),
                KMerData {
                    count: x.1,
                    adj: x.2,
                },
            )
        })
        .collect();
        assert_eq!(count_table, out1);

        let mut count_table = HashMap::new();
        let counter = KMerCounter { kmer_len: 3 };
        counter
            .count_kmers_from_fastq(&mut count_table, in1.as_bytes())
            .unwrap();
        for (k, v) in count_table.iter() {
            println!("{} {:?}", String::from_utf8_lossy(k), v);
        }

        let out1: CountTable = [
            ("ATT", 2, [2, 0, 0, 0, 1, 0, 0, 1]),
            ("TGC", 1, [0, 1, 0, 0, 0, 1, 0, 0]),
            ("GGC", 1, [0, 0, 0, 0, 0, 0, 1, 0]),
            ("TTG", 1, [1, 0, 0, 0, 0, 1, 0, 0]),
            ("TTT", 1, [1, 0, 0, 0, 0, 0, 0, 0]),
        ]
        .iter()
        .map(|x| {
            (
                new_kmer(x.0),
                KMerData {
                    count: x.1,
                    adj: x.2,
                },
            )
        })
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
        let hashes: Vec<usize> = kmers.iter().map(|kmer| hash_func(kmer, partition_num)).collect();
        let mut locality_counts = 0;
        for i in 0..hashes.len() - 1 {
            if hashes[i] != hashes[i + 1] {
                locality_counts += 1;
            }
        }

        let mut balance_count = vec![0usize; partition_num];
        hashes.iter().map(|p|balance_count[*p]+=1).collect::<()>();
        let mean = (balance_count.iter().sum::<usize>() / partition_num) as f64; 
        let mut v1 =0f64;
        for i in hashes.iter(){
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
        let (locality_counts1, var1) = test_locality(simplehash);
        println!("count by hash: {} {}", locality_counts1, var1);
        let (locality_counts2, var2) = test_locality(partition_by_hash);
        println!("count by locality hash: {} {}", locality_counts2, var2);
        assert!(locality_counts2 * 10 < locality_counts1);
        // benchmark code
        // for i in 5..30 {
        //     let (c, v) = test_locality(|a, b| _partition_by_hash(a, b, i));
        //     println!("{}, {}, {}", i, c, v);
        // }
    }

    #[test]
    fn packing_test() {
        let kmer: KMer = "ABCDE".as_bytes().to_vec();
        let data = KMerData {
            count: 1,
            adj: [1, 2, 3, 4, 5, 6, 7, 8],
        };
        let mut buf = vec![0u8; kmer.len()+KMER_DATA_SIZE];
        pack(&mut buf[..], &kmer, &data);
        let (kmer1, data1) = unpack(&buf[..]);
        assert_eq!(data, data1);
        assert_eq!(kmer1, kmer);
    }

    #[test]
    fn get_reverse_complement_test() {
        let k = new_kmer("ATGC");
        assert_eq!(new_kmer("GCAT"), get_reverse_complement(&k));
    }
}
