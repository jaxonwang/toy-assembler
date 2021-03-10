use fastq::{Parser, Record, RefRecord};
use std::collections::HashMap;
// use std::env;
// use std::fs;
use std::convert::TryFrom;
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
type CountTable = HashMap<KMer, CountNumber>;
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

mod mpiconst {
    pub const ROOT_RANK: mpi::topology::Rank = 0;
    pub const DONE_MARK: &[u8] = "Done".as_bytes();
}

const COUNT_NUMBER_SIZE: usize = std::mem::size_of::<CountNumber>();
type KCBytes = [u8; COUNT_NUMBER_SIZE];
union KC {
    count: CountNumber,
    bytes: KCBytes,
}

fn pack(kmer: &KMer, num: CountNumber) -> Vec<u8> {
    let mut buf = kmer.clone();
    let kc = KC { count: num };
    buf.append(&mut Vec::from(unsafe { kc.bytes }));
    return buf;
}

fn unpack(buf: Vec<u8>) -> (KMer, CountNumber) {
    let ha: [u8; 4] = <[u8; 4]>::try_from(&buf[buf.len() - COUNT_NUMBER_SIZE..]).unwrap();
    let kc = KC { bytes: ha };
    return (
        buf.iter()
            .take(buf.len() - COUNT_NUMBER_SIZE)
            .cloned()
            .collect(),
        unsafe { kc.count },
    );
}

fn update_count_table(count_table: &mut CountTable, kmer: KMer, count: CountNumber) {
    let update_count = match count_table.get(&kmer) {
        Some(c) => c + count,
        None => count,
    };
    count_table.insert(kmer, update_count);
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

fn get_canonical(kmer: KMer) -> KMer {
    use std::cmp::Ordering::*;
    let rev_comp = get_reverse_complement(&kmer);
    match kmer.cmp(&rev_comp) {
        Less => rev_comp,
        Greater => kmer,
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
        let read_slice = Record::seq(&record);

        // drop too short read
        if read_slice.len() < self.kmer_len {
            return true;
        }
        // drop read with unknown base
        for b in read_slice.iter() {
            // illumina unknown base
            if *b == 'N' as u8 || *b == '.' as u8 {
                return true;
            }
        }

        // TODO: further optimization
        for start_pos in 0..(read_slice.len() - self.kmer_len + 1) {
            let kmer = Vec::from(&read_slice[start_pos..start_pos + self.kmer_len]);

            update_count_table(count_table, get_canonical(kmer), 1);
        }
        return true;
    }

    pub fn distributed_count(&self, count_table: CountTable, world: &SystemCommunicator) {
        let myrank = world.rank();
        let p = world.size();
        let sample_num: usize = 10000 / p as usize; // total 10000 samples
        let root_process = world.process_at_rank(mpiconst::ROOT_RANK);

        let mut borders: Vec<KMer> = vec![];
        for _ in 0..p - 1 {
            borders.push(vec![0u8; self.kmer_len]); // here I know the KMEr is u8 vec
        }

        info!("start sampling.");

        if myrank == mpiconst::ROOT_RANK {
            let mut done = 0;
            let mut samples: Vec<KMer> = vec![];
            for (k, v) in count_table.iter().take(sample_num) {
                for _ in 0..*v {
                    samples.push(k.clone());
                }
            }

            while done != p - 1 {
                let (msg, s) = world.any_process().receive_vec::<u8>();
                debug!(
                    "from {}, msg: {} length:{}",
                    s.source_rank(),
                    String::from_utf8_lossy(&msg[..]),
                    msg.len()
                );
                if &msg[..] == mpiconst::DONE_MARK {
                    done += 1;
                    debug!("all samples from {} are sent", s.source_rank());
                    continue;
                }
                samples.push(msg);
            }
            debug!("all samples received");

            // calculate borders
            samples.sort_unstable_by(|a, b| a.cmp(b));
            let real_sample_num = samples.len();
            let mut partition_size = real_sample_num / p as usize;
            if partition_size * (p as usize) < real_sample_num {
                // adjust for a good size
                partition_size = partition_size + 1;
            }

            let mut pos = 0usize;
            let mut i = 0;
            loop {
                pos += partition_size;
                if pos >= real_sample_num {
                    break;
                }
                borders[i] = samples[pos].clone();
                i += 1;
            }
            assert!(
                borders.len() == p as usize - 1,
                format!("Bad border! {} {} ", borders.len(), p)
            );
        } else {
            for (k, v) in count_table.iter().take(sample_num) {
                for _ in 0..*v {
                    root_process.send(&k[..]);
                }
            }
            // send done
            root_process.send(mpiconst::DONE_MARK);
        }

        // bradcast the border
        debug!("broadcasting borders.");
        for k in borders.iter_mut() {
            root_process.broadcast_into(&mut k[..]);
        }
        debug!("my borders:{:?}", &borders);
        info!("sampling done. Start counting.");

        let mut new_table: CountTable = HashMap::new();
        let local_table = Arc::new(Mutex::new(HashMap::new())); // stores entries partitioned to local

        crossbeam::scope(|scope| {
            // dispatch kemrs by borders
            scope.spawn(|_| self.dispatcher(count_table, world, borders, local_table.clone()));

            let mut done = 0;
            while done != p {
                debug!("try_to receive");
                let (msg, s) = world.any_process().receive_vec::<u8>();
                if msg == mpiconst::DONE_MARK {
                    done += 1;
                    debug!("{} Done!", s.source_rank());
                    continue;
                }
                let (kmer, count) = unpack(msg);
                // debug!(
                //     "receive {} from {}",
                //     String::from_utf8_lossy(&kmer),
                //     s.source_rank()
                // );
                update_count_table(&mut new_table, kmer, count);
            }
        })
        .unwrap();
        //
        let guard = local_table.lock().unwrap();
        for (k, v) in guard.iter() {
            update_count_table(&mut new_table, k.clone(), *v);
        }
        info!("Count Done");
        for (k, v) in new_table.iter() {
            println!("{}: {}", String::from_utf8_lossy(k), v);
        }
    }

    fn dispatcher(
        &self,
        count_table: CountTable,
        world: &SystemCommunicator,
        borders: Vec<KMer>,
        local_table: Arc<Mutex<CountTable>>,
    ) {
        let myrank = world.rank();
        for (k, v) in count_table.iter() {
            let dst = partition_by_borders(k, &borders[..]);
            if dst == myrank as usize {
                // if local set localtable
                let mut guard = local_table.lock().unwrap();
                update_count_table(&mut *guard, k.clone(), *v);
                continue;
            }
            let dst_process = world.process_at_rank(dst as mpi::topology::Rank);
            let msg = pack(k, *v);
            // debug!(
            //     "send {} to {}",
            //     String::from_utf8_lossy(&msg),
            //     dst
            // );
            dst_process.send(&msg[..]);
        }
        for i in 0..world.size() {
            world.process_at_rank(i).send(mpiconst::DONE_MARK);
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

    let counter = KMerCounter { kmer_len: 4 };
    let in1 = "@
AAATTTCC
+
AAAAAEEE
";
    let f = in1.as_bytes();

    // use std::env;
    // use std::fs;
    // let fastq_path = match env::args().nth(1) {
    //     None => {
    //         println!("usage: a.out file");
    //         std::process::exit(1)
    //     }
    //     Some(p) => p,
    // };
    // let f = fs::File::open(fastq_path).unwrap();
    //
    let mut count_table = HashMap::new();
    counter.count_kmers_from_fastq(&mut count_table, f).unwrap();
    counter.distributed_count(count_table, &world);
}

#[cfg(test)]
mod tests {
    use super::*;

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
        let out1: CountTable = [("GGCAAT", 1), ("CAATTT", 1), ("GCAATT", 1)]
            .iter()
            .map(|x| (new_kmer(x.0), x.1))
            .collect();
        assert_eq!(count_table, out1);

        let mut count_table = HashMap::new();
        let counter = KMerCounter { kmer_len: 3 };
        counter
            .count_kmers_from_fastq(&mut count_table, in1.as_bytes())
            .unwrap();
        for (k, v) in count_table.iter() {
            println!("{} {}", String::from_utf8_lossy(k), v);
        }

        let out1: CountTable = [("ATT", 2), ("TGC", 1), ("GGC", 1), ("TTG", 1), ("TTT", 1)]
            .iter()
            .map(|x| (new_kmer(x.0), x.1))
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

    #[test]
    fn packing_test() {
        let kmer: KMer = "ABCDE".as_bytes().to_vec();
        let num = 123;
        let buf = pack(&kmer, num);
        let (kmer1, num1) = unpack(buf);
        assert_eq!(num, num1);
        assert_eq!(kmer1, kmer);
    }

    #[test]
    fn get_reverse_complement_test() {
        let k = new_kmer("ATGC");
        assert_eq!(new_kmer("GCAT"), get_reverse_complement(&k));
    }
}
