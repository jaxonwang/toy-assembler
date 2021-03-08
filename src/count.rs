use fastq::{Parser, Record, RefRecord};
use std::collections::HashMap;
// use std::env;
// use std::fs;
use std::io;
use std::sync::Arc;

extern crate fastq;
extern crate mpi;
extern crate log;
extern crate fern;
extern crate chrono;


use mpi::topology::SystemCommunicator;
use mpi::traits::*;
use log::{debug, info};

type CountNumber = u32;
type KMer = Vec<u8>;
type CountTable = HashMap<KMer, CountNumber>;

fn all_reads(count_table: &mut CountTable, kmer_len: usize, record: RefRecord) -> bool {
    let read_slice = Record::seq(&record);

    // drop too short read
    if read_slice.len() < kmer_len {
        return true;
    }

    for start_pos in 0..(read_slice.len() - kmer_len + 1) {
        let kmer = Vec::from(&read_slice[start_pos..start_pos + kmer_len]);

        let count = match count_table.get(&kmer) {
            Some(c) => c + 1,
            None => 1,
        };
        count_table.insert(kmer, count);
    }
    return true;
}

pub fn count_kmers_from_fastq(
    count_table: &mut CountTable,
    kmer_len: usize,
    r: impl io::Read,
) -> io::Result<()> {
    let parser = Parser::new(r);

    parser
        .each(|record| all_reads(count_table, kmer_len, record))
        .map(|_| ()) // I dont' care the bool in the result
}

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
}

// fn dispatcher(count_table: CountTable, world: &SystemCommunicator, borders: &[KMer]) {
//
// }
//

struct KmerCounter {
    kmer_len: usize,
}

impl KmerCounter {
    pub fn distributed_count(&self, count_table: CountTable, world: &SystemCommunicator) {

        let myrank = world.rank();
        let p = world.size();
        let sample_num: usize = 10000 / p as usize; // total 10000 samples
        let done_mark = "done".as_bytes();
        let root_process = world.process_at_rank(mpiconst::ROOT_RANK);

        let mut borders:Vec<KMer> = vec![];
        for _ in 0..p-1{
            borders.push(vec![0u8;self.kmer_len]); // here I know the KMEr is u8 vec
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
                if &msg[..] == done_mark {
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
            root_process.send(done_mark);
        }

        debug!("broadcasting borders.");
        for k in borders.iter_mut(){
            root_process.broadcast_into(&mut k[..]);
        }
        debug!("my borders:{:?}", &borders);
        info!("sampling done");
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


fn setup_logger(process_id:Arc<String>) -> Result<(), fern::InitError> {
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
        .level(log::LevelFilter::Debug)
        .chain(std::io::stdout())
        .apply()?;
    Ok(())
}

fn main() {

    let universe = mpi::initialize().unwrap();
    let world = universe.world();

    let process_id = Arc::new(format!("P{}", world.rank()));
    setup_logger(process_id.clone());

    let counter = KmerCounter { kmer_len: 3 };
    let in1 = "@
AAATTTCC
+
AAAAAEEE
";

    let mut count_table = HashMap::new();
    count_kmers_from_fastq(&mut count_table, 3, in1.as_bytes()).unwrap();
    counter.distributed_count(count_table, &world);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn count_kmers_from_fastq_test() {
        let in1 = "@
AAATTTCC
+
AAAAAEEE
";

        let mut count_table = HashMap::new();
        count_kmers_from_fastq(&mut count_table, 6, in1.as_bytes()).unwrap();
        let out1: CountTable = [("AAATTT", 1), ("AATTTC", 1), ("ATTTCC", 1)]
            .iter()
            .map(|x| (x.0.as_bytes().to_vec(), x.1))
            .collect();
        assert_eq!(count_table, out1);

        let mut count_table = HashMap::new();
        count_kmers_from_fastq(&mut count_table, 3, in1.as_bytes()).unwrap();

        let out1: CountTable = [
            ("AAT", 1),
            ("TTC", 1),
            ("AAA", 1),
            ("TCC", 1),
            ("TTT", 1),
            ("ATT", 1),
        ]
        .iter()
        .map(|x| (x.0.as_bytes().to_vec(), x.1))
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
}
