use fastq::{Parser, Record, RefRecord};
use std::collections::HashMap;
use std::env;
use std::fs;
use std::io;

extern crate fastq;

type KMer = Vec<u8>;
type CountTable = HashMap<KMer, u32>;

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

use mpi::topology::{Process, SystemCommunicator};
use mpi::topology;

fn partition_by_borders<O: Ord>(kmer: &O, borders: &[O]) -> usize{
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
            Less => j = middle
        };
    }
    if i == middle + 1 {
        return i
    }
    return  middle;
}


// fn dispatcher(count_table: CountTable, world: &SystemCommunicator, borders: &[KMer]) {}
//
// pub fn distributed_count(count_table: CountTable, world: &SystemCommunicator) {}
//
pub fn main() {
    let kmer_len= 30usize;

    let fastq_path = match env::args().nth(1) {
        None => {
            println!("usage: a.out file");
            std::process::exit(1)
        }
        Some(p) => p,
    };
    let f = fs::File::open(fastq_path).unwrap();

    let mut count_table = HashMap::new();
    count_kmers_from_fastq(&mut count_table, kmer_len, f).unwrap();

    for (key, v) in count_table.iter() {
        println!("{}: {}", String::from_utf8_lossy(key), v);
    }
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
    fn partition_by_borders_test(){
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

