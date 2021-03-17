use std::io::{Error, ErrorKind, Result};

use serde::{Deserialize, Serialize};
use serde::de::DeserializeOwned;

extern crate bincode;

pub trait Packed {
    fn pack(&self, buf: &mut Vec<u8>);
    fn unpack(buf: &[u8]) -> (Result<Self>, usize)
    where
        Self: Sized;
}

const USIZE_LEN: usize = std::mem::size_of::<usize>();
union UsizeBytes {
    n: usize,
    bytes: [u8; USIZE_LEN],
}

impl<T: Serialize + DeserializeOwned> Packed for T {
    fn pack(&self, buf: &mut Vec<u8>) {
        // let type_size:usize = std::mem::size_of::<T>();
        let s = bincode::serialize(&self).unwrap();
        let len = UsizeBytes { n: s.len() };

        buf.extend_from_slice(&unsafe { len.bytes }[..]);
        buf.extend_from_slice(&s[..]);
    }
    fn unpack(buf: & [u8]) -> (Result<T>, usize) {
        match buf.len() {
            0 => (Err(Error::new(ErrorKind::InvalidData, "Reaching End")), 0),
            l if l < USIZE_LEN => panic!("bad msg"),
            _ => {
                let mut bytes = [0u8; USIZE_LEN];
                bytes.copy_from_slice(&buf[..USIZE_LEN]);
                let len = unsafe { UsizeBytes { bytes }.n };
                (Ok(bincode::deserialize(&buf[USIZE_LEN..]).unwrap()), len + USIZE_LEN)
            }
        }
    }
}

#[cfg(test)]
mod tests {

    use super::*;

    #[derive(Serialize, Deserialize, Debug)]
    struct A {
        a: i32,
        b: f64,
        c: Vec<u8>,
    }
    impl PartialEq for A {
        fn eq(&self, other: &Self) -> bool {
            self.a == other.a && self.b == other.b && self.c == other.c
        }
    }

    #[test]
    fn pack_test() {
        let mut buf: Vec<u8> = vec![];
        let a = A {
            a: 123,
            b: 3.14159,
            c: vec![0u8, 1, 2, 3, 4, 5, 6, 7, 8, 9],
        };
        a.pack(&mut buf);
        if let (Ok(a1), _) = A::unpack(&buf[..]) {
            assert_eq!(a, a1);
        } else {
            assert!(false);
        }
    }
}
