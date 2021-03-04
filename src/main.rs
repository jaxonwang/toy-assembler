extern crate mpi;

// use mpi::request::WaitGuard;
use mpi::traits::*;
use std::thread;
use std::time;


fn main() {

    let universe = mpi::initialize().unwrap();
    let world = universe.world();
    // let size = world.size();
    let rank = world.rank();
    println!("Hello, world! {}", rank);
    thread::sleep(time::Duration::from_secs(10));
}
