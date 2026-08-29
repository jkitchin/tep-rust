//! The Rust snippet in `README.md`, so it cannot rot.
use tepsim::{Scenario, Simulation};

fn main() {
    let run = Simulation::new(Scenario::fault(4).with_hours(24.0)).run();
    println!("{:?} after {} samples", run.outcome, run.samples.len());
}
