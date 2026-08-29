//! Tennessee Eastman Process simulator.
//!
//! The high-level API: build a [`Scenario`], run a [`Simulation`], get a
//! [`Run`]. Everything below this is [`tepsim_core`] for the plant and
//! [`tepsim_control`] for the controllers, and a caller who wants to reach past
//! the facade can.
//!
//! # What this is a port of
//!
//! Downs and Vogel's 1993 Tennessee Eastman challenge problem, from the
//! original Fortran rather than from any later reimplementation. The port is
//! validated against that Fortran through a ten-tier ladder; the headline
//! result is that a complete 48-hour closed-loop run, 172,800 integrator steps,
//! is **bit-identical** in all 41 measurements and all 12 manipulated variables
//! when both are given the same `exp` and `pow`. See `book/src/deltas.md` for
//! every place the port knowingly differs, and why.
//!
//! # Example
//!
//! ```
//! use tepsim::{Scenario, Simulation};
//!
//! // Two hours of the fault-free plant, sampled every three minutes.
//! let run = Simulation::new(Scenario::baseline().with_hours(2.0)).run();
//! assert_eq!(run.samples.len(), 40);
//!
//! // With a disturbance, and its ground truth recorded alongside.
//! let run = Simulation::new(Scenario::fault(1).with_hours(2.0)).run();
//! assert!(run.samples[0].labels.faulted());
//! assert_eq!(run.samples[0].labels.faults().collect::<Vec<_>>(), vec![1]);
//! ```
//!
//! # Determinism
//!
//! A run is a pure function of its [`Scenario`]. No clock, no thread-local
//! state, no global. The same scenario gives bit-identical output on x86-64,
//! aarch64 and wasm32, which is what makes a recorded dataset reproducible from
//! its description rather than from a file.
//!
//! That is a claim, so it is measured. [`tier9`] holds a table of fixed
//! scenarios and the digest each is committed to produce, and every platform
//! checks itself against the same committed constants. `cargo xtask tier9`
//! runs the table natively and again on `wasm32-unknown-unknown`.

#![no_std]
#![forbid(unsafe_code)]

extern crate alloc;

pub mod integrator;
pub mod recorder;
pub mod run;
pub mod scenario;
pub mod sim;
pub mod text;
pub mod tier9;

pub use integrator::{Integrator, Step};
pub use recorder::{Columnar, Csv, CsvString, Decimating, Recorder, Ring, Selecting};
pub use run::{CHANNELS, Labels, MANIPULATED, MEASUREMENTS, Outcome, Run, Sample, channel_names};
pub use scenario::{DISTURBANCES, SCENARIO_VERSION, Scenario};
pub use sim::{Simulation, forced_disturbance_step};
pub use tepsim_scenario::{Action, Digest, Event, Invalid, Schedule};
pub use text::TextError;

pub use tepsim_control;
pub use tepsim_core;
pub use tepsim_scenario;
