//! High-level Tennessee Eastman Process simulator.
//!
//! The facade over [`tepsim_core`], [`tepsim_control`] and [`tepsim_scenario`]:
//! the run loop, integrators, recorders, and dataset writers.
//!
//! # Status
//!
//! Skeleton. Lands in phase 6; see `BACKLOG.org`.

#![forbid(unsafe_code)]

pub use tepsim_control;
pub use tepsim_core;
pub use tepsim_scenario;
