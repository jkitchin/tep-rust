//! The Tier 1 measuring apparatus: input pools, and a comparison that reports
//! a distribution instead of a verdict.
//!
//! Tier 1 is the bottom rung of the validation ladder. It proves that the
//! utility routines `TESUB1` through `TESUB4` agree with the original to near
//! machine precision over every composition and temperature the plant can
//! reach. Nothing above it means anything if it fails, since every unit
//! operation calls them.
//!
//! # Why this exists before the routines it measures
//!
//! Building the instrument first keeps the tolerance from being chosen to fit
//! whatever the port happens to produce. The pools, the metric and the report
//! format are fixed here, with no ported routine in the workspace to tune them
//! against, so the first number `TESUB1` produces is a measurement rather than
//! a negotiation.
//!
//! # Using it
//!
//! ```no_run
//! # #[cfg(feature = "oracle")] {
//! use tepsim_oracle::{Oracle, tier1::{Case, Comparison, Sweep}};
//!
//! let sweep = Sweep::SMOKE;
//! let mut fortran = Oracle::lock();
//! let mut comparison: Comparison<Case> = Comparison::new("TESUB1 ity=0");
//!
//! for case in sweep.cases() {
//!     let reference = fortran.tesub1(&case.z(), case.celsius, 0);
//!     let actual = reference; // the port goes here
//!     comparison.observe(case, actual, reference);
//! }
//!
//! println!("{comparison}");
//! comparison.assert_within(1e-13);
//! # }
//! ```
//!
//! The `println!` is not decoration. `CLAUDE.md` requires the measured numbers
//! in the log entry, and [`Comparison`]'s `Display` is the format they go in.
//!
//! # Not gated on the `oracle` feature
//!
//! Everything here is pure Rust: it generates inputs and compares numbers, and
//! never touches the Fortran itself. That means the harness can be tested, and
//! its generators proved correct, on a machine with no `gfortran` at all. Only
//! the differential test that feeds it the Fortran needs the toolchain.

mod compare;
mod sweeps;

pub use compare::{Comparison, EXACT_BUCKETS, relative_error, ulp_distance};
pub use sweeps::{
    BREAKPOINTS, Breakpoint, Case, Pool, Sampler, SimplexGrid, Sweep, TemperatureRange,
};
