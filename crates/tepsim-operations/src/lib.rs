//! Deciding what to do about the plant, rather than simulating it.
//!
//! Tiers 1 through 10 ask whether this port *is* the Fortran. This crate asks a
//! different kind of question: given drifting feed quality and drifting prices,
//! what does a run cost, what must the instruments have meant, and what should
//! the setpoints be. See `PLAN.org`, "Operations".
//!
//! # Why this is not in `tepsim-core`
//!
//! The core is `no_std`, `forbid(unsafe_code)` and bound by the determinism
//! invariant, and its bit-exactness against Fortran is the project's central
//! claim. An iterative optimiser has no business inside that. Nothing here
//! changes a single number the simulator produces; this crate only reads its
//! output and proposes [`tepsim::Action::Setpoint`] events, which the schedule
//! already knows how to apply.
//!
//! # Why it carries its own linear algebra
//!
//! `tepsim-stats` has the routines, and this crate may not use them:
//! `tepsim-stats` is development-only and `cargo xtask ci` asserts that no
//! shipped crate so much as names it. That is the right rule and this is the
//! price of it.
//!
//! # What is validated here and what is not
//!
//! The *structure* is validated: unit conversions are derived from
//! `teprob.f`'s own definitions of each channel and asserted, and every term
//! has a test for the direction it moves the answer.
//!
//! The *prices* are not, and cannot be. See [`economics::Prices`].

#![forbid(unsafe_code)]

pub mod economics;

pub use economics::{CostRate, Prices};
