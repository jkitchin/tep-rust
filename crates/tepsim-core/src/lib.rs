//! Tennessee Eastman Process plant model.
//!
//! This crate holds the process model itself: state, thermodynamics, reaction
//! kinetics, the stream network, and the measurement layer. It performs no I/O,
//! owns no global state, and spawns no threads, which is what lets the same code
//! serve native binaries, Python wheels, and a WebAssembly bundle.
//!
//! # Invariants
//!
//! These are not stylistic preferences. Validation tiers 9 and 10 depend on them.
//!
//! - **`no_std`.** Only [`alloc`] is used, and sparingly.
//! - **No `unsafe`.** Enforced by `#![forbid(unsafe_code)]`.
//! - **Determinism.** All arithmetic is `f64`. No `f32`, no SIMD, no reordered
//!   reductions, no parallelism. Transcendental functions come from a vendored
//!   pure-Rust `libm` so that x86-64, aarch64 and wasm32 agree bit for bit.
//! - **Provenance.** Every function ported from the original Fortran carries a
//!   claim naming the `teprob.f` line range it came from. A claim is its own
//!   comment line, marked `@port` and followed by the range, for example a
//!   line reading `@port teprob.f:505-522` above the reaction kinetics.
//!   `cargo xtask provenance` collects those claims and reports any part of the
//!   original that nothing accounts for. The marker must be anchored at the
//!   start of the comment: prose mentioning the convention, like this
//!   paragraph, must never be counted as coverage.
//!
//! # Status
//!
//! Skeleton. The model lands over phases 1 and 2; see `BACKLOG.org`.

#![no_std]
#![forbid(unsafe_code)]

extern crate alloc;

#[cfg(test)]
pub(crate) mod testing;

pub mod component;
pub mod constants;
pub mod rng;
pub mod stream;
pub mod variables;

pub use component::{ByComponent, Component, Composition};
pub use rng::TepRng;
pub use stream::Stream;
pub use variables::{Analyzer, MeasIndex, MvIndex, Unit, ValveId};
