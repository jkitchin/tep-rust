//! Development-only harness that drives the *original* Tennessee Eastman
//! Fortran so the Rust port can be differentially tested against it.
//!
//! This is the instrument the whole validation ladder rests on. It compiles the
//! unmodified `teprob.f` with gfortran and exposes `TEINIT`, `TEFUNC`, the eight
//! `TESUB` routines, and read/write access to every `COMMON` block, so a test
//! can force both implementations into an identical state, evaluate once, and
//! compare.
//!
//! # This crate is never shipped
//!
//! `publish = false`, and nothing in `tepsim`, `tepsim-py` or `tepsim-wasm` may
//! depend on it. `cargo xtask ci` asserts that. Enabling the `oracle` feature
//! requires gfortran and the vendored sources under `reference/`.
//!
//! # Reference numbers depend on the compiler
//!
//! Tier 1 and Tier 2 tolerances are measured against a specific gfortran build
//! with pinned flags. Changing either invalidates every recorded number, so it
//! is a deliberate, logged re-baseline rather than a casual edit. Never add
//! `-ffast-math` or anything else permitting reassociation.
//!
//! # Status
//!
//! Skeleton. The FFI shim lands in B-0004; see `BACKLOG.org`.

#![forbid(unsafe_code)]
