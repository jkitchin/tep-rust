//! Development-only harness that drives the *original* Tennessee Eastman
//! Fortran so the Rust port can be differentially tested against it.
//!
//! This is the instrument the whole validation ladder rests on. It compiles the
//! unmodified `reference/fortran/teprob.f` with gfortran and lets a test force
//! the Fortran into a chosen state, evaluate once, and compare against the Rust
//! port under identical conditions.
//!
//! # This crate is never shipped
//!
//! `publish = false`, and nothing in `tepsim`, `tepsim-py` or `tepsim-wasm` may
//! depend on it. `cargo xtask ci` asserts that. The `oracle` feature is off by
//! default; without it this crate is an empty shell that needs no Fortran
//! compiler, so a contributor with no gfortran can still build the workspace.
//!
//! # The Fortran has exactly one global state
//!
//! `COMMON` blocks are process-wide mutable globals, so the original supports
//! precisely one simulation per process and is not reentrant. That is a large
//! part of *why* the Rust port exists, and it is also a hazard here: Rust runs
//! tests on multiple threads by default. `Oracle` is therefore a process-wide
//! singleton guarded by a mutex, obtained through `Oracle::lock`. There is no
//! way to touch the Fortran without holding it.
//!
//! (Those two are code spans rather than links because the type only exists
//! under the `oracle` feature, and this page is built without it.)
//!
//! # Reference numbers depend on the compiler
//!
//! Tier 1 and Tier 2 tolerances are measured against a specific gfortran build
//! with pinned flags, asserted by [`build_info`]. Changing either invalidates
//! every number recorded in `LOG.org`, so it is a deliberate, logged
//! re-baseline. Never add `-ffast-math` or anything else permitting
//! reassociation.

// FFI to a Fortran library cannot be written in safe Rust. The unsafety is
// confined to `ffi` and to the accessors below, each of which states the
// argument for why it is sound. See CLAUDE.md for the project-wide rule.
#![deny(unsafe_op_in_unsafe_fn)]

pub mod build_info;

#[cfg(feature = "oracle")]
mod ffi;

#[cfg(feature = "oracle")]
pub use oracle::{N_STATES, Oracle};

#[cfg(feature = "oracle")]
mod oracle {
    use std::sync::{Mutex, MutexGuard, OnceLock};

    use crate::ffi;

    /// The original integrates 50 states. `teprob.f:24-26` invites callers to
    /// append their own beyond that; we never do, so this is exact.
    pub const N_STATES: usize = 50;

    /// Exclusive access to the Fortran's process-wide `COMMON` state.
    ///
    /// Held for as long as you need a consistent view. Every accessor takes
    /// `&mut self`, so the borrow checker prevents interleaving two logical
    /// operations even within one thread.
    #[derive(Debug)]
    pub struct Oracle {
        _private: (),
    }

    fn mutex() -> &'static Mutex<Oracle> {
        static LOCK: OnceLock<Mutex<Oracle>> = OnceLock::new();
        LOCK.get_or_init(|| Mutex::new(Oracle { _private: () }))
    }

    impl Oracle {
        /// Take exclusive control of the Fortran.
        ///
        /// Blocks until any other holder is done. Recovers from a poisoned
        /// mutex rather than panicking: a test that failed while holding the
        /// lock left the Fortran in an arbitrary state, not an unsound one, and
        /// every entry point here overwrites the state it depends on.
        pub fn lock() -> MutexGuard<'static, Oracle> {
            mutex()
                .lock()
                .unwrap_or_else(|poisoned| poisoned.into_inner())
        }

        /// `TEINIT`: load the nominal steady state. Returns `(time, states)`.
        pub fn init(&mut self) -> (f64, [f64; N_STATES]) {
            let nn = N_STATES as i32;
            let mut time = 0.0_f64;
            let mut yy = [0.0_f64; N_STATES];
            let mut yp = [0.0_f64; N_STATES];
            // SAFETY: `nn` matches the length of both arrays, so the Fortran
            // cannot index past them. We hold the lock, so nothing else is
            // touching the COMMON blocks concurrently.
            unsafe { ffi::teinit_(&nn, &mut time, yy.as_mut_ptr(), yp.as_mut_ptr()) };
            (time, yy)
        }

        /// `TEFUNC`: evaluate derivatives at `time` for state `yy`.
        ///
        /// Also advances the disturbance walks, draws measurement noise and
        /// ticks the sampled analysers, so calling it twice with the same
        /// arguments does not generally give the same answer. Reset
        /// [`Oracle::set_rng`] first if you need repeatability.
        pub fn derivatives(&mut self, time: f64, yy: &[f64; N_STATES]) -> [f64; N_STATES] {
            let nn = N_STATES as i32;
            let mut state = *yy;
            let mut yp = [0.0_f64; N_STATES];
            // SAFETY: as above. `state` is a private copy, so the Fortran
            // writing through `yy` cannot surprise the caller.
            unsafe { ffi::tefunc_(&nn, &time, state.as_mut_ptr(), yp.as_mut_ptr()) };
            yp
        }

        /// `XMEAS(1..41)`, the measurement vector, zero-based here.
        pub fn measurements(&mut self) -> [f64; 41] {
            // SAFETY: we hold the lock, and `Pv` mirrors the Fortran layout of
            // `COMMON/PV/`, so this reads initialised f64s written by Fortran.
            unsafe { (&raw const ffi::pv_.xmeas).read() }
        }

        /// `XMV(1..12)`, the manipulated variables, zero-based here.
        pub fn manipulated(&mut self) -> [f64; 12] {
            // SAFETY: as above.
            unsafe { (&raw const ffi::pv_.xmv).read() }
        }

        /// Overwrite `XMV(1..12)`. This is how a controller is simulated.
        pub fn set_manipulated(&mut self, xmv: &[f64; 12]) {
            // SAFETY: as above, writing a correctly typed and sized value.
            unsafe { (&raw mut ffi::pv_.xmv).write(*xmv) };
        }

        /// The twenty `IDV` disturbance flags.
        pub fn disturbances(&mut self) -> [i32; 20] {
            // SAFETY: as above.
            unsafe { (&raw const ffi::dvec_.idv).read() }
        }

        /// Set the twenty `IDV` flags. `TEFUNC` clamps anything positive to 1.
        pub fn set_disturbances(&mut self, idv: &[i32; 20]) {
            // SAFETY: as above.
            unsafe { (&raw mut ffi::dvec_.idv).write(*idv) };
        }

        /// The RNG state: `COMMON/RANDSD/ G`, a single `f64`.
        pub fn rng(&mut self) -> f64 {
            // SAFETY: as above.
            unsafe { (&raw const ffi::randsd_.g).read() }
        }

        /// Set the RNG state. Necessary before any comparison that must be
        /// reproducible, and the basis of the Tier 3 call-order diff.
        pub fn set_rng(&mut self, g: f64) {
            // SAFETY: as above.
            unsafe { (&raw mut ffi::randsd_.g).write(g) };
        }
    }
}
