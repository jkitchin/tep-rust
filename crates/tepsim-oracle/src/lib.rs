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
// Deliberately not gated on `oracle`: `cargo xtask fidelity` parses the golden
// trace on machines with no Fortran toolchain at all.
pub mod golden;
// Likewise: the Tier 1 pools and the ULP report are pure Rust, so they can be
// proved correct without a Fortran compiler present.
pub mod tier1;

#[cfg(feature = "oracle")]
mod ffi;

#[cfg(feature = "oracle")]
pub use ffi::{Const, Teproc, Wlk};
#[cfg(feature = "oracle")]
pub use oracle::{N_STATES, Oracle};

#[cfg(feature = "oracle")]
mod oracle {
    use std::sync::{Mutex, MutexGuard, OnceLock};

    use crate::ffi::{self, Const, Teproc, Wlk};

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

        /// A snapshot of `COMMON/TEPROC/`, the plant's whole working set.
        pub fn teproc(&mut self) -> Teproc {
            // SAFETY: as above.
            unsafe { (&raw const ffi::teproc_).read() }
        }

        /// Overwrite `COMMON/TEPROC/`. This is how a differential test forces
        /// the Fortran into a chosen state before evaluating.
        pub fn set_teproc(&mut self, v: &Teproc) {
            // SAFETY: as above.
            unsafe { (&raw mut ffi::teproc_).write(*v) };
        }

        /// A snapshot of `COMMON/WLK/`, the disturbance random-walk state.
        pub fn wlk(&mut self) -> Wlk {
            // SAFETY: as above.
            unsafe { (&raw const ffi::wlk_).read() }
        }

        /// Overwrite `COMMON/WLK/`. Needed to reproduce a disturbance
        /// trajectory, since the walks carry state between calls.
        pub fn set_wlk(&mut self, v: &Wlk) {
            // SAFETY: as above.
            unsafe { (&raw mut ffi::wlk_).write(*v) };
        }

        /// `COMMON/CONST/`, the thermodynamic coefficients. Read-only in
        /// practice: `TEINIT` sets them and nothing else writes them.
        pub fn constants(&mut self) -> Const {
            // SAFETY: as above.
            unsafe { (&raw const ffi::const_).read() }
        }

        /// The shutdown flag, non-zero when the plant has tripped.
        ///
        /// Only reachable because `build.rs` hoists `ISD` into
        /// `COMMON/SHUTDN/`; the original keeps it as a local. A test proves
        /// that rewrite changes no numbers.
        pub fn shutdown_flag(&mut self) -> i32 {
            // SAFETY: as above.
            unsafe { (&raw const ffi::shutdn_.isd).read() }
        }

        /// `TESUB1`: mixture enthalpy of composition `z` at `t` degrees C.
        ///
        /// `ity` selects the correlation: 0 for liquid, 1 for vapour, 2 for
        /// vapour with the ideal-gas correction subtracted.
        pub fn tesub1(&mut self, z: &[f64; 8], t: f64, ity: i32) -> f64 {
            let mut h = 0.0;
            // SAFETY: `z` is exactly the 8 elements the Fortran indexes.
            unsafe { ffi::tesub1_(z.as_ptr(), &t, &mut h, &ity) };
            h
        }

        /// `TESUB2`: temperature from enthalpy, by Newton iteration.
        ///
        /// `t` is both the initial guess and, on failure to converge in 100
        /// iterations, the value returned unchanged. The original reports no
        /// error in that case; quantifying how often it happens is B-0011.
        pub fn tesub2(&mut self, z: &[f64; 8], t_guess: f64, h: f64, ity: i32) -> f64 {
            let mut t = t_guess;
            // SAFETY: as above.
            unsafe { ffi::tesub2_(z.as_ptr(), &mut t, &h, &ity) };
            t
        }

        /// `TESUB3`: heat capacity, the temperature derivative of `TESUB1`.
        pub fn tesub3(&mut self, z: &[f64; 8], t: f64, ity: i32) -> f64 {
            let mut dh = 0.0;
            // SAFETY: as above.
            unsafe { ffi::tesub3_(z.as_ptr(), &t, &mut dh, &ity) };
            dh
        }

        /// `TESUB4`: liquid density of composition `x` at `t` degrees C.
        pub fn tesub4(&mut self, x: &[f64; 8], t: f64) -> f64 {
            let mut r = 0.0;
            // SAFETY: as above.
            unsafe { ffi::tesub4_(x.as_ptr(), &t, &mut r) };
            r
        }

        /// `TESUB6`: one measurement-noise sample of standard deviation `std`.
        ///
        /// Consumes twelve draws from the generator, summing them and
        /// subtracting 6 to approximate a Gaussian.
        pub fn tesub6(&mut self, std: f64) -> f64 {
            let mut x = 0.0;
            // SAFETY: no arrays involved.
            unsafe { ffi::tesub6_(&std, &mut x) };
            x
        }

        /// `TESUB7`: one draw. Negative `i` gives [-1,1), otherwise [0,1).
        pub fn tesub7(&mut self, i: i32) -> f64 {
            // SAFETY: no arrays involved.
            unsafe { ffi::tesub7_(&i) }
        }

        /// `TESUB8`: evaluate disturbance walk channel `i` (1-based) at time
        /// `t`, as a cubic in the time since the channel's last knot.
        pub fn tesub8(&mut self, i: i32, t: f64) -> f64 {
            // SAFETY: no arrays involved.
            unsafe { ffi::tesub8_(&i, &t) }
        }
    }
}
