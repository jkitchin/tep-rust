//! Compares the Rust [`TepRng`] against the live Fortran, draw by draw.
//!
//! Stronger than the committed vectors in two ways: it checks every draw rather
//! than a fold, and it interleaves the two output modes the way `TEFUNC`
//! actually does. Needs gfortran, so it runs in the non-fast gate.

#![cfg(feature = "oracle")]

use tepsim_core::TepRng;
use tepsim_oracle::Oracle;

/// Long enough to cover the first step's 264 noise draws many times over.
const DRAWS: usize = 200_000;

#[test]
fn every_draw_matches_the_fortran() {
    let mut oracle = Oracle::lock();

    for seed in [
        TepRng::DEFAULT_SEED,
        1_431_655_765.0,
        4_243_534_565.0,
        7_854_912_354.0,
    ] {
        oracle.set_rng(seed);
        let mut rust = TepRng::new(seed);

        for draw in 0..DRAWS {
            let fortran = oracle.tesub7(1);
            let ours = rust.unit();
            assert_eq!(
                ours.to_bits(),
                fortran.to_bits(),
                "seed {seed}, draw {draw}: Rust gives {ours:?}, Fortran gives \
                 {fortran:?}. The generator must be bit-exact or nothing \
                 downstream can be compared."
            );
            assert_eq!(
                rust.state().to_bits(),
                oracle.rng().to_bits(),
                "seed {seed}, draw {draw}: the raw state diverged"
            );
        }
    }
}

/// `TEFUNC` uses both modes, and `TESUB5` interleaves them within one call, so
/// check the interleaving rather than each mode in isolation.
#[test]
fn interleaved_modes_match_the_fortran() {
    let mut oracle = Oracle::lock();
    oracle.set_rng(TepRng::DEFAULT_SEED);
    let mut rust = TepRng::new(TepRng::DEFAULT_SEED);

    for draw in 0..50_000 {
        // An irregular pattern, so an implementation that happened to be right
        // only for alternating calls would still be caught.
        let signed = draw % 3 == 0 || draw % 7 == 0;
        let (ours, fortran) = if signed {
            (rust.signed(), oracle.tesub7(-1))
        } else {
            (rust.unit(), oracle.tesub7(1))
        };
        assert_eq!(
            ours.to_bits(),
            fortran.to_bits(),
            "draw {draw} (signed={signed}): {ours:?} vs {fortran:?}"
        );
    }
}

/// The sign of the argument is the mode selector, and zero counts as
/// non-negative. Easy to get backwards, and it would corrupt every disturbance.
#[test]
fn zero_selects_the_unit_scaling_not_the_signed_one() {
    let mut oracle = Oracle::lock();

    oracle.set_rng(TepRng::DEFAULT_SEED);
    let at_zero = oracle.tesub7(0);
    oracle.set_rng(TepRng::DEFAULT_SEED);
    let at_one = oracle.tesub7(1);
    oracle.set_rng(TepRng::DEFAULT_SEED);
    let at_minus_one = oracle.tesub7(-1);

    assert_eq!(
        at_zero.to_bits(),
        at_one.to_bits(),
        "TESUB7(0) must use the I.GE.0 branch"
    );
    assert_ne!(at_zero.to_bits(), at_minus_one.to_bits());

    let mut rust = TepRng::new(TepRng::DEFAULT_SEED);
    assert_eq!(rust.unit().to_bits(), at_zero.to_bits());
}
