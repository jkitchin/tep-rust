//! Proves the original Fortran can be driven from Rust: the milestone the whole
//! validation ladder depends on.
//!
//! If these fail, nothing downstream can be trusted, because every later tier
//! compares the Rust port against exactly this.

// Differential tests: exact comparisons are the property under test, and
// arithmetic is transcribed from `teprob.f` so a reader can check it against
// the listing line by line. Rearranging either would defeat the point.
#![allow(
    clippy::float_cmp,
    reason = "bit equality against the Fortran is the property under test"
)]
#![allow(
    clippy::suboptimal_flops,
    reason = "expressions are transcribed to be checkable against teprob.f"
)]
#![cfg(feature = "oracle")]

use tepsim_oracle::{N_STATES, Oracle};

/// Nominal steady state, from the original's documentation and from thirty
/// years of published TEP results. Loose bounds on purpose: this asserts "the
/// Fortran is wired up and running the plant we think it is", not a tolerance.
const NOMINAL_REACTOR_TEMP_C: (f64, f64) = (119.0, 122.0);
const NOMINAL_REACTOR_PRESSURE_KPA: (f64, f64) = (2600.0, 2800.0);

#[test]
fn teinit_loads_the_nominal_steady_state() {
    let mut oracle = Oracle::lock();
    let (time, yy) = oracle.init();

    assert_eq!(time, 0.0, "TEINIT sets TIME to zero");
    assert!(
        yy.iter().all(|v| v.is_finite()),
        "every one of the 50 initial states must be finite"
    );

    // The last twelve states are valve positions, as percentages.
    for (i, pos) in yy[38..50].iter().enumerate() {
        assert!(
            (0.0..=100.0).contains(pos),
            "valve {} initial position {pos} outside 0-100%",
            i + 1
        );
    }
}

#[test]
fn tefunc_returns_fifty_finite_derivatives() {
    let mut oracle = Oracle::lock();
    let (time, yy) = oracle.init();
    let yp = oracle.derivatives(time, &yy);

    assert_eq!(yp.len(), N_STATES);
    assert!(
        yp.iter().all(|v| v.is_finite()),
        "derivatives must all be finite at the nominal state"
    );
    assert!(
        yp.iter().any(|v| *v != 0.0),
        "an all-zero derivative vector would mean the plant is frozen, which \
         happens only on a shutdown trip and must not happen at t=0"
    );
}

#[test]
fn the_nominal_operating_point_is_the_published_one() {
    let mut oracle = Oracle::lock();
    let (time, yy) = oracle.init();
    let _ = oracle.derivatives(time, &yy);
    let xmeas = oracle.measurements();

    // XMEAS(7) reactor pressure, XMEAS(9) reactor temperature. One-based there,
    // zero-based here.
    let pressure = xmeas[6];
    let temperature = xmeas[8];

    assert!(
        (NOMINAL_REACTOR_PRESSURE_KPA.0..=NOMINAL_REACTOR_PRESSURE_KPA.1).contains(&pressure),
        "reactor pressure {pressure} kPa is not the published nominal ~2705"
    );
    assert!(
        (NOMINAL_REACTOR_TEMP_C.0..=NOMINAL_REACTOR_TEMP_C.1).contains(&temperature),
        "reactor temperature {temperature} C is not the published nominal ~120.4"
    );
}

/// The Fortran's generator, `teprob.f:1547-1555`.
///
/// Reproduced here only as a counting device, to work out how many draws a call
/// consumed. The real port, and the argument for why this arithmetic is
/// bit-reproducible despite overflowing 2^53, is B-0005.
fn lcg_next(g: f64) -> f64 {
    (g * 9228907.0) % 4294967296.0
}

/// How many draws separate two observed generator states, or `None` if `to` is
/// not reachable from `from` within `max` draws.
fn draws_between(from: f64, to: f64, max: usize) -> Option<usize> {
    let mut g = from;
    for n in 0..=max {
        if g == to {
            return Some(n);
        }
        g = lcg_next(g);
    }
    None
}

#[test]
fn the_rng_state_round_trips() {
    let mut oracle = Oracle::lock();
    const SEED: f64 = 4651207995.0; // teprob.f:1187
    oracle.set_rng(SEED);
    assert_eq!(
        oracle.rng(),
        SEED,
        "RNG state must round-trip through COMMON"
    );
}

/// At exactly t=0 the original draws nothing at all.
///
/// Three guards line up to make this true, and it is easy to assume otherwise.
/// `TEINIT` sets `TNEXT(I)=0.1` at `teprob.f:1362`, so both disturbance-walk
/// loops at `teprob.f:359` and `teprob.f:372` are skipped. The measurement
/// noise loop at `teprob.f:711` is guarded by `TIME.GT.0.0`. So the generator
/// is untouched, and the noise sequence begins on the *second* step.
#[test]
fn no_draws_happen_at_time_zero() {
    let mut oracle = Oracle::lock();
    const SEED: f64 = 4651207995.0;

    let (time, yy) = oracle.init();
    assert_eq!(time, 0.0);
    oracle.set_rng(SEED);
    let _ = oracle.derivatives(time, &yy);

    assert_eq!(
        oracle.rng(),
        SEED,
        "TEFUNC must not draw at t=0: measurement noise is gated on TIME > 0 \
         and the walk knots are not due until t=0.1"
    );
}

/// The first step past zero draws exactly 22 x 12 = 264 uniforms.
///
/// `TESUB6` sums twelve uniforms per measurement to approximate a Gaussian
/// (`teprob.f:1539-1546`), and the noise loop covers the 22 continuous
/// measurements. The composition analysers do not fire until t=0.1 and t=0.25,
/// so they contribute nothing yet. This exact count is what Tier 3's
/// call-order diff will have to reproduce.
#[test]
fn the_first_step_past_zero_draws_exactly_264() {
    let mut oracle = Oracle::lock();
    const SEED: f64 = 4651207995.0;
    const ONE_SECOND_H: f64 = 1.0 / 3600.0;

    let (_, yy) = oracle.init();
    oracle.set_rng(SEED);
    let _ = oracle.derivatives(ONE_SECOND_H, &yy);

    let draws = draws_between(SEED, oracle.rng(), 2000);
    assert_eq!(
        draws,
        Some(22 * 12),
        "expected 22 noisy measurements x 12 uniforms each"
    );
}

#[test]
fn setting_the_rng_makes_a_step_reproducible() {
    let mut oracle = Oracle::lock();
    const SEED: f64 = 4651207995.0;

    let (time, yy) = oracle.init();
    oracle.set_rng(SEED);
    let first = oracle.derivatives(time, &yy);
    let first_g = oracle.rng();

    let (time2, yy2) = oracle.init();
    oracle.set_rng(SEED);
    let second = oracle.derivatives(time2, &yy2);

    assert_eq!(
        first, second,
        "with the same state and the same RNG word, TEFUNC must be bit-identical; \
         without that, no differential test downstream can mean anything"
    );
    assert_eq!(
        first_g,
        oracle.rng(),
        "and must consume the same number of draws"
    );
}

#[test]
fn disturbance_flags_round_trip() {
    let mut oracle = Oracle::lock();
    oracle.set_disturbances(&[0; 20]);
    assert_eq!(oracle.disturbances(), [0; 20]);

    let mut idv = [0; 20];
    idv[3] = 1; // IDV(4): reactor cooling water inlet temperature step
    oracle.set_disturbances(&idv);
    assert_eq!(oracle.disturbances()[3], 1);

    oracle.set_disturbances(&[0; 20]);
}

#[test]
fn a_disturbance_actually_perturbs_the_plant() {
    let mut oracle = Oracle::lock();
    const SEED: f64 = 4651207995.0;

    let (time, yy) = oracle.init();
    oracle.set_disturbances(&[0; 20]);
    oracle.set_rng(SEED);
    let baseline = oracle.derivatives(time, &yy);

    let (time2, yy2) = oracle.init();
    let mut idv = [0; 20];
    idv[3] = 1; // IDV(4)
    oracle.set_disturbances(&idv);
    oracle.set_rng(SEED);
    let disturbed = oracle.derivatives(time2, &yy2);

    oracle.set_disturbances(&[0; 20]);

    assert_ne!(
        baseline, disturbed,
        "IDV(4) raises reactor cooling water inlet temperature by 5 degrees, so \
         it must change the derivative vector. Equality would mean the IDV \
         array is not reaching the Fortran."
    );
}

#[test]
fn manipulated_variables_round_trip() {
    let mut oracle = Oracle::lock();
    let (_, _) = oracle.init();
    let original = oracle.manipulated();

    let mut xmv = original;
    xmv[9] = 42.0; // XMV(10): reactor cooling water flow
    oracle.set_manipulated(&xmv);
    assert_eq!(oracle.manipulated()[9], 42.0);

    oracle.set_manipulated(&original);
    assert_eq!(oracle.manipulated(), original);
}

/// `TEINIT` does not reset the Newton warm starts, so two identical runs in
/// one process disagree in the last bits.
///
/// `TCR`, `TCS`, `TCC` and `TCV` are read as guesses and written as answers by
/// `TESUB2` (`teprob.f:460`, `465`, and the stripper and mixing-zone calls
/// beside them). Nothing in `teprob.f` ever assigns them otherwise. That makes
/// them process-global mutable state that survives `TEINIT`, and it makes any
/// test that starts from the nominal state order-dependent.
///
/// [`Oracle::init_cold`] is the fix. This pins both halves: that the problem
/// is real, and that the fix removes it.
#[test]
fn teinit_does_not_reset_the_newton_warm_starts() {
    let mut oracle = Oracle::lock();

    let temperatures = |o: &mut Oracle| {
        let c = o.teproc();
        [
            c.tcr.to_bits(),
            c.tcs.to_bits(),
            c.tcc.to_bits(),
            c.tcv.to_bits(),
        ]
    };
    let run = |o: &mut Oracle, steps: usize| {
        let (_, mut yy) = o.init_cold();
        o.set_disturbances(&[0; 20]);
        let mut t = 0.0;
        for _ in 0..steps {
            let yp = o.derivatives(t, &yy);
            for (slot, rate) in yy.iter_mut().zip(yp) {
                *slot += rate / 3600.0;
            }
            t += 1.0 / 3600.0;
        }
    };

    let reference = {
        let (_, _) = oracle.init_cold();
        temperatures(&mut oracle)
    };

    // A plain `TEINIT` after a run does *not* come back to it.
    run(&mut oracle, 2_000);
    let (_, _) = oracle.init();
    assert_ne!(
        temperatures(&mut oracle),
        reference,
        "TEINIT restored the warm starts, so init_cold is unnecessary and this \
         test is wrong about teprob.f"
    );

    // `init_cold` does, from any history.
    for steps in [1, 137, 5_000] {
        run(&mut oracle, steps);
        let (_, _) = oracle.init_cold();
        assert_eq!(
            temperatures(&mut oracle),
            reference,
            "init_cold gave a different warm start after {steps} steps"
        );
    }
}
