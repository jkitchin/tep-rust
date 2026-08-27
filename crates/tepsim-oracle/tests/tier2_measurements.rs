//! Tier 2 for the measurements and the shutdown detector:
//! `tepsim_core::measurements` against `teprob.f:679-710`, over all three
//! sampling pools.
//!
//! Run twice, default and `--features oracle,libm-system`; see
//! `tier2_equilibrium.rs`.
//!
//! # Getting a noise-free `XMEAS` out of the Fortran
//!
//! `teprob.f:711-716` adds a random draw to all twenty-two measurements, so
//! `COMMON/PV/XMEAS` after an ordinary call is not what this item computes.
//! The guard is the way in:
//!
//! ```fortran
//!       IF(TIME.GT.0.0.AND.ISD.EQ.0)THEN
//! ```
//!
//! At `TIME = 0` no noise is drawn at all, so every scenario here is forced at
//! time zero and the block is skipped. Subtracting the noise instead is not an
//! option: it is reproducible, but recovering it would mean reimplementing
//! `TESUB6` and the draw order, which is B-0024b's job and would make this
//! test depend on the thing it is meant to be independent of.
//!
//! ## What forcing at time zero costs, and why it is acceptable here
//!
//! `teprob.f:397-406` resets the whole walk state whenever `TIME` is zero, so
//! every disturbance channel sits at `SZERO` and the walk-driven inputs are all
//! nominal. `tier2/mod.rs` warns about exactly this, and the warning is about
//! pools that *silently* lose walk coverage.
//!
//! It is acceptable for this item and not for the others, for a specific
//! reason: this range is a pure unit-conversion layer over quantities the
//! earlier items already validated across the full walk range, and the eight
//! shutdown limits read `PTR`, `VLR`, `TCR`, `VLS` and `VLC`, none of which
//! any walk channel touches. State diversity is what this item needs, and all
//! 2,417 states are still here.
//!
//! The shutdown comparison inherits the same window for free, since `ISD` is
//! computed at `teprob.f:702-710`, before the noise block either way.

#![cfg(feature = "oracle")]

use std::collections::BTreeSet;

use tepsim_core::measurements::ShutdownCause;
use tepsim_core::{
    State, equilibrium, flows, heat, math, measurements, streams, stripper, vessels,
};
use tepsim_oracle::Oracle;
use tepsim_oracle::tier1::Comparison;
use tepsim_oracle::tier2::{Case, Pool, Pools, Scenario, Snapshot, adversarial, compare_field};

const DT: f64 = 1.0 / 3600.0;

/// `PLAN.org`, "Tier 2".
const TIER2_TOLERANCE: f64 = 1e-12;

/// How many continuous measurements this item covers. `XMEAS(23..41)` are the
/// analysers and belong to B-0024b.
const CONTINUOUS: usize = 22;

/// The same scenario with its clock at zero; see the module documentation.
fn at_time_zero(scenario: &Scenario) -> Scenario {
    Scenario {
        time: 0.0,
        ..scenario.clone()
    }
}

fn solve(oracle: &mut Oracle, scenario: &Scenario) -> Option<measurements::Measured> {
    let t = scenario.time;
    let raw = |n: usize| f64::from(scenario.disturbances[n - 1]);
    let a = oracle.tesub8(1, t) - raw(1) * 0.03;
    let b = oracle.tesub8(2, t) + raw(2) * 0.005;
    let feed = streams::FeedConditions {
        ac_feed_light: [a, b, 1.0 - a - b],
        d_feed_celsius: oracle.tesub8(3, t) + raw(3) * 5.0,
        ac_feed_celsius: oracle.tesub8(4, t),
    };
    let flow_drift = flows::FlowDrift {
        steam_capacity: oracle.tesub8(9, t),
        reactor_outlet: oracle.tesub8(12, t),
    };
    let heat_drift = heat::HeatDrift {
        reactor_coolant: oracle.tesub8(10, t),
        condenser_coolant: oracle.tesub8(11, t),
    };
    let seeds = vessels::TemperatureSeeds {
        reactor: scenario.common.tcr,
        separator: scenario.common.tcs,
        stripper: scenario.common.tcc,
        mixing: scenario.common.tcv,
    };
    let state = State::from_flat(&scenario.state);
    let unpacked = vessels::unpack(&state, seeds).ok()?;
    let eq = equilibrium::equilibrium(&unpacked);
    let mut table = streams::streams(&unpacked, &eq, &feed);
    let mut idv = [0.0; 20];
    for (slot, r) in idv.iter_mut().zip(scenario.disturbances) {
        *slot = f64::from(r.clamp(0, 1));
    }
    let mut flow = flows::flows(&state, &unpacked, &eq, &table, &idv, flow_drift);
    let _ = stripper::stripper(&mut table, &mut flow, unpacked.stripper.celsius);
    let hot = heat::heat_transfer(&state, &unpacked, &table, &flow, heat_drift);
    Some(measurements::measurements(
        &state,
        &unpacked,
        &table,
        &flow,
        &hot,
        (
            eq.reactor.pressure,
            eq.separator.pressure,
            eq.mixing_pressure,
        ),
    ))
}

/// Force at time zero, solve, and record. Returns the snapshot too, so the
/// shutdown flag can be compared from the same evaluation.
fn observe(
    oracle: &mut Oracle,
    scenario: &Scenario,
    pool: Pool,
    index: usize,
    comparison: &mut Comparison<Case>,
) -> Option<(Snapshot, measurements::Measured)> {
    let zeroed = at_time_zero(scenario);
    let snapshot = zeroed.force(oracle);
    let solved = solve(oracle, &zeroed)?;
    compare_field(
        comparison,
        pool,
        index,
        &solved.continuous,
        &snapshot.measurements[..CONTINUOUS],
    );
    Some((snapshot, solved))
}

/// Every scenario in the three pools, at time zero.
fn all_scenarios(
    oracle: &mut Oracle,
    steps: usize,
    perturbations: usize,
    seed: u64,
) -> Vec<(Pool, usize, Scenario)> {
    let pools = Pools::collect(oracle, steps, DT);
    let mut out = Vec::new();
    for index in 0..pools.trajectory.len() {
        out.push((Pool::Nominal, index, pools.nominal_case(index)));
    }
    let mut sampler = tepsim_oracle::tier1::Sampler::new(seed);
    for index in 0..perturbations {
        out.push((
            Pool::Perturbed,
            index,
            pools.perturbed_case(index, &mut sampler),
        ));
    }
    let (boundaries, missed) = adversarial::build(oracle, &pools.nominal_case(0));
    assert!(missed.is_empty(), "the catalogue regressed: {missed:?}");
    for (index, boundary) in boundaries.iter().enumerate() {
        out.push((Pool::Adversarial, index, boundary.scenario.clone()));
    }
    out
}

#[test]
fn the_measurements_match_the_fortran_over_all_three_pools() {
    let mut oracle = Oracle::lock();
    let scenarios = all_scenarios(&mut oracle, 400, 2_000, 0x7E2_0024);
    let mut comparison: Comparison<Case> = Comparison::new("XMEAS(1..22), noise-free");
    let mut skipped = 0;

    for (pool, index, scenario) in &scenarios {
        if observe(&mut oracle, scenario, *pool, *index, &mut comparison).is_none() {
            skipped += 1;
        }
    }

    println!(
        "transcendentals come from the {} libm",
        if math::USES_SYSTEM_LIBM {
            "platform"
        } else {
            "vendored"
        }
    );
    println!("{comparison}");
    println!("states skipped because the port would not converge: {skipped}");
    assert_eq!(
        skipped, 0,
        "{skipped} states failed to converge in the port"
    );
    comparison.assert_within(TIER2_TOLERANCE);
}

/// With the transcendentals taken out of the comparison, the algebra must be
/// bit-identical. See `tier2_equilibrium.rs`.
#[test]
#[cfg(feature = "libm-system")]
fn the_algebra_is_bit_identical_once_exp_and_pow_agree() {
    let mut oracle = Oracle::lock();
    let scenarios = all_scenarios(&mut oracle, 200, 500, 7);
    let mut comparison: Comparison<Case> = Comparison::new("XMEAS(1..22), noise-free");

    for (pool, index, scenario) in &scenarios {
        observe(&mut oracle, scenario, *pool, *index, &mut comparison);
    }
    println!("{comparison}");
    assert_eq!(
        comparison.max_ulp(),
        0,
        "the measurements are not bit-identical under the platform libm.\n\
         Either the conversions in `tepsim_core::measurements` do not \
         associate the way `teprob.f:679-701` does, which is a port bug; or \
         gfortran no longer resolves the transcendentals to the same code as \
         the platform libm."
    );
}

/// The shutdown detector must agree with `ISD` on every state, in both
/// directions.
///
/// This is the one place a disagreement would be invisible in the derivative
/// comparison rather than obvious: `teprob.f:807-811` zeroes all fifty
/// derivatives on a trip, so a port that tripped when the Fortran did not
/// would return zeros that look like a perfectly valid frozen plant.
#[test]
fn the_shutdown_detector_agrees_with_the_fortran_on_every_state() {
    let mut oracle = Oracle::lock();
    let scenarios = all_scenarios(&mut oracle, 400, 2_000, 0x7E2_0024);
    let mut comparison: Comparison<Case> = Comparison::new("ISD as 0 or 1");
    let (mut tripped, mut healthy) = (0, 0);

    for (pool, index, scenario) in &scenarios {
        let Some((snapshot, solved)) = observe(
            &mut oracle,
            scenario,
            *pool,
            *index,
            &mut Comparison::new("discarded"),
        ) else {
            continue;
        };
        if snapshot.tripped {
            tripped += 1;
        } else {
            healthy += 1;
        }
        comparison.observe(
            Case {
                pool: *pool,
                index: *index,
                component: 1,
            },
            f64::from(u8::from(solved.shutdown.is_tripped())),
            f64::from(u8::from(snapshot.tripped)),
        );
    }

    println!("{comparison}");
    println!("{tripped} states trip, {healthy} do not");
    assert_eq!(
        comparison.max_ulp(),
        0,
        "the port and the Fortran disagree about whether the plant is down"
    );
    assert!(tripped > 0, "no state tripped, so the detector is untested");
    assert!(healthy > 0, "every state tripped, which cannot be right");
}

/// Every one of the eight causes must be reached by the pool.
///
/// The original throws the reason away, so the oracle cannot confirm *which*
/// condition fired, only that one did. What can be checked, and is, is that
/// the port never claims a cause on a state the Fortran calls healthy, and
/// never calls a state healthy that the Fortran trips. Coverage of the eight
/// individual causes is then meaningful rather than decorative.
#[test]
fn the_pool_reaches_several_distinct_shutdown_causes() {
    let mut oracle = Oracle::lock();
    let scenarios = all_scenarios(&mut oracle, 400, 2_000, 0x7E2_0024);
    let mut seen: BTreeSet<&'static str> = BTreeSet::new();

    for (pool, index, scenario) in &scenarios {
        let Some((snapshot, solved)) = observe(
            &mut oracle,
            scenario,
            *pool,
            *index,
            &mut Comparison::new("discarded"),
        ) else {
            continue;
        };
        assert_eq!(
            solved.shutdown.is_tripped(),
            snapshot.tripped,
            "disagreement at {pool}#{index}"
        );
        for cause in solved.shutdown.causes() {
            seen.insert(cause.describe());
        }
    }

    // Which boundary produces which cause, since the oracle's bare flag
    // cannot say and the mapping is not obvious.
    let pools = Pools::collect(&mut oracle, 100, DT);
    let (boundaries, _) = adversarial::build(&mut oracle, &pools.nominal_case(0));
    for boundary in &boundaries {
        let zeroed = at_time_zero(&boundary.scenario);
        let snapshot = zeroed.force(&mut oracle);
        let Some(solved) = solve(&mut oracle, &zeroed) else {
            continue;
        };
        if snapshot.tripped {
            let causes: Vec<&str> = solved.shutdown.causes().map(|c| c.describe()).collect();
            println!("{:52} {causes:?}", boundary.target.name);
        }
    }

    println!("shutdown causes reached: {seen:?}");
    for cause in ShutdownCause::ALL {
        assert!(
            seen.contains(cause.describe()),
            "'{}' was never reached by any state in the pool, so that limit is \
             implemented and never validated against the oracle. Every limit \
             needs a state *past* it, not only on it: all eight comparisons \
             are strict. B-0024a added four such states for exactly this \
             reason; if one no longer lands, the catalogue regressed.",
            cause.describe()
        );
    }
}

/// Forcing at time zero really does skip the noise, and the noise really is
/// there otherwise.
///
/// The whole comparison rests on that guard, so it is checked directly rather
/// than assumed. If `teprob.f:711` ever changed, every number in this file
/// would silently become a comparison against noisy data.
#[test]
fn time_zero_is_what_suppresses_the_noise() {
    let mut oracle = Oracle::lock();
    let pools = Pools::collect(&mut oracle, 50, DT);

    // A state part-way along the trajectory, so its natural time is non-zero.
    let scenario = pools.nominal_case(30);
    assert!(scenario.time > 0.0, "the pool did not advance the clock");

    let noisy = scenario.force(&mut oracle);
    let clean = at_time_zero(&scenario).force(&mut oracle);
    assert!(!clean.tripped, "a tripped state skips the noise too");

    let differing = (0..CONTINUOUS)
        .filter(|i| noisy.measurements[*i].to_bits() != clean.measurements[*i].to_bits())
        .count();
    println!("{differing} of {CONTINUOUS} measurements differ with the clock running");
    assert!(
        differing >= CONTINUOUS - 2,
        "only {differing} measurements changed when the clock ran, so \
         teprob.f:711-716 is not adding noise the way this file assumes"
    );

    // And two evaluations at time zero agree exactly, which is what makes the
    // comparison reproducible.
    let again = at_time_zero(&scenario).force(&mut oracle);
    for i in 0..CONTINUOUS {
        assert_eq!(
            clean.measurements[i].to_bits(),
            again.measurements[i].to_bits(),
            "XMEAS({}) is not reproducible at time zero",
            i + 1
        );
    }
}
