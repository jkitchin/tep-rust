//! Which parts of `TEFUNC` draw from the generator, in what order, and how
//! many times.
//!
//! B-0027. Tier 3 exists because every consumer shares one generator word, so
//! a port can reproduce the right *distribution* from the wrong *stream* and
//! look correct in every statistic. That failure is invisible to Tier 2, which
//! pins the generator before each evaluation, and it is invisible to Tier 5,
//! which compares distributions. Only the order catches it.
//!
//! # The five consumers, read from the source
//!
//! | Line | Consumer | Draws per firing |
//! |---|---|---|
//! | 367 | `TESUB5`, walk channels 1-9 | 3 |
//! | 388 | the spike rule's dwell, channels 10-12 | 1 |
//! | 713 | measurement noise, `XMEAS(1..22)` | 12 each, so 264 |
//! | 747 | the gas analysers, `XMEAS(23..36)` | 12 each, so 168 |
//! | 756 | the product analyser, `XMEAS(37..41)` | 12 each, so 60 |
//!
//! `TESUB5` draws three because it calls `TESUB7` at `teprob.f:1528`, `1529`
//! and `1530`; `TESUB6` draws twelve, at `1542-1543`.
//!
//! Every one of those sites is *conditional*, which is the whole difficulty:
//!
//! - 367 fires only when `TIME >= TNEXT(I)`, per channel.
//! - 388 fires only when the interpolated value has fallen to 0.1 or below.
//! - 713 is skipped entirely at `TIME = 0` and on a tripped plant
//!   (`teprob.f:711`).
//! - 747 and 756 fire on their own schedules, 0.1 and 0.25 hours.
//!
//! So the draw count per evaluation is not a constant, and this file measures
//! it rather than asserting a number read off the listing.
//!
//! # How the count is recovered without instrumenting anything
//!
//! `TESUB7` is a multiplicative congruential generator with no state beyond
//! `G` (`teprob.f:1551`). So stepping a port-side generator forward from the
//! value before a call until it matches the value after gives the number of
//! draws exactly, with no edit to the Fortran and no counter to keep in sync.

// Differential tests: exact comparisons are the property under test, and
// arithmetic is transcribed from `teprob.f` so that a reader can check it
// against the listing line by line. Rearranging either would defeat the point.
#![allow(
    clippy::float_cmp,
    reason = "bit equality against the Fortran is the property under test"
)]
#![allow(
    clippy::suboptimal_flops,
    reason = "expressions are transcribed to be checkable against teprob.f"
)]
#![cfg(feature = "oracle")]

use tepsim_core::TepRng;
use tepsim_oracle::Oracle;
use tepsim_oracle::tier2::{Pools, Scenario};

const DT: f64 = 1.0 / 3600.0;

/// How many draws took `before` to `after`, or `None` if more than `limit`.
///
/// Exact, because the generator is a pure function of its own word.
fn draws_between(before: f64, after: f64, limit: usize) -> Option<usize> {
    if before == after {
        return Some(0);
    }
    let mut rng = TepRng::new(before);
    for step in 1..=limit {
        let _ = rng.unit();
        if rng.state() == after {
            return Some(step);
        }
    }
    None
}

/// Force a scenario and report how many draws the evaluation made.
fn draws_for(oracle: &mut Oracle, scenario: &Scenario) -> usize {
    oracle.set_teproc(&scenario.common);
    oracle.set_wlk(&scenario.walk);
    oracle.set_rng(scenario.rng);
    oracle.set_manipulated(&scenario.manipulated);
    oracle.set_disturbances(&scenario.disturbances);
    let before = oracle.rng();
    let _ = oracle.derivatives(scenario.time, &scenario.state);
    let after = oracle.rng();
    draws_between(before, after, 2_000).expect("fewer than two thousand draws")
}

/// The counting method itself must be right before anything is measured with
/// it. `TESUB7` advances `G` by exactly one step per call.
#[test]
fn the_draw_counter_agrees_with_the_oracle_generator() {
    let mut oracle = Oracle::lock();
    for seed in [TepRng::DEFAULT_SEED, 1.0, 1_431_655_765.0, 4_294_967_295.0] {
        for expected in [1_usize, 2, 12, 37] {
            oracle.set_rng(seed);
            for _ in 0..expected {
                let _ = oracle.tesub7(1);
            }
            let after = oracle.rng();
            assert_eq!(
                draws_between(seed, after, 100),
                Some(expected),
                "counting from {seed} after {expected} draws"
            );
        }
    }
}

/// At `TIME = 0` the noise block is skipped and the walks are reset, so the
/// only draws are whatever the walk advance itself makes.
///
/// This is the quiet case, and it is what `tier2_measurements.rs` relies on to
/// get a noise-free `XMEAS`.
#[test]
fn a_time_zero_evaluation_draws_far_less_than_a_running_one() {
    let mut oracle = Oracle::lock();
    let pools = Pools::collect(&mut oracle, 60, DT);

    let quiet = draws_for(
        &mut oracle,
        &Scenario {
            time: 0.0,
            ..pools.nominal_case(30)
        },
    );
    let running = draws_for(&mut oracle, &pools.nominal_case(30));

    println!("draws at TIME=0: {quiet}; at the scenario's own time: {running}");
    assert!(
        running >= quiet + 200,
        "a running evaluation drew {running} and a time-zero one {quiet}; the \
         noise block alone should account for 264, so teprob.f:711 is not \
         being skipped the way tier2_measurements.rs assumes"
    );
}

/// The measurement noise is 22 measurements times 12 draws each.
///
/// Measured as a difference between two evaluations that differ only in
/// whether the noise block runs, so nothing has to be instrumented.
#[test]
fn the_measurement_noise_costs_exactly_264_draws() {
    let mut oracle = Oracle::lock();
    let pools = Pools::collect(&mut oracle, 60, DT);

    // A time just past zero, so the noise block runs but no analyser schedule
    // has come due yet: TGAS starts at 0.1 and TPROD at 0.25.
    let scenario = pools.nominal_case(30);
    let just_after = Scenario {
        time: 1.0e-6,
        ..scenario.clone()
    };
    let at_zero = Scenario {
        time: 0.0,
        ..scenario
    };

    let with_noise = draws_for(&mut oracle, &just_after);
    let without = draws_for(&mut oracle, &at_zero);
    println!(
        "with noise {with_noise}, without {without}, difference {}",
        with_noise - without
    );
    assert_eq!(
        with_noise - without,
        22 * 12,
        "the noise block should be 22 measurements at 12 draws each"
    );
}

/// The analysers fire on their own schedules and each costs 12 draws per
/// composition, so the count jumps at 0.1 and again at 0.25 hours.
#[test]
fn the_analysers_add_their_own_draws_on_schedule() {
    let mut oracle = Oracle::lock();
    let pools = Pools::collect(&mut oracle, 60, DT);
    let base = pools.nominal_case(30);

    // `TGAS` and `TPROD` come out of the scenario's own COMMON, which
    // `Pools::collect` captured at t = 0, so they are 0.1 and 0.25.
    let at = |time: f64| Scenario {
        time,
        ..base.clone()
    };

    let before_gas = draws_for(&mut oracle, &at(0.05));
    let after_gas = draws_for(&mut oracle, &at(0.15));
    let after_product = draws_for(&mut oracle, &at(0.30));

    println!("draws: t=0.05 {before_gas}, t=0.15 {after_gas}, t=0.30 {after_product}");

    // The product analyser is clean: nothing else comes due between 0.15 and
    // 0.30, so the whole difference is XMEAS(37..41) at twelve draws each.
    assert_eq!(
        after_product - after_gas,
        5 * 12,
        "past 0.25 the product analyser adds XMEAS(37..41), five compositions"
    );

    // The gas analysers are *not* clean, and this is the finding. `TNEXT`
    // starts at 0.1 for all twelve walk channels (`teprob.f:1362`) and `TGAS`
    // starts at 0.1 as well (`teprob.f:741`), so the first gas sample and the
    // first walk advance come due on the same tick, forever after in step.
    //
    // 198 is not even a multiple of twelve, which is what gives it away: the
    // analysers alone could only ever contribute a multiple of twelve.
    let gas_and_walk = after_gas - before_gas;
    assert_eq!(gas_and_walk, 198);
    assert_ne!(
        gas_and_walk % 12,
        0,
        "the step is a clean multiple of twelve, so the walk no longer \
         coincides with the analyser and this test's premise has changed"
    );
    assert_eq!(
        gas_and_walk,
        14 * 12 + 30,
        "168 for the gas, 30 for the walk"
    );
}

/// The walk advance, isolated from the analyser it coincides with.
///
/// Pushing `TNEXT` past the horizon stops all twelve channels from firing,
/// which is the only way to separate them: `TNEXT` and `TGAS` both start at
/// 0.1 and stay in step.
///
/// The 30 draws split as 27 + 3: nine channels through `TESUB5` at three each
/// (`teprob.f:1528-1530`), plus one each from the three spike channels taking
/// the `SWLK <= 0.1` branch at `teprob.f:388`.
///
/// **That last part is data-dependent**, and it is the reason Tier 3 needs a
/// trace rather than a count. A spike channel above 0.1 draws *nothing*
/// (`teprob.f:381-385`), so the number of draws per evaluation depends on the
/// walk's own value, which depends on earlier draws.
#[test]
fn the_walk_advance_costs_thirty_draws_of_which_three_are_conditional() {
    let mut oracle = Oracle::lock();
    let pools = Pools::collect(&mut oracle, 60, DT);
    let base = pools.nominal_case(30);

    let frozen_walk = |time: f64| {
        let mut s = Scenario {
            time,
            ..base.clone()
        };
        // Nothing comes due before the horizon, so no channel advances.
        s.walk.tnext = [1.0e6; 12];
        s
    };

    let mut with_walk = Scenario {
        time: 0.15,
        ..base.clone()
    };
    with_walk.walk.tnext = base.walk.tnext;

    let without = draws_for(&mut oracle, &frozen_walk(0.15));
    let with = draws_for(&mut oracle, &with_walk);
    println!("t=0.15 with the walk frozen: {without}; running: {with}");

    assert_eq!(
        without,
        (22 + 14) * 12,
        "with the walk frozen only the noise and the gas analysers draw"
    );
    assert_eq!(with - without, 30, "the walk advance costs 30 draws here");
    assert_eq!(
        with - without,
        9 * 3 + 3,
        "nine channels through TESUB5 at three draws each, plus one each from \
         the three spike channels taking the dwell branch"
    );
}

/// A tripped plant draws nothing for noise: `teprob.f:711` guards on `ISD`.
///
/// So the *stream position* depends on whether the plant is down, which means
/// a port that got the shutdown detector wrong would desynchronise the
/// generator and every subsequent draw would differ. B-0024a validated the
/// detector on every state; this is why that mattered beyond the derivative.
#[test]
fn a_tripped_plant_skips_the_noise_and_so_the_stream_position_depends_on_it() {
    let mut oracle = Oracle::lock();
    let pools = Pools::collect(&mut oracle, 60, DT);
    let (boundaries, _) =
        tepsim_oracle::tier2::adversarial::build(&mut oracle, &pools.nominal_case(0));

    let mut tripped_counts = Vec::new();
    let mut healthy_counts = Vec::new();
    for boundary in &boundaries {
        // Give it a time past the analyser schedules so the counts are
        // comparable, and past zero so the noise block is eligible.
        let scenario = Scenario {
            time: 0.5,
            ..boundary.scenario.clone()
        };
        let snapshot = scenario.force(&mut oracle);
        let count = draws_for(&mut oracle, &scenario);
        if snapshot.tripped {
            tripped_counts.push(count);
        } else {
            healthy_counts.push(count);
        }
    }

    println!("tripped: {tripped_counts:?}");
    println!("healthy: {healthy_counts:?}");
    assert!(!tripped_counts.is_empty() && !healthy_counts.is_empty());
    let max_tripped = *tripped_counts.iter().max().expect("non-empty");
    let min_healthy = *healthy_counts.iter().min().expect("non-empty");
    assert_eq!(
        min_healthy - max_tripped,
        22 * 12,
        "a tripped evaluation drew {max_tripped} and a healthy one \
         {min_healthy}; the difference should be exactly the 264-draw noise \
         block"
    );

    // The part a port would get wrong: a trip does *not* silence the
    // generator. `teprob.f:711` guards only the continuous-measurement noise.
    // The three analyser blocks at `744-761` have no `ISD` guard at all, so
    // they keep drawing, and at t = 0.5 a tripped evaluation still makes
    // 30 walk + 168 gas + 60 product = 258 draws.
    assert_eq!(
        max_tripped,
        30 + 14 * 12 + 5 * 12,
        "a tripped plant should still advance the walks and sample both \
         analysers; only teprob.f:711-716 is guarded"
    );
    assert_ne!(max_tripped, 0, "a trip must not silence the generator");
}

/// The draw count is not constant across a run, so a port cannot get the
/// stream right by counting.
///
/// This is the finding that justifies Tier 3 having a *trace* rather than a
/// count: the walk channels fire on twelve independent schedules, and the
/// analysers on two more.
#[test]
fn the_draw_count_varies_across_the_trajectory() {
    let mut oracle = Oracle::lock();
    let pools = Pools::collect(&mut oracle, 400, DT);
    let mut counts = std::collections::BTreeMap::new();

    for index in 0..pools.trajectory.len() {
        let scenario = pools.nominal_case(index);
        *counts
            .entry(draws_for(&mut oracle, &scenario))
            .or_insert(0_usize) += 1;
    }

    println!("draw counts over 400 nominal states: {counts:?}");
    assert!(
        counts.len() > 1,
        "every evaluation drew the same {} times, so either the walks never \
         fire in this window or the counting is wrong",
        counts.keys().next().expect("non-empty")
    );
}
