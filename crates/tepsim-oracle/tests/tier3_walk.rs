//! Tier 3 for the walk advance: `tepsim_core::walk` against
//! `teprob.f:340-406`, comparing both the walk *state* and the *draw trace*.
//!
//! # Two claims, not one
//!
//! B-0028 demonstrated that a port can produce identical values from a
//! different stream, and B-0029 built the instrument that sees it. So this
//! file checks both:
//!
//! - **State.** After each advance, every channel's four cubic coefficients
//!   and its two times must match `COMMON/WLK/` bit for bit.
//! - **Trace.** The draws made getting there must match, in order, with the
//!   same scaling.
//!
//! Either alone would pass a port the other catches.
//!
//! # Driving the comparison
//!
//! The oracle's walk advance is buried inside `TEFUNC`, so it cannot be called
//! on its own. Instead both sides are put into the same walk state, `TEFUNC`
//! is evaluated once, and `COMMON/WLK/` is read back. The port's `advance` is
//! run against the same inputs and the two states compared.
//!
//! That means every comparison also drags the whole plant model through, and
//! its draws land in the same trace. The measurement noise is 264 of them and
//! it comes *after* the walk advance in the stream, so the walk's draws are a
//! prefix: comparing the first `n` is comparing the walk alone.

#![cfg(feature = "oracle")]

use tepsim_core::disturbance::Draw;
use tepsim_core::walk::{CHANNELS, Walks, advance};
use tepsim_core::{Segment, TepRng, TracingRng};
use tepsim_oracle::tier2::{Pools, Scenario};
use tepsim_oracle::{Oracle, Wlk, tier3};

const DT: f64 = 1.0 / 3600.0;

/// Read the twelve channels out of `COMMON/WLK/`.
fn channels_of(wlk: &Wlk) -> [(Segment, f64); CHANNELS] {
    core::array::from_fn(|i| {
        (
            Segment {
                constant: wlk.adist[i],
                linear: wlk.bdist[i],
                quadratic: wlk.cdist[i],
                cubic: wlk.ddist[i],
                until: wlk.tnext[i],
            },
            wlk.tlast[i],
        )
    })
}

/// Put the port into the same walk state the scenario carries.
fn walks_from(wlk: &Wlk) -> Walks {
    let mut walks = Walks::default();
    for (index, channel) in walks.channels.iter_mut().enumerate() {
        channel.segment = Segment {
            constant: wlk.adist[index],
            linear: wlk.bdist[index],
            quadratic: wlk.cdist[index],
            cubic: wlk.ddist[index],
            until: wlk.tnext[index],
        };
        channel.since = wlk.tlast[index];
    }
    walks.flags = wlk.idvwlk;
    walks
}

/// Advance both sides once from the same state and compare state and trace.
///
/// Returns how many draws the walk made, so the caller can tell an advance
/// that did something from one that did not.
fn compare(oracle: &mut Oracle, scenario: &Scenario, label: &str) -> usize {
    let disturbances: [f64; 20] = core::array::from_fn(|i| f64::from(scenario.disturbances[i]));

    // The port first, so its trace is the walk's alone.
    let mut walks = walks_from(&scenario.walk);
    let mut rng = TracingRng::new(TepRng::new(scenario.rng));
    advance(&mut walks, &mut rng, scenario.time, &disturbances);
    let ours: Vec<Draw> = rng.draws().to_vec();

    // Then the Fortran, whose trace has the whole evaluation in it.
    oracle.set_teproc(&scenario.common);
    oracle.set_wlk(&scenario.walk);
    oracle.set_rng(scenario.rng);
    oracle.set_manipulated(&scenario.manipulated);
    oracle.set_disturbances(&scenario.disturbances);
    tier3::clear(oracle);
    let _ = oracle.derivatives(scenario.time, &scenario.state);
    let theirs = tier3::trace(oracle);
    let after = oracle.wlk();

    // The walk advance runs first in `TEFUNC`, so its draws are a prefix.
    assert!(
        theirs.len() >= ours.len(),
        "{label}: the port drew {} and the whole Fortran evaluation only {}",
        ours.len(),
        theirs.len()
    );
    if let Some(divergence) = tier3::diff(&ours, &theirs[..ours.len()]) {
        panic!("{label}: the walk's draws diverge. {divergence}");
    }

    // And the state it produced.
    let expected = channels_of(&after);
    for (index, channel) in walks.channels.iter().enumerate() {
        let (segment, since) = expected[index];
        let pairs = [
            ("ADIST", channel.segment.constant, segment.constant),
            ("BDIST", channel.segment.linear, segment.linear),
            ("CDIST", channel.segment.quadratic, segment.quadratic),
            ("DDIST", channel.segment.cubic, segment.cubic),
            ("TNEXT", channel.segment.until, segment.until),
            ("TLAST", channel.since, since),
        ];
        for (name, a, b) in pairs {
            assert_eq!(
                a.to_bits(),
                b.to_bits(),
                "{label}: {name}({}) is {a}, the Fortran has {b}",
                index + 1
            );
        }
    }
    assert_eq!(walks.flags, after.idvwlk, "{label}: IDVWLK");
    ours.len()
}

#[test]
fn the_walk_advance_matches_the_fortran_along_the_nominal_trajectory() {
    let mut oracle = Oracle::lock();
    let pools = Pools::collect(&mut oracle, 400, DT);
    let mut advanced = 0;
    let mut quiet = 0;

    for index in 0..pools.trajectory.len() {
        let scenario = pools.nominal_case(index);
        let drew = compare(&mut oracle, &scenario, &format!("nominal#{index}"));
        if drew > 0 {
            advanced += 1;
        } else {
            quiet += 1;
        }
    }
    println!("{advanced} evaluations advanced a channel, {quiet} did not");
    assert!(
        advanced > 0,
        "no channel ever came due, so nothing was tested"
    );
    assert!(
        quiet > 0,
        "every evaluation advanced, which cannot be right"
    );
}

/// Every disturbance that drives a channel, one at a time, over long enough
/// for the channel to re-segment many times.
///
/// This is where the spike rule's *fire* branch is reached: B-0027 measured
/// that all three spike channels dwell at the nominal point, so an active
/// disturbance is the only way in.
#[test]
fn every_walk_disturbance_matches_the_fortran_over_a_long_run() {
    let mut oracle = Oracle::lock();
    let pools = Pools::collect(&mut oracle, 10, DT);
    let base = pools.nominal_case(0);

    // The ten disturbances that drive a walk channel (`teprob.f:347-358`).
    for fault in [8, 9, 10, 11, 12, 13, 16, 17, 18, 20] {
        let mut scenario = base.clone();
        scenario.disturbances[fault - 1] = 1;

        // Walk the state forward by hand, re-forcing at each step, so the
        // channels re-segment many times.
        let mut t = 0.0;
        let mut drew_total = 0;
        for step in 0..600 {
            scenario.time = t;
            drew_total += compare(&mut oracle, &scenario, &format!("IDV({fault}) step {step}"));
            // Carry the Fortran's own walk state forward, so the chain
            // advances rather than restarting from the same point.
            scenario.walk = oracle.wlk();
            scenario.rng = oracle.rng();
            t += 0.05;
        }
        println!("IDV({fault}): {drew_total} walk draws over 30 hours");
        assert!(
            drew_total > 50,
            "IDV({fault}) made only {drew_total} walk draws in 30 hours"
        );
    }
}

/// The spike channels must actually fire somewhere in the sweep, or the
/// `teprob.f:381-385` branch is implemented and never evaluated.
#[test]
fn the_spike_branch_is_reached_by_an_active_disturbance() {
    let mut oracle = Oracle::lock();
    let pools = Pools::collect(&mut oracle, 10, DT);
    let mut scenario = pools.nominal_case(0);
    scenario.disturbances[16] = 1; // IDV(17), channel 10

    let mut fired = 0;
    let mut dwelled = 0;
    let mut t = 0.0;
    for _ in 0..1_500 {
        scenario.time = t;
        let mut walks = walks_from(&scenario.walk);
        let mut rng = TracingRng::new(TepRng::new(scenario.rng));
        advance(
            &mut walks,
            &mut rng,
            t,
            &core::array::from_fn(|i| f64::from(scenario.disturbances[i])),
        );

        let channel = &walks.channels[9];
        if channel.segment.cubic != 0.0 {
            fired += 1;
        } else if channel.segment.quadratic != 0.0 {
            dwelled += 1;
        }

        oracle.set_teproc(&scenario.common);
        oracle.set_wlk(&scenario.walk);
        oracle.set_rng(scenario.rng);
        oracle.set_disturbances(&scenario.disturbances);
        let _ = oracle.derivatives(t, &scenario.state);
        scenario.walk = oracle.wlk();
        scenario.rng = oracle.rng();
        t += 0.05;
    }
    println!("channel 10 over 75 hours: {fired} fired segments, {dwelled} dwells");
    assert!(
        fired > 0,
        "the spike branch at teprob.f:381-385 was never taken"
    );
    assert!(
        dwelled > 0,
        "the dwell branch at teprob.f:387-393 was never taken"
    );
}

/// The `TIME = 0` reset, against the Fortran.
///
/// The port runs the advance and *then* discards it, because `teprob.f:397`
/// comes after both loops. A port that returned early would have the right
/// state and the wrong stream.
#[test]
fn the_time_zero_reset_matches_including_its_draws() {
    let mut oracle = Oracle::lock();
    let pools = Pools::collect(&mut oracle, 10, DT);
    let mut scenario = pools.nominal_case(0);
    scenario.time = 0.0;

    // Make every channel due, so the advance has something to do before the
    // reset throws it away.
    for index in 0..CHANNELS {
        scenario.walk.tnext[index] = 0.0;
    }
    let drew = compare(&mut oracle, &scenario, "the t=0 reset");
    println!("the t=0 reset still drew {drew} times");
    assert!(
        drew > 0,
        "nothing drew, so the reset is short-circuiting the advance on both \
         sides and this test cannot tell the two implementations apart"
    );
}
