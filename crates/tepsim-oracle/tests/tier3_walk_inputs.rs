//! `Plant::advance_discrete` against `teprob.f:340-416` and `793-804`.
//!
//! # This item exists because six copies agreed and were all wrong
//!
//! Every Tier 2 differential from B-0020 onward reimplemented
//! `teprob.f:407-416` in its own harness, to supply the port's inputs. All six
//! came out 0 ULP, which looked like six independent confirmations of the
//! reading.
//!
//! They were one reading copied six times, and it was incomplete.
//! `teprob.f:407-408` is
//!
//! ```fortran
//!       XST(1,4)=TESUB8(1,TIME)-IDV(1)*0.03D0
//!      .-IDV(2)*2.43719D-3
//! ```
//!
//! and every copy carried the `IDV(1)` term and dropped the `IDV(2)` one. It
//! was invisible because `Pools::collect` starts from `TEINIT`, where every
//! `IDV` is zero, and neither the perturbed nor the adversarial pool switches
//! one on. So the term was multiplied by zero in every one of the ~2,400
//! states each of those files sweeps.
//!
//! The lesson is not about that term. It is that a comparison can be exact
//! over thousands of states and still say nothing about a branch none of them
//! enter, and that copying a harness copies its blind spots. This file
//! therefore sweeps with disturbances *on*, one at a time and in combination.

#![cfg(feature = "oracle")]

use tepsim_core::walk::{CHANNELS, Walks};
use tepsim_core::{Inputs, Plant, Segment, SimTime, State};
use tepsim_oracle::tier2::{Pools, Scenario};
use tepsim_oracle::{Oracle, Wlk};

const DT: f64 = 1.0 / 3600.0;

/// Put the port into the walk state a scenario carries.
fn channels_from(wlk: &Wlk) -> Walks {
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

/// Run one step on both sides from the same condition and compare every
/// walk-driven quantity the pure phase reads.
fn compare(oracle: &mut Oracle, scenario: &Scenario, label: &str) {
    let mut plant = Plant::new();
    plant.set_channels(channels_from(&scenario.walk));
    plant.set_rng(scenario.rng);
    plant.set_valve_command(scenario.common.vcv);
    let inputs = Inputs {
        manipulated: scenario.manipulated,
        disturbances: core::array::from_fn(|i| f64::from(scenario.disturbances[i])),
    };
    plant.advance_discrete(SimTime(scenario.time), &inputs);

    let snapshot = scenario.force(oracle);
    let after = oracle.wlk();
    let w = plant.walk_inputs();

    // teprob.f:407-410, the mixed feed composition.
    let pairs = [
        (
            "XST(1,4)",
            w.feed.ac_feed_light[0],
            snapshot.common.xst[3][0],
        ),
        (
            "XST(2,4)",
            w.feed.ac_feed_light[1],
            snapshot.common.xst[3][1],
        ),
        (
            "XST(3,4)",
            w.feed.ac_feed_light[2],
            snapshot.common.xst[3][2],
        ),
        ("TST(1)", w.feed.d_feed_celsius, snapshot.common.tst[0]),
        ("TST(4)", w.feed.ac_feed_celsius, snapshot.common.tst[3]),
        ("TCWR", w.coolant.reactor, snapshot.common.tcwr),
        ("TCWS", w.coolant.condenser, snapshot.common.tcws),
    ];
    for (name, ours, theirs) in pairs {
        assert_eq!(
            ours.to_bits(),
            theirs.to_bits(),
            "{label}: {name} is {ours}, the Fortran has {theirs}"
        );
    }

    // The four drifts are TEFUNC locals, so they are checked through the
    // channels they come from: if the walk state matches, `TESUB8` of it does.
    for index in 0..CHANNELS {
        let channel = &plant.channels().channels[index];
        let expected = [
            after.adist[index],
            after.bdist[index],
            after.cdist[index],
            after.ddist[index],
            after.tnext[index],
            after.tlast[index],
        ];
        let ours = [
            channel.segment.constant,
            channel.segment.linear,
            channel.segment.quadratic,
            channel.segment.cubic,
            channel.segment.until,
            channel.since,
        ];
        for (slot, (a, b)) in ours.iter().zip(expected).enumerate() {
            assert_eq!(
                a.to_bits(),
                b.to_bits(),
                "{label}: channel {} field {slot} is {a}, the Fortran has {b}",
                index + 1
            );
        }
    }

    // The generator is deliberately *not* compared here. A full `TEFUNC`
    // evaluation also draws 264 times for measurement noise and up to 228 more
    // for the analysers (`teprob.f:711-761`), and none of that is
    // `advance_discrete`'s: it is B-0024b's, in the post-phase. So the port's
    // word is legitimately behind the Fortran's after this call.
    //
    // The stream is covered where it can be: `tier3_walk.rs` compares the
    // walk's draws as a prefix of the evaluation's trace, which is the same
    // claim without the confound.

    // The valve latch, hoisted here from teprob.f:793-804.
    for (index, (ours, theirs)) in plant
        .valve_command()
        .iter()
        .zip(snapshot.common.vcv)
        .enumerate()
    {
        assert_eq!(
            ours.to_bits(),
            theirs.to_bits(),
            "{label}: VCV({}) is {ours}, the Fortran has {theirs}",
            index + 1
        );
    }
}

#[test]
fn the_walk_driven_inputs_match_along_the_nominal_trajectory() {
    let mut oracle = Oracle::lock();
    let pools = Pools::collect(&mut oracle, 400, DT);
    for index in 0..pools.trajectory.len() {
        compare(
            &mut oracle,
            &pools.nominal_case(index),
            &format!("nominal#{index}"),
        );
    }
}

/// Every disturbance, one at a time, over a long enough run to matter.
///
/// This is the sweep the six copied harnesses never ran, and the only one that
/// can see a term multiplied by an `IDV`.
#[test]
fn every_disturbance_matches_over_a_long_run() {
    let mut oracle = Oracle::lock();
    let pools = Pools::collect(&mut oracle, 10, DT);
    let base = pools.nominal_case(0);

    for fault in 1..=20 {
        let mut scenario = base.clone();
        scenario.disturbances[fault - 1] = 1;
        let mut t = 0.0;
        for step in 0..200 {
            scenario.time = t;
            compare(&mut oracle, &scenario, &format!("IDV({fault}) step {step}"));
            scenario.walk = oracle.wlk();
            scenario.rng = oracle.rng();
            scenario.common = oracle.teproc();
            t += 0.05;
        }
    }
    println!("all twenty disturbances matched over 10 simulated hours each");
}

/// The two composition step faults, specifically.
///
/// `IDV(1)` shifts A down; `IDV(2)` shifts A down *and* B up, through two
/// separate terms on two separate lines. The `IDV(2)` term on `XST(1,4)` is
/// the one every earlier harness dropped, so it gets its own check with a
/// stated magnitude rather than only appearing inside a sweep.
#[test]
fn the_two_composition_faults_move_different_amounts() {
    let mut oracle = Oracle::lock();
    let pools = Pools::collect(&mut oracle, 10, DT);
    let base = pools.nominal_case(0);

    let run = |oracle: &mut Oracle, fault: Option<usize>| -> [f64; 3] {
        let mut scenario = base.clone();
        if let Some(f) = fault {
            scenario.disturbances[f - 1] = 1;
        }
        let snapshot = scenario.force(oracle);
        [
            snapshot.common.xst[3][0],
            snapshot.common.xst[3][1],
            snapshot.common.xst[3][2],
        ]
    };

    let nominal = run(&mut oracle, None);
    let one = run(&mut oracle, Some(1));
    let two = run(&mut oracle, Some(2));

    println!("nominal {nominal:?}");
    println!("IDV(1)  {one:?}");
    println!("IDV(2)  {two:?}");

    // IDV(1): A down by 0.03, B untouched.
    assert!((nominal[0] - one[0] - 0.03).abs() < 1e-12, "IDV(1) on A");
    assert!(
        (nominal[1] - one[1]).abs() < 1e-15,
        "IDV(1) must not touch B"
    );

    // IDV(2): A down by 2.43719e-3 and B up by 0.005. Dropping the A term is
    // exactly the bug the six copied harnesses carried.
    assert!(
        (nominal[0] - two[0] - 2.43719e-3).abs() < 1e-12,
        "IDV(2) should move A by 2.43719e-3, moved it by {}",
        nominal[0] - two[0]
    );
    assert!((two[1] - nominal[1] - 0.005).abs() < 1e-12, "IDV(2) on B");

    // And C absorbs both, since teprob.f:410 makes it the remainder.
    for case in [one, two] {
        let sum = case[0] + case[1] + case[2];
        assert!(
            (sum - 1.0).abs() < 1e-15,
            "the three should sum to one: {sum}"
        );
    }
}

/// The generator moves exactly once per step, however many times the
/// derivative is evaluated.
///
/// That is the promise the three-phase split makes, and it is what lets an
/// RK4 driver exist at all. With the generator now living in `Plant`, it is
/// checkable directly.
#[test]
fn the_generator_moves_once_per_step_not_once_per_evaluation() {
    let mut oracle = Oracle::lock();
    let pools = Pools::collect(&mut oracle, 10, DT);
    let scenario = pools.nominal_case(0);

    let mut plant = Plant::new();
    plant.set_channels(channels_from(&scenario.walk));
    plant.set_rng(scenario.rng);
    let inputs = Inputs {
        manipulated: scenario.manipulated,
        disturbances: [0.0; 20],
    };
    let state = State::from_flat(&scenario.state);

    plant.advance_discrete(SimTime(0.15), &inputs);
    let after_advance = plant.rng();

    // Four evaluations, as RK4 would make.
    for _ in 0..4 {
        let _ = plant
            .derivatives(SimTime(0.15), &state, &inputs)
            .expect("converges");
    }
    assert_eq!(
        plant.rng().to_bits(),
        after_advance.to_bits(),
        "evaluating the derivative moved the generator, so the pure phase is \
         not pure and an RK4 driver would draw four times per step"
    );
}
