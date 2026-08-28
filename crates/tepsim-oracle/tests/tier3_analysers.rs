//! Tier 3 for the noise and the analysers: `tepsim_core::analysers` against
//! `teprob.f:711-761`, comparing all forty-one measurements *and* the draw
//! trace of a whole step.
//!
//! # The whole step, at last
//!
//! Every earlier file drives one piece of `TEFUNC` and compares it. This one
//! can drive the entire thing: `Plant::advance_discrete`, `derivatives` and
//! `sample_measurements` together are all of `teprob.f:340-816`, so the port's
//! generator must end each step in exactly the place the Fortran's does.
//!
//! That is the strongest statement the project can make short of Tier 4, and
//! it is a statement no individual tier makes: Tier 1 pins the generator per
//! call, Tier 2 pins it per evaluation, and only here does it run free.
//!
//! # Why the guards matter more than the arithmetic
//!
//! `teprob.f:711` skips the continuous noise at `TIME = 0` and on a tripped
//! plant, while the analyser blocks at `744-761` have no such guard. B-0027
//! measured a tripped evaluation at 258 draws against a healthy one at 522. A
//! port that silenced everything on a trip would leave the generator 264 steps
//! behind and every later draw would differ, while producing measurements that
//! look entirely reasonable.

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

use tepsim_core::walk::Walks;
use tepsim_core::{Inputs, Plant, Segment, SimTime, State, vessels};
use tepsim_oracle::tier2::{Pools, Scenario};
use tepsim_oracle::{Oracle, Wlk, tier3};

const DT: f64 = 1.0 / 3600.0;

/// `PLAN.org`, "Tier 2". The measurements are compared against the same gate.
const TOLERANCE: f64 = 1e-12;

fn walks_of(wlk: &Wlk) -> Walks {
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

/// Put the port into a scenario's condition, analyser latch included.
fn configured(scenario: &Scenario) -> (Plant, Inputs) {
    let mut plant = Plant::new();
    plant.set_channels(walks_of(&scenario.walk));
    plant.set_rng(scenario.rng);
    plant.set_valve_command(scenario.common.vcv);
    plant.set_seeds(vessels::TemperatureSeeds {
        reactor: scenario.common.tcr,
        separator: scenario.common.tcs,
        stripper: scenario.common.tcc,
        mixing: scenario.common.tcv,
    });
    plant.set_analysers(tepsim_core::Analysers {
        // `XDEL(23..41)` is stored from index 22 in the Fortran's 41-long
        // array; the port's is 19 long and starts at `XMEAS(23)`.
        stored: core::array::from_fn(|i| scenario.common.xdel[22 + i]),
        // `XMEAS(23..41)` as the previous step left it. The Fortran holds
        // these between samples, so a harness that did not restore them would
        // be comparing against whatever the last unrelated call left behind.
        reported: core::array::from_fn(|i| scenario.measurements[22 + i]),
        next_gas: scenario.common.tgas,
        next_product: scenario.common.tprod,
    });
    let inputs = Inputs {
        manipulated: scenario.manipulated,
        disturbances: core::array::from_fn(|i| f64::from(scenario.disturbances[i])),
    };
    (plant, inputs)
}

/// Run one whole step on both sides and compare measurements, analyser state
/// and the entire draw trace.
fn compare(oracle: &mut Oracle, scenario: &Scenario, label: &str) -> (usize, f64) {
    // The port: all three phases, in the order `euler_step` runs them.
    let (mut plant, inputs) = configured(scenario);
    let state = State::from_flat(&scenario.state);
    let t = SimTime(scenario.time);
    plant.advance_discrete(t, &inputs);
    let (_, signals) = plant.derivatives(t, &state, &inputs).expect("converges");
    let ours = plant.sample_measurements(t, &signals);

    // The Fortran.
    oracle.set_teproc(&scenario.common);
    oracle.set_wlk(&scenario.walk);
    oracle.set_rng(scenario.rng);
    oracle.set_manipulated(&scenario.manipulated);
    oracle.set_disturbances(&scenario.disturbances);
    tier3::clear(oracle);
    let _ = oracle.derivatives(scenario.time, &scenario.state);
    let theirs = oracle.measurements();
    let after = oracle.teproc();
    let draws = tier3::trace(oracle).len();

    // All forty-one measurements.
    let mut worst = 0.0_f64;
    for (index, (a, b)) in ours.as_array().iter().zip(theirs).enumerate() {
        if b == 0.0 {
            assert_eq!(
                a.to_bits(),
                b.to_bits(),
                "{label}: XMEAS({}) is {a}, the Fortran has exactly zero",
                index + 1
            );
            continue;
        }
        let relative = (a - b).abs() / b.abs();
        worst = worst.max(relative);
        assert!(
            relative < TOLERANCE,
            "{label}: XMEAS({}) is {a}, the Fortran has {b} ({relative:e})",
            index + 1
        );
    }

    // The analyser latch and both schedules.
    //
    // `XDEL` is a composition, so under the vendored libm it carries the same
    // one-ULP `exp` difference everything downstream of the equilibrium does;
    // it is compared against the Tier 2 gate rather than bit for bit. The two
    // *schedules* are pure arithmetic on `0.1` and `0.25` and must be exact.
    for index in 0..19 {
        let (ours, theirs) = (plant.analysers().stored[index], after.xdel[22 + index]);
        if theirs == 0.0 {
            assert_eq!(
                ours.to_bits(),
                theirs.to_bits(),
                "{label}: XDEL({})",
                index + 23
            );
        } else {
            let relative = (ours - theirs).abs() / theirs.abs();
            assert!(
                relative < TOLERANCE,
                "{label}: XDEL({}) is {ours}, the Fortran has {theirs} \
                 ({relative:e})",
                index + 23
            );
        }
    }
    assert_eq!(
        plant.analysers().next_gas.to_bits(),
        after.tgas.to_bits(),
        "{label}: TGAS"
    );
    assert_eq!(
        plant.analysers().next_product.to_bits(),
        after.tprod.to_bits(),
        "{label}: TPROD"
    );

    // And the generator, after a whole step. Nothing else in the project
    // asserts this: it is only true once every consumer is ported.
    assert_eq!(
        plant.rng().to_bits(),
        oracle.rng().to_bits(),
        "{label}: after one whole step the port's generator is at {} and the \
         Fortran's at {}. Some consumer draws a different number of times.",
        plant.rng(),
        oracle.rng()
    );
    (draws, worst)
}

#[test]
fn a_whole_step_matches_the_fortran_along_the_nominal_trajectory() {
    let mut oracle = Oracle::lock();
    let pools = Pools::collect(&mut oracle, 400, DT);
    let mut total = 0;
    let mut worst = 0.0_f64;

    for index in 0..pools.trajectory.len() {
        let (draws, relative) = compare(
            &mut oracle,
            &pools.nominal_case(index),
            &format!("nominal#{index}"),
        );
        total += draws;
        worst = worst.max(relative);
    }
    println!("{total} draws over 400 whole steps; worst measurement {worst:e}");
    assert!(total > 100_000, "only {total} draws");
}

/// Every disturbance, over a run long enough for both analyser schedules to
/// come due many times.
#[test]
fn a_whole_step_matches_with_every_disturbance_active() {
    let mut oracle = Oracle::lock();
    let pools = Pools::collect(&mut oracle, 10, DT);
    let base = pools.nominal_case(0);
    let mut worst = 0.0_f64;

    for fault in 1..=20 {
        let mut scenario = base.clone();
        scenario.disturbances[fault - 1] = 1;
        let mut t = 0.0;
        for step in 0..120 {
            scenario.time = t;
            let (_, relative) =
                compare(&mut oracle, &scenario, &format!("IDV({fault}) step {step}"));
            worst = worst.max(relative);
            scenario.walk = oracle.wlk();
            scenario.rng = oracle.rng();
            scenario.common = oracle.teproc();
            // The analysers hold `XMEAS(23..41)` between samples, so this is
            // part of the state and has to be carried like the rest.
            scenario.measurements = oracle.measurements();
            t += 0.05;
        }
    }
    println!("all twenty disturbances, 6 simulated hours each; worst {worst:e}");
}

/// A tripped plant draws for the analysers and not for the continuous noise.
///
/// The guard at `teprob.f:711` is on `ISD`; the blocks at `744-761` are not.
/// A port that silenced everything would desynchronise the stream while its
/// measurements still looked reasonable.
#[test]
fn a_trip_silences_the_continuous_noise_and_not_the_analysers() {
    let mut oracle = Oracle::lock();
    let pools = Pools::collect(&mut oracle, 60, DT);
    let (boundaries, _) =
        tepsim_oracle::tier2::adversarial::build(&mut oracle, &pools.nominal_case(0));

    let mut tripped = Vec::new();
    let mut healthy = Vec::new();
    for boundary in &boundaries {
        let scenario = Scenario {
            time: 0.5,
            ..boundary.scenario.clone()
        };
        let snapshot = scenario.force(&mut oracle);
        let (draws, _) = compare(&mut oracle, &scenario, boundary.target.name);
        if snapshot.tripped {
            tripped.push(draws);
        } else {
            healthy.push(draws);
        }
    }

    let max_tripped = *tripped.iter().max().expect("some state trips");
    let min_healthy = *healthy.iter().min().expect("some state does not");
    println!("tripped {max_tripped} draws, healthy {min_healthy}");
    assert_eq!(
        min_healthy - max_tripped,
        22 * 12,
        "the difference should be exactly the 264-draw continuous noise block"
    );
    assert_eq!(
        max_tripped,
        30 + 14 * 12 + 5 * 12,
        "a tripped plant should still advance the walks and sample both \
         analysers"
    );
}

/// Many steps in a row, carrying the port's own state forward rather than
/// re-forcing it from the Fortran each time.
///
/// Every other test in the project restores the port from `COMMON` before each
/// comparison, which hides an error that only shows once the port's own state
/// has been carried. This one carries it.
#[test]
fn the_port_stays_in_step_when_it_carries_its_own_state() {
    let mut oracle = Oracle::lock();
    let pools = Pools::collect(&mut oracle, 10, DT);
    let scenario = pools.nominal_case(0);

    let (mut plant, inputs) = configured(&scenario);
    oracle.set_teproc(&scenario.common);
    oracle.set_wlk(&scenario.walk);
    oracle.set_rng(scenario.rng);
    oracle.set_manipulated(&scenario.manipulated);
    oracle.set_disturbances(&scenario.disturbances);

    let mut state = State::from_flat(&scenario.state);
    let mut fortran = scenario.state;
    let mut t = 0.0;

    for step in 0..2_000 {
        let time = SimTime(t);
        plant.advance_discrete(time, &inputs);
        let (derivative, signals) = plant.derivatives(time, &state, &inputs).expect("converges");
        let ours = plant.sample_measurements(time, &signals);

        let yp = oracle.derivatives(t, &fortran);
        let theirs = oracle.measurements();

        assert_eq!(
            plant.rng().to_bits(),
            oracle.rng().to_bits(),
            "step {step} at t={t}: the generators parted company"
        );
        for (index, (a, b)) in ours.as_array().iter().zip(theirs).enumerate() {
            if b == 0.0 {
                continue;
            }
            let relative = (a - b).abs() / b.abs();
            assert!(
                relative < 1e-9,
                "step {step}: XMEAS({}) is {a}, the Fortran has {b}",
                index + 1
            );
        }

        state = state.step(DT, &derivative);
        for (slot, rate) in fortran.iter_mut().zip(yp) {
            *slot += DT * rate;
        }
        t += DT;
    }
    println!("2,000 steps carried on both sides, generators in lockstep throughout");
}
