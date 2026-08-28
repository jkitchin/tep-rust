//! The facade's own properties, without the oracle.
//!
//! What the port computes is settled elsewhere: `facade_equivalence.rs` in
//! `tepsim-oracle` proves this loop is bit-identical to the one validated
//! against the Fortran. These are about the *API*: that a scenario determines
//! its run, that the builders build what they say, that a trip is reported
//! rather than thrown, and that ground truth is recorded.

// Exact comparisons throughout are the point: a run is a pure function of its
// scenario, so equality is exact or the property does not hold.
#![allow(
    clippy::float_cmp,
    reason = "determinism means the values are exactly equal or the test fails"
)]

use tepsim::{Outcome, Scenario, Simulation, forced_disturbance_step};

/// Short enough to be quick, long enough to cross an analyser sample and
/// several control periods.
const SHORT: f64 = 0.5;

// ---------------------------------------------------------------------------
// Determinism
// ---------------------------------------------------------------------------

/// A run is a pure function of its scenario.
///
/// This is the property everything downstream rests on: a dataset that can be
/// regenerated from its description does not have to be stored, shipped or
/// trusted.
#[test]
fn the_same_scenario_gives_the_same_run() {
    let scenario = Scenario::fault(6).with_hours(SHORT);
    let first = Simulation::new(scenario).run();
    let second = Simulation::new(scenario).run();
    assert_eq!(first, second);

    // Interleaving another simulation changes nothing: there is no shared
    // state. The original cannot say this, because its whole working set lives
    // in six COMMON blocks.
    let mut a = Simulation::new(scenario);
    let mut b = Simulation::new(Scenario::fault(13).with_hours(SHORT));
    let mut interleaved = Vec::new();
    loop {
        let left = a.step();
        let right = b.step();
        if let Some(sample) = left {
            interleaved.push(sample);
        }
        if right.is_none() && a.is_halted() {
            break;
        }
    }
    assert_eq!(interleaved, first.samples);
}

#[test]
fn different_seeds_give_different_runs() {
    let base = Scenario::baseline().with_hours(SHORT);
    let first = Simulation::new(base.with_seed(4_651_207_995.0)).run();
    let second = Simulation::new(base.with_seed(1_431_655_765.0)).run();
    assert_ne!(first.samples, second.samples);

    // A real difference, not a last-bit one: the disturbance walks are
    // different realisations.
    let worst = first
        .samples
        .iter()
        .zip(&second.samples)
        .flat_map(|(a, b)| a.row().into_iter().zip(b.row()))
        .filter(|(_, y)| *y != 0.0)
        .map(|(x, y)| (x - y).abs() / y.abs())
        .fold(0.0_f64, f64::max);
    println!("two seeds differ by {worst:.3e} at worst");
    assert!(worst > 1e-6, "the seed barely mattered: {worst:.3e}");
}

// ---------------------------------------------------------------------------
// The scenario builders
// ---------------------------------------------------------------------------

#[test]
fn the_baseline_is_the_run_the_original_driver_does() {
    let scenario = Scenario::baseline();
    assert_eq!(scenario.hours, 48.0);
    // NPTS = 172800 at a one-second step.
    assert_eq!(scenario.steps(), 172_800);
    // Sampled every 180 steps, as temain_mod.f:401 does.
    assert_eq!(scenario.samples(), 960);
    assert!(scenario.controlled);
    assert_eq!(scenario.active_faults().count(), 0);
    // The seed compiled into teprob.f:1187.
    assert_eq!(scenario.seed, 4_651_207_995.0);
}

#[test]
fn the_builders_build_what_they_say() {
    let scenario = Scenario::fault(7)
        .with_hours(3.0)
        .with_seed(99.0)
        .sampling_every(60)
        .open_loop()
        .with_fault(12);

    assert_eq!(scenario.hours, 3.0);
    assert_eq!(scenario.seed, 99.0);
    assert_eq!(scenario.sample_every, 60);
    assert!(!scenario.controlled);
    assert_eq!(scenario.active_faults().collect::<Vec<_>>(), vec![7, 12]);
    assert_eq!(scenario.steps(), 10_800);
    assert_eq!(scenario.samples(), 180);

    let idv = scenario.disturbance_vector();
    assert_eq!(idv[6], 1.0);
    assert_eq!(idv[11], 1.0);
    assert_eq!(idv.iter().sum::<f64>(), 2.0);
}

/// The step count comes out exact rather than one short.
///
/// `hours / step_hours` is not always exact, and truncating the quotient then
/// loses a step. Whole hours happen to be safe at the default one-second step:
/// `48.0 / (1.0 / 3600.0)` is exactly 172800. **4.1 hours is not**, it is
/// 14759.99999999999818, and truncation gives 14,759.
///
/// So the rounding is load bearing but not for the reason one would guess, and
/// the case below is the one that actually discriminates. A silent off-by-one
/// would put a sample count one out and break a comparison against a published
/// file for no visible reason.
#[test]
fn the_step_count_does_not_lose_a_step_to_rounding() {
    for hours in [1.0, 2.0, 8.0, 24.0, 48.0, 100.0] {
        let scenario = Scenario::baseline().with_hours(hours);
        assert_eq!(
            scenario.steps(),
            (hours * 3600.0) as usize,
            "{hours} h came out at {} steps",
            scenario.steps()
        );
    }

    // The discriminating case: truncation really would lose a step here.
    let quotient = 4.1_f64 / (1.0 / 3600.0);
    assert_eq!(quotient as usize, 14_759, "4.1 h no longer truncates short");
    assert_eq!(Scenario::baseline().with_hours(4.1).steps(), 14_760);
}

#[test]
#[should_panic(expected = "out of range")]
fn a_fault_index_out_of_range_is_rejected() {
    let _ = Scenario::fault(21);
}

// ---------------------------------------------------------------------------
// Outcomes
// ---------------------------------------------------------------------------

#[test]
fn the_controlled_plant_survives_and_the_open_loop_one_does_not() {
    let closed = Simulation::new(Scenario::baseline().with_hours(4.0)).run();
    assert_eq!(closed.outcome, Outcome::Completed);
    assert_eq!(closed.tripped_at(), None);

    let open = Simulation::new(Scenario::baseline().with_hours(4.0).open_loop()).run();
    let Outcome::Tripped { step, hours, cause } = open.outcome else {
        panic!("the open-loop plant did not trip: {:?}", open.outcome);
    };
    println!("open loop tripped at step {step} ({hours:.3} h) on {cause:?}");
    assert!(hours < 4.0);
}

/// A trip freezes the plant rather than stopping the run, and the frozen
/// samples keep coming.
///
/// `teprob.f:807-811` zeroes all fifty derivatives. Discarding what follows
/// would hide the difference between a port that trips where the original does
/// and one that does not, and every published dataset of a tripped run contains
/// exactly these constant rows.
#[test]
fn a_trip_freezes_the_plant_rather_than_ending_the_run() {
    let scenario = Scenario::baseline().with_hours(4.0).open_loop();
    let run = Simulation::new(scenario).run();
    let trip = run.tripped_at().expect("it trips");

    assert_eq!(
        run.samples.len(),
        scenario.samples(),
        "the run stopped early instead of freezing"
    );

    // After the trip the state is frozen, so the derivative is zero and the
    // measurements stop moving except for the noise added afterwards.
    let after: Vec<_> = run
        .samples
        .iter()
        .filter(|s| s.step > trip + scenario.sample_every)
        .collect();
    assert!(
        after.len() > 2,
        "not enough samples after the trip to check"
    );
    let pressure: Vec<f64> = after.iter().map(|s| s.measurements[6]).collect();
    let spread = pressure.iter().fold(f64::NEG_INFINITY, |a, b| a.max(*b))
        - pressure.iter().fold(f64::INFINITY, |a, b| a.min(*b));
    println!(
        "{} samples after the trip, reactor pressure spread {spread:.6}",
        after.len()
    );
    assert!(
        spread < 1.0,
        "the plant kept moving after the trip: pressure spread {spread}"
    );
}

/// And with the Class C fix on, it stops instead. Delta D-007.
#[test]
fn the_quirk_fix_ends_the_run_at_a_trip() {
    let mut scenario = Scenario::baseline().with_hours(4.0).open_loop();
    scenario.quirks.trip_ends_the_run = true;
    let run = Simulation::new(scenario).run();

    assert!(matches!(run.outcome, Outcome::Tripped { .. }));
    assert!(
        run.samples.len() < scenario.samples(),
        "the run continued past the trip with the fix on"
    );
    println!(
        "with trip_ends_the_run: {} rows instead of {}",
        run.samples.len(),
        scenario.samples()
    );
}

// ---------------------------------------------------------------------------
// Ground truth
// ---------------------------------------------------------------------------

#[test]
fn a_requested_fault_is_labelled_from_the_first_sample() {
    let run = Simulation::new(Scenario::fault(4).with_hours(SHORT)).run();
    let first = &run.samples[0];
    assert!(first.labels.faulted());
    assert_eq!(first.labels.faults().collect::<Vec<_>>(), vec![4]);
    // Onset at the start, so the elapsed time is the sample's own time.
    let since = first.labels.since_onset[3].expect("IDV(4) has an onset");
    assert!((since - first.hours).abs() < 1e-12);

    // And the fault-free plant is labelled fault-free.
    let clean = Simulation::new(Scenario::baseline().with_hours(SHORT)).run();
    assert!(!clean.samples[0].labels.faulted());
    assert_eq!(clean.samples[0].labels.since_onset[3], None);
}

/// The driver's unasked-for `IDV(12)` is labelled, with its own onset time.
///
/// This is the whole argument for recording ground truth rather than assuming
/// it. A run labelled "fault-free" in the literature is fault-free for eight
/// hours and then is not, and nothing in the published data says so. Delta
/// D-011.
#[test]
fn the_drivers_forced_disturbance_is_labelled_with_its_own_onset() {
    let scenario = Scenario::baseline().with_hours(9.0);
    let run = Simulation::new(scenario).run();

    let before = run
        .samples
        .iter()
        .find(|s| s.step == forced_disturbance_step() - 180)
        .expect("a sample before the mark");
    let after = run
        .samples
        .iter()
        .find(|s| s.step == forced_disturbance_step() + 180)
        .expect("a sample after the mark");

    assert!(!before.labels.faulted(), "labelled faulted before the mark");
    assert_eq!(after.labels.faults().collect::<Vec<_>>(), vec![12]);
    let since = after.labels.since_onset[11].expect("IDV(12) has an onset");
    println!(
        "IDV(12) onset at step {}, {since:.4} h before sample at step {}",
        forced_disturbance_step(),
        after.step
    );
    assert!(since > 0.0 && since < 0.2);

    // Turning the quirk off leaves the run genuinely fault-free.
    let mut honest = scenario;
    honest.driver_forces_idv12 = false;
    let clean = Simulation::new(honest).run();
    assert!(
        clean.samples.iter().all(|s| !s.labels.faulted()),
        "IDV(12) still fired with the forcing off"
    );
    assert_ne!(
        clean.samples.last().map(|s| s.row()),
        run.samples.last().map(|s| s.row()),
        "the forced disturbance changed nothing, so this test proves nothing"
    );
}

// ---------------------------------------------------------------------------
// The output shape
// ---------------------------------------------------------------------------

#[test]
fn the_columns_are_the_fifty_three_channels_in_order() {
    let run = Simulation::new(Scenario::baseline().with_hours(SHORT)).run();
    assert_eq!(run.columns().len(), tepsim::CHANNELS);
    assert_eq!(tepsim::channel_names().len(), tepsim::CHANNELS);

    // Column 6 is XMEAS(7), reactor pressure. Column 41 is XMV(1).
    let pressure = run.measurement(7);
    assert_eq!(pressure, run.column(6));
    assert_eq!(run.manipulated(1), run.column(41));

    // Every column is as long as the run.
    for channel in 0..tepsim::CHANNELS {
        assert_eq!(run.column(channel).len(), run.samples.len());
    }

    // Reactor pressure sits where Downs and Vogel say it does.
    assert!(
        pressure.iter().all(|p| (2600.0..2800.0).contains(p)),
        "reactor pressure left its nominal band"
    );
}

#[test]
fn sampling_every_step_records_every_step() {
    let scenario = Scenario::baseline().with_hours(0.01).sampling_every(1);
    let run = Simulation::new(scenario).run();
    assert_eq!(run.samples.len(), 36);
    for (index, sample) in run.samples.iter().enumerate() {
        assert_eq!(sample.step, index + 1);
    }
}
