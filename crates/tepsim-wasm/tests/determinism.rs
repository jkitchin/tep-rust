//! Determinism and transport invariants for the browser bindings.
//!
//! Cross-platform bit-identical output is a hard invariant of this project
//! (`PLAN.org`, "Numerics and determinism"; Tier 9). The bindings are the
//! easiest place to lose it by accident, because they are the one layer that
//! touches a foreign runtime: an `f32` in a buffer type, a `js_sys::Math` call
//! resolving to the host's libm, a `Date` in a progress calculation. None of
//! those exist and these tests are what keeps it that way.
//!
//! # The pinned digest
//!
//! `run_digest_is_pinned` hard-codes the digest of a fixed run. That number was
//! not chosen; it was measured, then written down. If a change moves it, the
//! change altered the numbers a browser produces, and the correct response is
//! to find out why rather than to update the constant. Moving it is a logged
//! re-baseline, the same as a `gfortran` or toolchain change.
//!
//! These run on the host. Comparing a browser's `selfCheckDigest` against the
//! same constant is the wasm half of Tier 9;
//! `runner::tepsim_wasm_self_check_digest` exports it to any WebAssembly
//! runtime without needing `wasm-bindgen` glue, so that half can be automated
//! before the browser app exists.

use tepsim::run::Outcome;
use tepsim::{Integrator, Scenario, Simulation};
use tepsim_wasm::channels;
use tepsim_wasm::digest::Fnv1a64;
use tepsim_wasm::runner::{
    ConfigError, MAX_SAMPLES, MAX_STEPS, ROW_WIDTH, Runner, hours_since_onset, scenario_digest,
    self_check_digest, self_check_scenario, tepsim_wasm_self_check_digest, validate,
};

/// Measured, then written down. See the module docs before changing it.
///
/// FNV-1a 64 over the 20 rows of [`self_check_scenario`], `ROW_WIDTH` values
/// each, hashed as IEEE 754 bit patterns in emission order. Produced on
/// rustc 1.97.1, aarch64-apple-darwin, 2026-08-28, and verified identical from
/// a `wasm32-unknown-unknown` build of the same commit.
const PINNED_RUN_DIGEST: u64 = 0xc8a2_6889_992f_1719;

/// Bit-exact comparison of two `f64`s.
///
/// `assert_eq!` on floats trips `clippy::float_cmp`, and rightly: approximate
/// equality is nearly always what is meant. Here it is not. These are pinned
/// constants and a one-ULP difference in any of them is exactly the failure
/// worth catching.
fn assert_bits(actual: f64, expected: f64, what: &str) {
    assert_eq!(
        actual.to_bits(),
        expected.to_bits(),
        "{what}: {actual} is not bit-identical to {expected}"
    );
}

/// A short run that still exercises the controllers and the measurement layer.
fn short() -> Scenario {
    Scenario::baseline().with_hours(0.5)
}

fn runner(scenario: Scenario) -> Runner {
    Runner::new(scenario).expect("the scenario is valid")
}

/// The self-check scenario is spelled out here as well as in the crate, so an
/// edit to one has to be a deliberate edit to both. The crate's copy is what a
/// browser runs; this copy is what says what a browser is supposed to run.
#[test]
fn the_self_check_scenario_is_the_one_documented() {
    let scenario = self_check_scenario();
    assert_bits(scenario.seed, 4_651_207_995.0, "the teprob.f:1187 seed");
    assert_bits(scenario.hours, 1.0, "duration");
    assert_bits(scenario.step_hours, 1.0 / 3600.0, "the one-second step");
    assert_eq!(scenario.sample_every, 180, "the 180-second output cadence");
    assert_eq!(scenario.disturbances, [false; 20], "fault free");
    assert!(scenario.controlled, "closed loop");
    assert!(scenario.driver_forces_idv12, "the driver's IDV(12), D-011");
    assert!(
        !scenario.quirks.trip_ends_the_run,
        "D-007 off, as published"
    );
    assert_eq!(
        scenario.integrator,
        Integrator::Euler,
        "the faithful method"
    );
    assert!(scenario.integrator.is_faithful());

    let plan = validate(&scenario).expect("valid");
    assert_eq!(plan.steps, 3600);
    assert_eq!(plan.samples, 20);
}

#[test]
fn run_digest_is_pinned() {
    let mut run = runner(self_check_scenario());
    let values = run.run_to_end();

    assert_eq!(values.len(), 20 * ROW_WIDTH);
    assert_eq!(
        run.checksum(),
        PINNED_RUN_DIGEST,
        "the browser bindings no longer produce the run they used to. \
         Find out what changed before touching this constant."
    );
}

/// The glue-free export and the Rust function must agree, or a runtime checking
/// the module without `wasm-bindgen` would be checking something else.
#[test]
fn the_glue_free_export_returns_the_pinned_digest() {
    assert_eq!(self_check_digest(), PINNED_RUN_DIGEST);
    assert_eq!(tepsim_wasm_self_check_digest(), PINNED_RUN_DIGEST);
}

/// The single most important test here. These bindings must not be a second,
/// slightly different simulator: chunking is a transport concern and must not
/// touch a number.
#[test]
fn chunked_output_matches_the_facades_own_run() {
    let scenario = short();
    let expected = Simulation::new(scenario).run();

    let mut run = runner(scenario);
    let values = run.run_to_end();

    assert_eq!(values.len(), expected.samples.len() * ROW_WIDTH);
    for (i, sample) in expected.samples.iter().enumerate() {
        let row = &values[i * ROW_WIDTH..(i + 1) * ROW_WIDTH];
        assert_bits(row[0], sample.hours, &format!("sample {i} time"));
        for (c, value) in sample.row().iter().enumerate() {
            assert_bits(
                row[1 + c],
                *value,
                &format!("sample {i} channel {c} ({})", channels::column_ids()[1 + c]),
            );
        }
    }
}

/// The property the pinned digest stands in for: the output is a function of
/// the scenario, not of how the caller chose to collect it.
#[test]
fn chunk_size_changes_nothing() {
    let scenario = short();
    let whole = {
        let mut run = runner(scenario);
        let values = run.run_to_end();
        (values, run.checksum())
    };

    // Sizes that straddle the ends: one sample at a time, one that does not
    // divide the run, exactly the run, and more than the run.
    for chunk in [1, 3, 10, 11, 4096] {
        let mut run = runner(scenario);
        let mut values = Vec::new();
        while !run.is_finished() {
            let batch = run.step_chunk(chunk);
            assert!(
                !batch.is_empty(),
                "step_chunk({chunk}) returned nothing while the run was unfinished, \
                 which would spin a worker forever"
            );
            assert_eq!(
                batch.len() % ROW_WIDTH,
                0,
                "step_chunk({chunk}) returned a partial row"
            );
            values.extend_from_slice(&batch);
        }
        assert_eq!(values, whole.0, "chunk size {chunk} changed the samples");
        assert_eq!(
            run.checksum(),
            whole.1,
            "chunk size {chunk} changed the digest"
        );
    }
}

#[test]
fn stepping_past_the_end_yields_nothing_and_does_not_panic() {
    let mut run = runner(short());
    let _ = run.run_to_end();
    assert!(run.is_finished());
    for _ in 0..3 {
        assert!(run.step_chunk(100).is_empty());
    }
    assert_eq!(run.emitted_samples(), run.plan().samples);
    assert_eq!(run.outcome_name(), Some("completed"));
}

#[test]
fn asking_for_zero_samples_is_a_no_op() {
    let mut run = runner(short());
    let before = run.checksum();
    assert!(run.step_chunk(0).is_empty());
    assert_eq!(run.emitted_samples(), 0);
    assert_eq!(run.steps_taken(), 0);
    assert_eq!(run.checksum(), before);
}

/// A different seed must produce a different run. Trivial to state, and exactly
/// what breaks if the seed is dropped on the floor between the scenario and the
/// generator.
#[test]
fn the_seed_reaches_the_generator() {
    let mut a = runner(short());
    let mut b = runner(short().with_seed(1_234_567.0));
    assert_ne!(a.run_to_end(), b.run_to_end());
    assert_ne!(a.checksum(), b.checksum());
}

/// A disturbance must change the trajectory and must be reported as ground
/// truth. `IDV(1)` steps the mixed feed's A fraction down by 0.03.
#[test]
fn a_fault_changes_the_run_and_is_labelled() {
    let clean = short();
    let faulted = clean.with_fault(1);

    assert_ne!(
        scenario_digest(&clean),
        scenario_digest(&faulted),
        "the scenario digest must distinguish fault sets"
    );

    let mut a = runner(clean);
    let mut b = runner(faulted);
    assert_ne!(
        a.run_to_end(),
        b.run_to_end(),
        "IDV(1) did not reach the plant"
    );

    assert!(
        a.labels().faults().next().is_none(),
        "the clean run is clean"
    );
    assert_eq!(
        b.labels().faults().collect::<Vec<_>>(),
        vec![1],
        "the faulted run says so"
    );
    assert!(
        hours_since_onset(b.labels(), 1).is_some(),
        "IDV(1) has an onset time"
    );
    assert_eq!(
        hours_since_onset(b.labels(), 2),
        None,
        "IDV(2) never came on"
    );
    assert_eq!(hours_since_onset(b.labels(), 0), None, "out of range");
    assert_eq!(hours_since_onset(b.labels(), 21), None, "out of range");
}

/// The driver forces `IDV(12)` on at eight hours whatever the scenario asked
/// for. Delta D-011. A browser plotting ground truth has to see it happen, and
/// this is the only place these bindings surface it.
#[test]
fn the_drivers_forced_idv12_shows_up_in_the_labels() {
    let scenario = Scenario::baseline().with_hours(9.0);
    let mut run = runner(scenario);
    assert!(
        run.scenario().active_faults().next().is_none(),
        "asked for none"
    );

    let _ = run.run_to_end();
    assert_eq!(
        run.labels().faults().collect::<Vec<_>>(),
        vec![12],
        "the driver switched IDV(12) on unasked at hour eight"
    );
    let since = hours_since_onset(run.labels(), 12).expect("an onset");
    assert!(
        (0.0..1.5).contains(&since),
        "IDV(12) came on about an hour before the end, not {since} h ago"
    );
}

#[test]
fn faults_are_one_based_and_bounded() {
    let mut scenario = short();
    assert_eq!(
        tepsim_wasm::runner::set_fault(&mut scenario, 0, true),
        Err(ConfigError::FaultIndex)
    );
    assert_eq!(
        tepsim_wasm::runner::set_fault(&mut scenario, 21, true),
        Err(ConfigError::FaultIndex)
    );
    assert!(!tepsim_wasm::runner::fault(&scenario, 0));
    assert!(!tepsim_wasm::runner::fault(&scenario, 21));

    assert_eq!(
        tepsim_wasm::runner::set_fault(&mut scenario, 1, true),
        Ok(())
    );
    assert_eq!(
        tepsim_wasm::runner::set_fault(&mut scenario, 20, true),
        Ok(())
    );
    assert!(tepsim_wasm::runner::fault(&scenario, 1));
    assert!(tepsim_wasm::runner::fault(&scenario, 20));
    assert!(!tepsim_wasm::runner::fault(&scenario, 2));
}

/// A mid-run fault request must not silently do nothing, and must not silently
/// change the run either. It records the request and says a rebuild is needed.
#[test]
fn a_mid_run_fault_request_asks_for_a_rebuild() {
    let mut run = runner(short());
    assert!(!run.pending_restart());
    let _ = run.step_chunk(2);

    run.request_fault(6, true).expect("IDV(6) is in range");
    assert!(run.pending_restart(), "the request needs a rebuild");
    assert!(
        !run.scenario().disturbances[5],
        "and it must not have altered the run already under way"
    );
    assert!(
        run.requested_scenario().disturbances[5],
        "but the requested scenario carries it"
    );

    assert_eq!(run.request_fault(21, true), Err(ConfigError::FaultIndex));
}

/// A browser can send any number for any field. Every one of these would
/// otherwise reach a division by zero or an `f64` to `usize` cast that
/// saturates, and in wasm a panic aborts the module.
#[test]
fn invalid_scenarios_are_rejected_with_the_reason() {
    let bad = |mutate: fn(&mut Scenario)| {
        let mut scenario = short();
        mutate(&mut scenario);
        Runner::new(scenario).err()
    };

    assert_eq!(bad(|s| s.seed = 0.0), Some(ConfigError::Seed));
    assert_eq!(bad(|s| s.seed = -1.0), Some(ConfigError::Seed));
    assert_eq!(bad(|s| s.seed = f64::NAN), Some(ConfigError::Seed));
    assert_eq!(bad(|s| s.seed = f64::INFINITY), Some(ConfigError::Seed));

    assert_eq!(bad(|s| s.hours = 0.0), Some(ConfigError::Hours));
    assert_eq!(bad(|s| s.hours = -1.0), Some(ConfigError::Hours));
    assert_eq!(bad(|s| s.hours = f64::NAN), Some(ConfigError::Hours));
    assert_eq!(bad(|s| s.hours = f64::INFINITY), Some(ConfigError::Hours));

    assert_eq!(bad(|s| s.step_hours = 0.0), Some(ConfigError::StepHours));
    assert_eq!(bad(|s| s.step_hours = -1.0), Some(ConfigError::StepHours));
    assert_eq!(
        bad(|s| s.step_hours = f64::NAN),
        Some(ConfigError::StepHours)
    );

    // The one that would panic rather than misbehave: `Scenario::samples`
    // divides by this.
    assert_eq!(bad(|s| s.sample_every = 0), Some(ConfigError::SampleEvery));

    assert_eq!(
        bad(|s| {
            // `Scenario::steps` rounds, and rounding is half away from zero, so
            // a step of exactly twice the duration would round back up to one.
            s.hours = 1.0;
            s.step_hours = 4.0;
        }),
        Some(ConfigError::NoSteps),
        "less than one step"
    );
    assert_eq!(
        bad(|s| s.sample_every = 100_000_000),
        Some(ConfigError::NoSamples),
        "a cadence longer than the run"
    );
    assert_eq!(
        bad(|s| s.hours = 1e9),
        Some(ConfigError::TooLong),
        "a run past the step limit"
    );

    for error in [
        ConfigError::Seed,
        ConfigError::Hours,
        ConfigError::StepHours,
        ConfigError::SampleEvery,
        ConfigError::NoSteps,
        ConfigError::NoSamples,
        ConfigError::TooLong,
        ConfigError::FaultIndex,
    ] {
        assert!(!error.message().is_empty());
    }
}

#[test]
fn the_limits_are_enforced_exactly() {
    let mut at_cap = Scenario::baseline();
    at_cap.step_hours = 1.0;
    at_cap.sample_every = 1;
    at_cap.hours = MAX_STEPS as f64;
    assert_eq!(
        validate(&at_cap).map(|p| p.steps),
        Err(ConfigError::TooLong),
        "MAX_STEPS steps would also be MAX_STEPS samples, over the sample cap"
    );

    let mut sample_capped = at_cap;
    sample_capped.hours = (MAX_SAMPLES + 1) as f64;
    sample_capped.sample_every = 1;
    assert_eq!(validate(&sample_capped), Err(ConfigError::TooLong));

    let mut just_under = at_cap;
    just_under.hours = MAX_SAMPLES as f64;
    assert_eq!(
        validate(&just_under).map(|p| p.samples),
        Ok(MAX_SAMPLES),
        "exactly at the sample cap is allowed"
    );
}

/// The scenario digest must respond to every field a run's numbers depend on,
/// or two different runs share an identifier and a shared link means nothing.
#[test]
fn the_scenario_digest_covers_every_field() {
    let base = Scenario::baseline();
    let moved: Vec<(&str, Scenario)> = vec![
        ("seed", base.with_seed(7.0)),
        ("hours", base.with_hours(1.0)),
        ("step_hours", {
            let mut s = base;
            s.step_hours = 0.5 / 3600.0;
            s
        }),
        ("sample_every", base.sampling_every(90)),
        ("disturbances", base.with_fault(3)),
        ("controlled", base.open_loop()),
        ("driver_forces_idv12", {
            let mut s = base;
            s.driver_forces_idv12 = false;
            s
        }),
        ("trip_ends_the_run", {
            let mut s = base;
            s.quirks.trip_ends_the_run = true;
            s
        }),
        ("integrator", base.with_integrator(Integrator::Rk4)),
    ];

    let baseline = scenario_digest(&base);
    for (field, scenario) in moved {
        assert_ne!(
            scenario_digest(&scenario),
            baseline,
            "changing {field} left the scenario digest unmoved"
        );
    }
    assert_eq!(scenario_digest(&base), baseline, "and it is stable");
}

/// The integrator changes every number, which is the point of offering it and
/// the reason it has to be in the digest.
#[test]
fn a_different_integrator_is_a_different_run() {
    let euler = short();
    let rk4 = euler.with_integrator(Integrator::Rk4);
    assert!(euler.integrator.is_faithful());
    assert!(!rk4.integrator.is_faithful());

    let mut a = runner(euler);
    let mut b = runner(rk4);
    assert_ne!(a.run_to_end(), b.run_to_end());
}

/// Open loop is the diagnostic mode: the plant trips on reactor pressure after
/// about three hours. That the bindings report a trip rather than pretending
/// the run completed is what a browser needs in order to say what happened.
#[test]
fn an_open_loop_run_trips_and_says_so() {
    let mut run = runner(Scenario::baseline().with_hours(6.0).open_loop());
    let _ = run.run_to_end();

    match run.outcome() {
        Some(Outcome::Tripped { hours, cause, .. }) => {
            assert!(
                (1.0..6.0).contains(&hours),
                "an open-loop plant trips within the run, not at {hours} h"
            );
            assert!(cause.is_some(), "and the cause is recorded");
        }
        other => panic!("expected an open-loop trip, got {other:?}"),
    }
    assert_eq!(run.outcome_name(), Some("tripped"));

    // The default is D-007 off: the plant freezes and keeps reporting, exactly
    // as teprob.f:807-811 does, so every planned sample still arrives.
    assert_eq!(
        run.emitted_samples(),
        run.plan().samples,
        "a trip must not truncate the run while D-007 is off"
    );
}

/// Row geometry is what a browser slices buffers with. If it moves, every chart
/// silently plots the wrong channel.
#[test]
fn row_geometry_matches_the_facade() {
    assert_eq!(channels::measurement_count(), 41);
    assert_eq!(channels::manipulated_count(), 12);
    assert_eq!(channels::channel_count(), 53);
    assert_eq!(ROW_WIDTH, 54, "a time column, then the 53 channels");
    assert_eq!(channels::column_ids().len(), ROW_WIDTH);
    assert_eq!(channels::column_labels().len(), ROW_WIDTH);
    assert_eq!(channels::column_units().len(), ROW_WIDTH);
}

/// The identifiers come from the facade and must stay aligned with it, or a
/// CSV header names one channel and holds another.
#[test]
fn the_column_ids_are_the_facades_channel_names() {
    let ids = channels::column_ids();
    assert_eq!(ids[0], "time_hours");
    for (i, name) in tepsim::channel_names().iter().enumerate() {
        assert_eq!(&ids[1 + i], name, "column {} drifted", i + 1);
    }
}

/// Duplicate names make a legend ambiguous. "Component A" is the description of
/// both `XMEAS(23)` and `XMEAS(29)`, which is why the labels carry their index.
#[test]
fn column_names_are_unique() {
    for (what, names) in [
        ("ids", channels::column_ids()),
        ("labels", channels::column_labels()),
    ] {
        let mut sorted = names.clone();
        sorted.sort();
        sorted.dedup();
        assert_eq!(sorted.len(), names.len(), "{what} collide");
    }
}

#[test]
fn the_sampled_columns_are_the_analyser_channels() {
    let sampled = channels::sampled_columns();
    // Column 0 is time, so XMEAS(23..=41) sits at offsets 23 through 41.
    assert_eq!(sampled, (23u32..=41).collect::<Vec<_>>());
    assert_eq!(sampled.len(), 19);
    let ids = channels::column_ids();
    for offset in sampled {
        assert!(
            ids[offset as usize].starts_with("XMEAS_"),
            "offset {offset} is not a measurement"
        );
    }
}

/// The digest must be a function of bit patterns, not of values, or `-0.0` and
/// `0.0` collide and a sign flip in a channel goes unnoticed.
#[test]
fn the_digest_distinguishes_signed_zero_and_is_order_sensitive() {
    let digest = |values: &[f64]| {
        let mut hash = Fnv1a64::new();
        hash.write_slice(values);
        hash.finish()
    };
    assert_ne!(digest(&[0.0]), digest(&[-0.0]));
    assert_ne!(digest(&[1.0, 2.0]), digest(&[2.0, 1.0]));
    assert_ne!(digest(&[]), digest(&[0.0]));
    assert_eq!(digest(&[]), Fnv1a64::new().finish());
}

/// FNV-1a is a published algorithm with published test vectors. Checking one
/// pins the constants, so a mistyped prime cannot make every digest in the
/// project self-consistently wrong.
#[test]
fn the_hash_matches_the_published_fnv1a_vector() {
    let mut hash = Fnv1a64::new();
    for byte in b"foobar" {
        hash.write_u8(*byte);
    }
    assert_eq!(
        hash.finish(),
        0x8594_4171_f739_67e8,
        "FNV-1a 64 of \"foobar\""
    );
    assert_eq!(
        Fnv1a64::new().finish(),
        0xcbf2_9ce4_8422_2325,
        "the offset basis"
    );
}
