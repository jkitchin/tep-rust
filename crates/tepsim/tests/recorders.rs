//! The recorder sinks.
//!
//! What these check is that a sink keeps *exactly* what it claims and that the
//! wrappers compose, because a sink that silently dropped or duplicated a
//! sample would corrupt every downstream number without ever failing.

#![allow(
    clippy::float_cmp,
    reason = "the run is deterministic, so recorded values are exactly equal"
)]

use tepsim::{
    CHANNELS, Columnar, Csv, Decimating, Outcome, Recorder, Ring, Scenario, Selecting, Simulation,
};

fn short() -> Scenario {
    Scenario::baseline().with_hours(1.0)
}

/// Recording into a sink and collecting into a `Run` see the same samples.
#[test]
fn a_recorder_sees_exactly_what_run_collects() {
    let scenario = short();
    let collected = Simulation::new(scenario).run();

    let mut columnar = Columnar::with_capacity(scenario.samples());
    let outcome = Simulation::new(scenario).run_into(&mut columnar);

    assert_eq!(outcome, collected.outcome);
    assert_eq!(columnar.len(), collected.samples.len());
    for channel in 0..CHANNELS {
        assert_eq!(
            columnar.column(channel),
            collected.column(channel).as_slice(),
            "channel {channel}"
        );
    }
    assert_eq!(
        columnar.steps(),
        collected
            .samples
            .iter()
            .map(|s| s.step)
            .collect::<Vec<_>>()
            .as_slice()
    );
}

/// The unit sink keeps nothing, which is what makes it the right way to time a
/// run or to check that it completes.
#[test]
fn the_unit_sink_records_nothing_and_still_runs() {
    let outcome = Simulation::new(short()).run_into(&mut ());
    assert_eq!(outcome, Outcome::Completed);
}

// ---------------------------------------------------------------------------
// Ring
// ---------------------------------------------------------------------------

#[test]
fn the_ring_keeps_the_last_samples_in_order() {
    let scenario = short();
    let all = Simulation::new(scenario).run();
    let capacity = 7;

    let mut ring = Ring::new(capacity);
    let _ = Simulation::new(scenario).run_into(&mut ring);

    assert_eq!(ring.len(), capacity);
    assert_eq!(ring.seen(), all.samples.len());

    // Oldest first, and they are the last `capacity` samples of the run.
    let held: Vec<_> = ring.iter().copied().collect();
    let expected = &all.samples[all.samples.len() - capacity..];
    assert_eq!(held.as_slice(), expected);

    // Step numbers ascend, which is the property a naive ring gets wrong by
    // returning the buffer in physical rather than logical order.
    for pair in held.windows(2) {
        assert!(pair[1].step > pair[0].step);
    }
}

#[test]
fn a_ring_larger_than_the_run_keeps_everything() {
    let scenario = short();
    let all = Simulation::new(scenario).run();
    let mut ring = Ring::new(all.samples.len() * 2);
    let _ = Simulation::new(scenario).run_into(&mut ring);
    assert_eq!(ring.len(), all.samples.len());
    assert_eq!(ring.iter().copied().collect::<Vec<_>>(), all.samples);
}

#[test]
#[should_panic(expected = "records nothing")]
fn a_ring_of_zero_capacity_is_refused() {
    let _ = Ring::new(0);
}

// ---------------------------------------------------------------------------
// Decimating
// ---------------------------------------------------------------------------

/// One sample in `factor`, *starting with the first*.
///
/// Keeping the `factor`-th instead would start the series late, which shifts
/// every time axis by a factor-dependent amount and is the kind of error that
/// shows up much later as a mysterious offset.
#[test]
fn decimating_keeps_one_in_n_starting_at_the_first() {
    let scenario = short();
    let all = Simulation::new(scenario).run();

    for factor in [1_usize, 2, 5, 10] {
        let mut sink = Decimating::new(Columnar::new(), factor);
        let _ = Simulation::new(scenario).run_into(&mut sink);
        let kept = sink.into_inner();

        let expected = all.samples.len().div_ceil(factor);
        assert_eq!(kept.len(), expected, "factor {factor}");

        let steps: Vec<usize> = all.samples.iter().step_by(factor).map(|s| s.step).collect();
        assert_eq!(kept.steps(), steps.as_slice(), "factor {factor}");
        assert_eq!(
            kept.steps().first(),
            all.samples.first().map(|s| &s.step),
            "factor {factor} did not start at the first sample"
        );
    }
}

#[test]
#[should_panic(expected = "keep nothing")]
fn decimating_by_zero_is_refused() {
    let _ = Decimating::new(Columnar::new(), 0);
}

// ---------------------------------------------------------------------------
// Selecting
// ---------------------------------------------------------------------------

/// Selected channels keep their values, the rest are zeroed, and every index
/// still means the same variable.
#[test]
fn selecting_blanks_rather_than_reshapes() {
    let scenario = short();
    let all = Simulation::new(scenario).run();

    // XMEAS(7) is channel 6, XMV(1) is channel 41.
    let wanted = [6_usize, 41];
    let mut sink = Selecting::new(Columnar::new(), &wanted);
    let _ = Simulation::new(scenario).run_into(&mut sink);
    let kept = sink.into_inner();

    assert_eq!(kept.len(), all.samples.len());
    for channel in 0..CHANNELS {
        if wanted.contains(&channel) {
            assert_eq!(
                kept.column(channel),
                all.column(channel).as_slice(),
                "channel {channel} was selected and should be intact"
            );
        } else {
            assert!(
                kept.column(channel).iter().all(|v| *v == 0.0),
                "channel {channel} was not selected and is not zero"
            );
        }
    }
    // And the selected ones are not accidentally all zero too.
    assert!(kept.column(6).iter().any(|v| *v > 2000.0));
}

#[test]
#[should_panic(expected = "records nothing")]
fn selecting_no_channels_is_refused() {
    let _ = Selecting::new(Columnar::new(), &[]);
}

// ---------------------------------------------------------------------------
// Composition
// ---------------------------------------------------------------------------

#[test]
fn the_wrappers_compose() {
    let scenario = short();
    let all = Simulation::new(scenario).run();

    let mut sink = Decimating::new(Selecting::new(Ring::new(4), &[6]), 3);
    let _ = Simulation::new(scenario).run_into(&mut sink);

    let ring = sink.into_inner().into_inner();
    assert_eq!(ring.len(), 4);
    assert_eq!(ring.seen(), all.samples.len().div_ceil(3));
    // Selected channel intact, everything else blanked.
    for sample in ring.iter() {
        assert!(sample.measurements[6] > 2000.0);
        assert_eq!(sample.measurements[0], 0.0);
    }
}

// ---------------------------------------------------------------------------
// CSV
// ---------------------------------------------------------------------------

#[test]
fn the_csv_sink_writes_a_header_and_a_row_per_sample() {
    let scenario = short().sampling_every(1200);
    let mut sink = Csv::new(String::new());
    let outcome = Simulation::new(scenario).run_into(&mut sink);
    assert_eq!(outcome, Outcome::Completed);
    assert_eq!(sink.error(), None);

    let text = sink.into_inner();
    let lines: Vec<&str> = text.lines().collect();
    assert_eq!(lines.len(), scenario.samples() + 1, "header plus samples");
    assert!(lines[0].starts_with("step,hours,XMEAS_1_A_feed"));
    // 53 channels plus step and hours.
    assert_eq!(lines[0].split(',').count(), CHANNELS + 2);
    for line in &lines[1..] {
        assert_eq!(line.split(',').count(), CHANNELS + 2);
    }
}

/// Seventeen significant digits round-trip an `f64` exactly.
///
/// This is what makes a recorded dataset reproducible rather than
/// approximately reproducible, and it is worth a test because a shorter format
/// looks identical to a reader and silently loses the last bits.
#[test]
fn the_csv_values_round_trip_exactly() {
    let scenario = short().sampling_every(1200);
    let run = Simulation::new(scenario).run();
    let mut sink = Csv::new(String::new());
    let _ = Simulation::new(scenario).run_into(&mut sink);
    let text = sink.into_inner();

    for (line, sample) in text.lines().skip(1).zip(&run.samples) {
        let fields: Vec<&str> = line.split(',').collect();
        for (index, value) in sample.row().into_iter().enumerate() {
            let parsed: f64 = fields[index + 2].parse().expect("a number");
            assert_eq!(
                parsed.to_bits(),
                value.to_bits(),
                "channel {index} did not round-trip: {} against {value}",
                fields[index + 2]
            );
        }
    }
}

#[test]
fn the_csv_label_columns_appear_only_when_asked() {
    let scenario = Scenario::fault(4).with_hours(0.2).sampling_every(180);

    let mut plain = Csv::new(String::new());
    let _ = Simulation::new(scenario).run_into(&mut plain);
    let plain = plain.into_inner();
    assert!(!plain.lines().next().expect("a header").contains("fault"));

    let mut labelled = Csv::new(String::new()).with_labels();
    let _ = Simulation::new(scenario).run_into(&mut labelled);
    let labelled = labelled.into_inner();
    let header = labelled.lines().next().expect("a header");
    assert!(header.ends_with(",fault,hours_since_onset"));
    // Every row names the fault that is on.
    for line in labelled.lines().skip(1) {
        let fields: Vec<&str> = line.split(',').collect();
        assert_eq!(fields[CHANNELS + 2], "4", "in {line}");
    }
}

/// A write error stops the sink rather than being swallowed.
#[test]
fn a_csv_write_error_is_reported() {
    struct Failing(usize);
    impl core::fmt::Write for Failing {
        fn write_str(&mut self, _s: &str) -> core::fmt::Result {
            if self.0 == 0 {
                return Err(core::fmt::Error);
            }
            self.0 -= 1;
            Ok(())
        }
    }

    let mut sink = Csv::new(Failing(3));
    let _ = Simulation::new(short().sampling_every(1200)).run_into(&mut sink);
    assert!(
        sink.error().is_some(),
        "the sink swallowed a write failure and reported success"
    );
}
