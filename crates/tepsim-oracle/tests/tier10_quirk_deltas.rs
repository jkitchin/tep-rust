//! Tier 10: every Class C quirk fix, with its delta measured.
//!
//! B-0057. `PLAN.org` requires that a Class C fix be implemented behind a flag,
//! that its effect be *measured* rather than argued about, and that it stay off
//! until someone signs it off. This is the measurement.
//!
//! There are two Class C deltas:
//!
//! - **D-007**, whether a shutdown ends the run or freezes the plant.
//!   `tepsim_core::QuirkFixes::trip_ends_the_run`.
//! - **D-011**, whether the driver forces `IDV(12)` on at eight hours whatever
//!   the scenario asked for. `tepsim::Scenario::driver_forces_idv12`.
//!
//! # Why the Tier 5 battery cannot measure D-007
//!
//! The backlog assumed both deltas could be measured by running the battery
//! with the fix on and off. That works for D-011 and **does not work for
//! D-007**, and the reason is worth stating because it is not a limitation of
//! the battery.
//!
//! Every statistic in the battery compares two ensembles of runs of the same
//! shape. Turning D-007 on makes a tripped run *stop*, so it produces fewer
//! samples than its counterpart, and there is no correspondence between the two
//! sample sets to compare. A Kolmogorov-Smirnov statistic between a
//! 960-sample run and a 220-sample one is not zero and not meaningful: it is
//! measuring the truncation.
//!
//! So D-007's delta is measured in the terms the fix actually changes: which
//! scenarios trip, at which step, and how much of the run is lost.

#![cfg(feature = "oracle")]
// Exact comparisons: a run either is bit-identical to another or it is not,
// and that is the property under test throughout.
#![allow(
    clippy::float_cmp,
    reason = "bit equality between two runs is the property under test"
)]

use tepsim::{Scenario, Simulation};
use tepsim_oracle::tier5::Scenario as TierScenario;

/// Long enough to cross the eight-hour mark that D-011 turns on, and long
/// enough for the open-loop plant to trip for D-007.
const HOURS: f64 = 10.0;

/// One scenario's worth of D-007 evidence.
#[derive(Debug)]
struct TripDelta {
    label: String,
    /// The step it tripped at, if it did.
    tripped_at: Option<usize>,
    /// Samples with the fix off: the plant freezes and keeps reporting.
    frozen_samples: usize,
    /// Samples with the fix on: the run stops.
    stopped_samples: usize,
}

impl TripDelta {
    /// How much of the run the fix discards.
    fn lost_fraction(&self) -> f64 {
        if self.frozen_samples == 0 {
            return 0.0;
        }
        1.0 - (self.stopped_samples as f64 / self.frozen_samples as f64)
    }
}

/// **D-007's delta.** Measured as truncation, not as a distributional distance.
#[test]
fn d007_the_trip_fix_discards_the_frozen_tail() {
    let mut rows = Vec::new();

    // Open loop trips reliably; several closed-loop faults do not trip at all,
    // and both cases belong in the table.
    let cases: Vec<(String, Scenario)> = core::iter::once((
        "nominal, open loop".to_string(),
        Scenario::baseline().with_hours(HOURS).open_loop(),
    ))
    .chain((1..=6).map(|n| {
        (
            format!("IDV({n}), closed loop"),
            Scenario::fault(n).with_hours(HOURS),
        )
    }))
    .collect();

    for (label, scenario) in cases {
        let frozen = Simulation::new(scenario).run();
        let mut fixed = scenario;
        fixed.quirks.trip_ends_the_run = true;
        let stopped = Simulation::new(fixed).run();

        rows.push(TripDelta {
            label,
            tripped_at: frozen.tripped_at(),
            frozen_samples: frozen.samples.len(),
            stopped_samples: stopped.samples.len(),
        });
    }

    println!(
        "D-007 over {HOURS} h. `frozen` is the default (teprob.f:807-811 zeroes \
         the derivatives and the plant keeps reporting); `stopped` is the fix.\n"
    );
    println!(
        "  {:<22} {:>10} {:>9} {:>9} {:>8}",
        "scenario", "trip step", "frozen", "stopped", "lost"
    );
    let mut tripped = 0;
    for row in &rows {
        println!(
            "  {:<22} {:>10} {:>9} {:>9} {:>7.1}%",
            row.label,
            row.tripped_at
                .map_or_else(|| "-".to_string(), |s| s.to_string()),
            row.frozen_samples,
            row.stopped_samples,
            row.lost_fraction() * 100.0
        );
        if row.tripped_at.is_some() {
            tripped += 1;
        }
    }

    // The table has to contain both outcomes or it measures nothing.
    assert!(
        tripped > 0,
        "no scenario tripped, so D-007 changed nothing here"
    );
    assert!(
        tripped < rows.len(),
        "every scenario tripped, so the table cannot show that the fix is a \
         no-op where there is no trip"
    );

    for row in &rows {
        if row.tripped_at.is_some() {
            assert!(
                row.stopped_samples < row.frozen_samples,
                "{}: the fix kept as many samples as the default",
                row.label
            );
        } else {
            // Where nothing trips the fix is exactly a no-op, and that is the
            // other half of the claim.
            assert_eq!(
                row.stopped_samples, row.frozen_samples,
                "{}: the fix changed a run that never tripped",
                row.label
            );
        }
    }
}

/// Up to the trip, the two settings are bit-identical.
///
/// That is what makes D-007 a *truncation* rather than a different simulation:
/// the fix does not change any number, it decides how many of them there are.
#[test]
fn d007_changes_nothing_before_the_trip() {
    let scenario = Scenario::baseline().with_hours(HOURS).open_loop();
    let frozen = Simulation::new(scenario).run();
    let mut fixed = scenario;
    fixed.quirks.trip_ends_the_run = true;
    let stopped = Simulation::new(fixed).run();

    assert!(stopped.samples.len() > 2, "nothing to compare");
    for (index, (a, b)) in stopped.samples.iter().zip(&frozen.samples).enumerate() {
        assert_eq!(
            a.row(),
            b.row(),
            "sample {index} differs before the fix can have acted"
        );
    }
    println!(
        "the first {} samples are identical; the default then reports {} more \
         from the frozen plant",
        stopped.samples.len(),
        frozen.samples.len() - stopped.samples.len()
    );
}

/// **D-011's delta**, measured with the battery, which does work here because
/// both settings produce runs of the same shape.
///
/// The numbers are B-0040's, recomputed: the two part at the walk's next
/// segment boundary rather than at the eight-hour mark, and the worst
/// difference is around a tenth.
#[test]
fn d011_the_forced_disturbance_moves_the_plant_by_a_tenth() {
    use tepsim_oracle::Oracle;
    use tepsim_oracle::tier5::battery::compare;
    use tepsim_oracle::tier5::{run_port, seed, start};

    let mut oracle = Oracle::lock();
    let start = start(&mut oracle);
    drop(oracle);

    let seeds = 6;
    let hours = HOURS as usize;

    // `run_port` always forces IDV(12); build the unforced side directly.
    let forced: Vec<_> = (0..seeds)
        .map(|s| run_port(&start, TierScenario::NOMINAL, seed(s), hours))
        .collect();
    let honest: Vec<_> = (0..seeds)
        .map(|s| {
            let mut scenario = Scenario::baseline()
                .with_hours(hours as f64)
                .with_seed(seed(s));
            scenario.driver_forces_idv12 = false;
            let run = Simulation::new(scenario).run();
            tepsim_oracle::tier5::Run {
                scenario: TierScenario::NOMINAL,
                seed: seed(s),
                samples: run.samples.iter().map(tepsim::Sample::row).collect(),
                tripped: None,
            }
        })
        .collect();

    let report = compare(TierScenario::NOMINAL, &forced, &honest);

    let worst = report
        .variables
        .iter()
        .filter(|v| !v.constant)
        .map(|v| {
            let d = libm::fabs(v.paired.welch.difference);
            let m = v.paired.margin;
            (if m > 0.0 { d / m } else { 0.0 }, v.variable + 1)
        })
        .fold((0.0_f64, 0), |a, b| if b.0 > a.0 { b } else { a });

    println!(
        "D-011 over {hours} h at {seeds} seeds: worst paired mean difference is \
         {:.2} margins, on channel {} (a margin is a tenth of that channel's \
         own standard deviation)",
        worst.0, worst.1
    );

    // A *mean over the whole run* understates this considerably, and the
    // number is reported with that caveat rather than presented as the
    // headline. `IDV(12)` arrives at hour eight and its effect does not reach
    // the plant until channel 6's next walk boundary, so of these ten hours it
    // is active for under two, and the eight quiet hours dilute the mean by
    // roughly a factor of five.
    //
    // The undiluted figure is B-0040's, measured instantaneously rather than
    // as a mean: worst 1.092e-1 relative at XMEAS(37) over the run, 3.262e-2
    // at the end. That is what a sign-off should weigh.
    //
    // What is asserted here is only that the quirk is not cosmetic. A third of
    // the equivalence margin is a difference that, between the port and the
    // Fortran, this project would call a failure.
    assert!(
        worst.0 > 0.1,
        "forcing IDV(12) shifted no channel's mean by even a tenth of a \
         margin, so either it is cosmetic after all or this comparison is not \
         reaching it"
    );
}

/// Both Class C fixes are off by default, which is the standing rule.
#[test]
fn every_class_c_fix_is_off_by_default() {
    let scenario = Scenario::baseline();
    assert!(
        !scenario.quirks.trip_ends_the_run,
        "D-007's fix is on by default"
    );
    assert!(
        scenario.driver_forces_idv12,
        "D-011's quirk is not being reproduced by default"
    );
    // And the extension, which is not a fix but is also off.
    assert!(!scenario.extensions.continuous_disturbances);

    assert_eq!(
        Simulation::new(scenario).run().outcome,
        Simulation::new(Scenario::default()).run().outcome
    );
}
