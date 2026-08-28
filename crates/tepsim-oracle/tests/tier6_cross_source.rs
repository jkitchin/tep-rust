//! Tier 6: the cross-source detector experiment.
//!
//! B-0050. Train a PCA fault detector on one simulator's fault-free data,
//! evaluate it on both simulators' faulted data, and ask whether the difference
//! between the two evaluations is larger than the difference the reference
//! simulator shows against itself. Then swap the training source and do it
//! again.
//!
//! Every detector parameter is fixed in [`tepsim_oracle::tier6`], with its
//! reason, and none of them was chosen after a number was on the screen.
//!
//! `cargo test` runs [`Battery::SMOKE`], which **cannot conclude anything**:
//! six test seeds admit ten distinct half-splits, the smallest attainable
//! permutation p-value is `1/11 = 0.09`, and every verdict comes back
//! undecided. What it does prove is that the whole pipeline runs, that both
//! sources produce runs of the same shape, that the detector separates faulted
//! from fault-free data at all, and that no metric is `NaN`. `TEP_TIER6=full`
//! runs the same code at a size that gates.

#![cfg(feature = "oracle")]

use tepsim_core::math;
use tepsim_oracle::Oracle;
use tepsim_oracle::tier5::{Scenario, VARIABLES, run_port, seed, start};
use tepsim_oracle::tier6::{
    ALPHA, Battery, CALIBRATION_SPLITS, CONFIDENCE, Experiment, LAG_COUNTS, PERSISTENCE, Source,
    calibrate, evaluate, fit, half_splits, score, separation,
};
use tepsim_stats::Summary;

/// Two groups whose means differ by exactly one, at two group sizes.
///
/// The property [`separation`] exists for: doubling both groups leaves the
/// difference of means alone and multiplies the statistic by `sqrt(2)`, so a
/// half-split of `n/2` against `n/2` and a cross-source comparison of `n`
/// against `n` land on the same scale. Without this the null would be wider
/// than the statistic it judges and every verdict would drift toward passing.
#[test]
fn separation_puts_unequal_group_sizes_on_one_scale() {
    let small = separation(&[1.0, 1.0], &[0.0, 0.0]);
    let large = separation(&[1.0; 4], &[0.0; 4]);
    // 1 / sqrt(1/2 + 1/2) and 1 / sqrt(1/4 + 1/4).
    assert!(
        (small - 1.0).abs() < 1e-15,
        "two against two should give exactly one, got {small}"
    );
    assert!(
        (large - core::f64::consts::SQRT_2).abs() < 1e-15,
        "four against four should give sqrt(2), got {large}"
    );
    // The raw difference of means is one in both cases, which is precisely why
    // a raw difference cannot be compared across group sizes.
    assert!(separation(&[], &[1.0]).is_nan());
}

/// The split generator returns **distinct** splits, and no more than exist.
///
/// Six items have `C(6,3)/2 = 10` balanced half-splits. Asking for twenty must
/// return ten, not twenty with repeats: a permutation p-value whose denominator
/// counts the same split eleven times is not a p-value. Twenty items have
/// 92,378, so all twenty requested are distinct.
#[test]
fn half_splits_are_distinct_and_balanced() {
    let six = half_splits(6, CALIBRATION_SPLITS);
    assert_eq!(
        six.len(),
        10,
        "six items have exactly ten balanced half-splits"
    );
    let twenty = half_splits(20, CALIBRATION_SPLITS);
    assert_eq!(twenty.len(), CALIBRATION_SPLITS);
    for (left, right) in &twenty {
        assert_eq!(left.len(), 10);
        assert_eq!(right.len(), 10);
        let mut all: Vec<usize> = left.iter().chain(right).copied().collect();
        all.sort_unstable();
        assert_eq!(
            all,
            (0..20).collect::<Vec<_>>(),
            "the halves must partition"
        );
    }
    // And the same call twice gives the same splits, so a verdict that moved is
    // the port moving and not the calibration.
    assert_eq!(twenty, half_splits(20, CALIBRATION_SPLITS));
}

/// The smoke battery cannot gate, and says so.
///
/// This is asserted rather than commented because it is the one thing a reader
/// of a green `cargo test` most needs to know.
#[test]
fn the_smoke_battery_declines_to_decide() {
    let splits = half_splits(Battery::SMOKE.test_seeds, CALIBRATION_SPLITS);
    assert!(
        splits.len() + 1 < (1.0 / ALPHA) as usize,
        "the smoke battery must not be able to reject at {ALPHA}, but it has \
         {} splits",
        splits.len()
    );
    assert!(
        half_splits(Battery::FULL.test_seeds, CALIBRATION_SPLITS).len() + 1
            >= (1.0 / ALPHA) as usize,
        "the full battery must be able to reject"
    );
}

/// `XMEAS(7)`, reactor pressure, zero-based.
const PRESSURE: usize = 6;

/// How the deviations of one variable are distorted, for the control below.
#[derive(Clone, Copy, Debug)]
enum Distortion {
    /// Add this many of the variable's own standard deviations to every sample.
    Shift(f64),
    /// Multiply the deviations about the run mean by this factor.
    Scale(f64),
}

impl Distortion {
    fn label(self) -> String {
        match self {
            Self::Shift(k) => format!("mean +{k} sd"),
            Self::Scale(k) => format!("sd x {k}"),
        }
    }

    /// Apply it to one variable of a copy of the runs.
    fn apply(
        self,
        runs: &[tepsim_oracle::tier5::Run],
        variable: usize,
    ) -> Vec<tepsim_oracle::tier5::Run> {
        let mut out = runs.to_vec();
        for run in &mut out {
            let column: Vec<f64> = run.samples.iter().map(|row| row[variable]).collect();
            let summary = Summary::of(&column);
            let (mean, sd) = (summary.mean(), summary.sd());
            assert!(sd > 0.0, "variable {variable} never moved");
            for row in &mut run.samples {
                row[variable] = match self {
                    Self::Shift(k) => row[variable] + k * sd,
                    Self::Scale(k) => mean + k * (row[variable] - mean),
                };
            }
        }
        out
    }
}

/// The pipeline catches a source that really is different, and how different it
/// has to be.
///
/// The cross-source statistics in the main test come out at exactly zero, which
/// is a pass and, on its own, an uninformative one. This is the control: the
/// same fit-score-threshold-calibrate path is fed one source and a distorted
/// copy of it, and the cross-source value has to land strictly above every
/// within-source draw.
///
/// The gate is a mean shift in one variable, swept until it is caught. That is
/// a difference the detector's own metrics can express: a shifted operating
/// point breaks the correlation structure the model learned and lands in the
/// residual subspace, which is what SPE is for. Asserting that *some* shift up
/// to eight standard deviations moves the false alarm rate outside the
/// reference's own spread proves the whole path works, and the smallest shift
/// that does it is a sensitivity number worth recording.
///
/// The variance collapse `PLAN.org` names as the existing Python port's real
/// failure, reactor pressure at a standard deviation of 8.10 against the
/// Fortran's 61.48, is run too, and it is **not** asserted on. That is a
/// finding rather than a concession: shrinking a variable's spread makes the
/// data less extreme, so it can only ever *lower* an alarm rate, and there is
/// not much room below a false alarm rate of a few percent. Tier 5's TOST on
/// log variances is the gate for that failure mode and it catches it with two
/// orders to spare. A tier that measures downstream task performance is not the
/// instrument for a second-moment defect, and pretending otherwise here would
/// put a gate where it cannot bite.
///
/// The port alone, so this costs twelve four-hour runs and no Fortran.
#[test]
fn a_source_that_really_differs_is_caught() {
    let hours = 4;
    let start = {
        let mut oracle = Oracle::lock();
        start(&mut oracle)
    };
    let training: Vec<_> = (0..6)
        .map(|s| run_port(&start, Scenario::NOMINAL, seed(s), hours))
        .collect();
    let clean: Vec<_> = (6..12)
        .map(|s| run_port(&start, Scenario::NOMINAL, seed(s), hours))
        .collect();

    let model = fit(&training, 0);
    let limits = model.limits(CONFIDENCE);
    let reference = evaluate(&score(&model, &clean, 0), &limits, false);
    let splits = half_splits(clean.len(), CALIBRATION_SPLITS);

    println!(
        "\n-- control: distorting XMEAS({}) of one source --\n  {:<14} {:>10} {:>12} {:>12} {:>8}",
        PRESSURE + 1,
        "distortion",
        "statistic",
        "cross",
        "within max",
        "caught"
    );
    let caught = |distortion: Distortion| -> bool {
        let candidate = evaluate(
            &score(&model, &distortion.apply(&clean, PRESSURE), 0),
            &limits,
            false,
        );
        let mut any = false;
        for (name, a, b) in [
            ("T2", &reference.t_squared, &candidate.t_squared),
            ("SPE", &reference.spe, &candidate.spe),
        ] {
            let calibrated = calibrate("false alarm rate", &a.rate, &b.rate, &splits);
            let hit = calibrated.cross > calibrated.within_max();
            any |= hit;
            println!(
                "  {:<14} {name:>10} {:>12.4e} {:>12.4e} {:>8}",
                distortion.label(),
                calibrated.cross,
                calibrated.within_max(),
                if hit { "yes" } else { "no" }
            );
        }
        any
    };

    let shifts = [0.5, 1.0, 2.0, 4.0, 8.0];
    let smallest = shifts.into_iter().find(|k| caught(Distortion::Shift(*k)));

    // Run the variance collapse too, and report it. Not a gate; see the doc
    // comment for why Tier 5 owns that failure mode.
    let collapse = 8.10 / 61.48;
    let collapse_caught = caught(Distortion::Scale(collapse));
    println!(
        "  the variance collapse PLAN.org names (reactor pressure sd 8.10 \
         against 61.48, a factor of {collapse:.4}) is {} by a false alarm rate. \
         Tier 5's TOST on log variances is the gate for it. A smaller spread \
         can only lower an alarm rate, and there is little room below a few \
         percent.",
        if collapse_caught {
            "caught"
        } else {
            "NOT caught"
        }
    );

    assert!(
        smallest.is_some(),
        "no mean shift up to {} standard deviations of XMEAS({}) moved the \
         false alarm rate outside the reference's own run-to-run spread. The \
         pipeline cannot see a difference it must be able to see, so a passing \
         cross-source verdict from it would mean nothing.",
        shifts.last().copied().unwrap_or_default(),
        PRESSURE + 1
    );
    println!(
        "  smallest single-variable mean shift this detector notices: {} sd",
        smallest.unwrap_or_default()
    );
    assert_eq!(VARIABLES, 53, "the recorded variable count moved");
}

fn print_parameters(battery: Battery) {
    println!(
        "\nTier 6 cross-source detector experiment\n\
         =======================================\n\
         battery      : {} faults, {} training seeds, {} test seeds, {} h per run\n\
         runs         : {} per source, {} in total\n\
         detectors    : PCA at lags {:?}, trained on each source, {} confidence\n\
         retention    : smallest k explaining 90% of the variance\n\
         persistence  : {} consecutive alarms for a detection delay\n\
         libm         : {}",
        battery.faults,
        battery.train_seeds,
        battery.test_seeds,
        battery.hours,
        battery.runs_per_source(),
        2 * battery.runs_per_source(),
        LAG_COUNTS,
        CONFIDENCE,
        PERSISTENCE,
        if math::USES_SYSTEM_LIBM {
            "platform"
        } else {
            "vendored"
        }
    );
}

fn print_models(experiment: &Experiment) {
    println!(
        "\n-- fitted detectors --\n  {:<20} {:>6} {:>7} {:>9} {:>10} {:>12} {:>12}",
        "detector", "cols", "rows", "kept", "explained", "T2 limit", "SPE limit"
    );
    for detector in &experiment.detectors {
        println!(
            "  {:<20} {:>6} {:>7} {:>9} {:>10.4} {:>12.4} {:>12.4}",
            detector.label(),
            detector.columns,
            detector.training_rows,
            detector.components,
            detector.explained_variance,
            detector.limits.t_squared,
            detector.limits.spe
        );
    }
}

fn print_comparisons(experiment: &Experiment) {
    println!(
        "\n-- every comparison: the detector is fixed, the evaluation source is not --\n  \
         {:<20} {:<10} {:<22} {:>9} {:>9} {:>9} {:>9} {:>9} {:>7} {:>9}",
        "detector",
        "scenario",
        "metric",
        "fortran",
        "rust",
        "diff",
        "cross",
        "within",
        "p",
        "verdict"
    );
    for detector in &experiment.detectors {
        for (scenario, c) in detector.comparisons() {
            println!(
                "  {:<20} {:<10} {:<22} {:>9.4} {:>9.4} {:>9.2e} {:>9.2e} {:>9.2e} {:>7.3} {:>9}",
                detector.label(),
                scenario.label(),
                c.metric,
                c.fortran.mean(),
                c.port.mean(),
                c.difference,
                c.calibrated.cross,
                c.calibrated.within_max(),
                c.calibrated.p_value(),
                match c.passes() {
                    Some(true) => "pass",
                    Some(false) => "FAIL",
                    None => "undecided",
                }
            );
        }
    }
}

/// Why the cross-source statistics come out where they do.
///
/// Printed separately from the verdicts because it is the diagnostic, not the
/// gate: how far apart the two sources' monitoring statistics got, against how
/// close either got to its control limit.
fn print_agreement(experiment: &Experiment) {
    println!(
        "\n-- could any alarm have differed? --\n  {:<20} {:<10} {:>12} {:>12} {:>12} {:>12} {:>10}",
        "detector", "scenario", "T2 gap", "T2 margin", "SPE gap", "SPE margin", "gap/margin"
    );
    for detector in &experiment.detectors {
        for (scenario, a) in detector.agreements() {
            println!(
                "  {:<20} {:<10} {:>12.3e} {:>12.3e} {:>12.3e} {:>12.3e} {:>10.3e}",
                detector.label(),
                scenario.label(),
                a.t_squared_gap,
                a.t_squared_approach,
                a.spe_gap,
                a.spe_approach,
                a.decisive()
            );
        }
    }

    println!(
        "\n-- raw sample agreement, detector-independent --\n  {:<10} {:>14} {:>16} {:>10} {:>10}",
        "scenario", "max absolute", "max / variable sd", "variable", "sample"
    );
    for (scenario, t) in &experiment.trajectories {
        println!(
            "  {:<10} {:>14.4e} {:>16.4e} {:>10} {:>10}",
            scenario.label(),
            t.max_absolute,
            t.max_scaled,
            t.variable + 1,
            t.sample
        );
    }
}

/// The four train/test combinations, side by side, for the detection rate.
///
/// This is the table `PLAN.org` describes, and it is printed separately from
/// the statistics because the statistics answer "is the gap larger than the
/// noise" while this answers "what were the numbers".
fn print_four_combinations(experiment: &Experiment) {
    for lags in LAG_COUNTS {
        let Some(f) = experiment
            .detectors
            .iter()
            .find(|d| d.lags == lags && d.trained_on == Source::Fortran)
        else {
            continue;
        };
        let Some(r) = experiment
            .detectors
            .iter()
            .find(|d| d.lags == lags && d.trained_on == Source::Port)
        else {
            continue;
        };
        println!(
            "\n-- four combinations, {} --\n  {:<10} {:<22} {:>12} {:>12} {:>12} {:>12}",
            if lags == 0 {
                "PCA".to_string()
            } else {
                format!("DPCA({lags})")
            },
            "scenario",
            "metric",
            "F train/F eval",
            "F train/R eval",
            "R train/R eval",
            "R train/F eval"
        );
        for (metric, index) in [("T2 false alarm rate", 0), ("SPE false alarm rate", 1)] {
            let (a, b) = (&f.false_alarm[index], &r.false_alarm[index]);
            println!(
                "  {:<10} {:<22} {:>12.4} {:>12.4} {:>12.4} {:>12.4}",
                "nominal",
                metric,
                a.fortran.mean(),
                a.port.mean(),
                b.port.mean(),
                b.fortran.mean()
            );
        }
        for (fault_f, fault_r) in f.faults.iter().zip(&r.faults) {
            for index in 0..fault_f.metrics.len() {
                let (a, b) = (&fault_f.metrics[index], &fault_r.metrics[index]);
                println!(
                    "  {:<10} {:<22} {:>12.4} {:>12.4} {:>12.4} {:>12.4}",
                    fault_f.scenario.label(),
                    a.metric,
                    a.fortran.mean(),
                    a.port.mean(),
                    b.port.mean(),
                    b.fortran.mean()
                );
            }
        }
    }
}

/// The largest detection rate anywhere, on each source.
///
/// Used to prove the detector detects. A cross-source comparison of two
/// detectors that never fire would pass every gate and mean nothing.
fn best_detection_rates(experiment: &Experiment) -> (f64, f64) {
    let mut fortran = 0.0_f64;
    let mut port = 0.0_f64;
    for detector in &experiment.detectors {
        for fault in &detector.faults {
            for metric in &fault.metrics {
                if metric.metric.ends_with("detection rate") {
                    fortran = fortran.max(metric.fortran.mean());
                    port = port.max(metric.port.mean());
                }
            }
        }
    }
    (fortran, port)
}

/// The whole experiment, at whatever size `TEP_TIER6` selects.
#[test]
fn a_detector_cannot_tell_which_simulator_trained_it() {
    let battery = Battery::selected();
    print_parameters(battery);

    let mut oracle = Oracle::lock();
    let experiment = tepsim_oracle::tier6::run(&mut oracle, battery);

    print_models(&experiment);
    print_four_combinations(&experiment);
    print_comparisons(&experiment);
    print_agreement(&experiment);

    // Structural facts, true at any size.
    for detector in &experiment.detectors {
        assert!(
            detector.components >= 1,
            "{} retained no components",
            detector.label()
        );
        assert!(
            detector.limits.t_squared.is_finite() && detector.limits.spe.is_finite(),
            "{} produced a non-finite control limit",
            detector.label()
        );
        for (scenario, comparison) in detector.comparisons() {
            for (side, summary) in [("fortran", comparison.fortran), ("rust", comparison.port)] {
                assert!(
                    summary.mean().is_finite(),
                    "{} {} {} on {side} is not finite",
                    detector.label(),
                    scenario.label(),
                    comparison.metric
                );
            }
            assert!(
                comparison.calibrated.cross.is_finite(),
                "{} {} {} produced a non-finite cross-source statistic",
                detector.label(),
                scenario.label(),
                comparison.metric
            );
        }
        // A plant that trips in one simulator and not in the other at the same
        // seed is a difference no statistic should have to find. This holds at
        // any battery size, so it is asserted at every size.
        for fault in &detector.faults {
            assert_eq!(
                fault.fortran_tripped,
                fault.port_tripped,
                "{}: {} of the Fortran's runs tripped and {} of the port's",
                fault.scenario.label(),
                fault.fortran_tripped,
                fault.port_tripped
            );
        }
    }

    // The detector has to actually detect, or the comparison is vacuous.
    let (fortran_best, port_best) = best_detection_rates(&experiment);
    println!("\nbest detection rate anywhere: fortran {fortran_best:.4}, rust {port_best:.4}");
    assert!(
        fortran_best >= 0.5 && port_best >= 0.5,
        "no detector reached a 50% detection rate on any fault \
         (fortran {fortran_best:.4}, rust {port_best:.4}); a cross-source \
         comparison of detectors that never fire proves nothing"
    );

    let (zero, total) = experiment.identical();
    let (gap, margin, decisive) = experiment.worst_agreement();
    println!(
        "\n{zero} of {total} comparisons had a cross-source statistic of \
         exactly zero."
    );
    println!(
        "largest gap between the two sources' monitoring statistics: {gap:.4e}\n\
         closest either came to its control limit:                   {margin:.4e}\n\
         ratio:                                                      {decisive:.4e} \
         ({})",
        if decisive < 1.0 {
            "no alarm on any sample could have been decided differently"
        } else {
            "at least one alarm could have flipped"
        }
    );
    if let Some((scenario, t)) = experiment.worst_trajectory() {
        println!(
            "worst raw sample disagreement: {:.4e} absolute, {:.4e} of the \
             variable's own sd, on variable {} of {}",
            t.max_absolute,
            t.max_scaled,
            t.variable + 1,
            scenario.label()
        );
    }
    println!(
        "pairing: the largest seed-to-seed correlation between the two \
         sources' per-run values is {:.4}. Near zero means the paired runs have \
         separated by the time the metric is taken, so the unpaired null is the \
         right one; near one means the cross-source comparison is a paired one \
         and the unpaired null is conservative.",
        experiment.worst_pairing()
    );

    match experiment.worst() {
        Some((detector, scenario, metric, p)) => println!(
            "worst permutation p-value: {p:.4} at {detector} / {} / {metric}",
            scenario.label()
        ),
        None => println!("worst permutation p-value: none finite"),
    }

    let failures = experiment.failures();
    if experiment.gated() {
        println!(
            "\nGATING on {} distinct half-splits: a statistic fails when it is \
             the strict maximum of the {} draws.",
            experiment.splits,
            experiment.splits + 1
        );
        for (detector, scenario, metric) in &failures {
            println!("  FAIL {detector} / {} / {metric}", scenario.label());
        }
        assert!(
            failures.is_empty(),
            "{} cross-source statistic(s) sit outside the reference's own \
             run-to-run spread: {failures:?}",
            failures.len()
        );
    } else {
        println!(
            "\nREPORTING ONLY: {} distinct half-splits, so the smallest \
             attainable p-value is {:.3} and nothing can be rejected at \
             {ALPHA}. Set {}=full for the gating run.",
            experiment.splits,
            1.0 / (experiment.splits + 1) as f64,
            Battery::ENV
        );
        assert!(
            failures.is_empty(),
            "an ungated experiment reported a failure, which should be \
             unreachable: {failures:?}"
        );
    }
}
