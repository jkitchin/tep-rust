//! The Tier 5 battery, run at whatever size `TEP_TIER5` selects.
//!
//! `cargo test` runs the smoke battery, which *reports*: four seeds give three
//! distinct half-splits and a permutation test with three draws cannot reject
//! anything at 0.05, so a passing smoke run is evidence that the machinery
//! works and nothing more. `cargo xtask validate --tiers 5` sets
//! `TEP_TIER5=full` and the same code *gates*.
//!
//! That distinction is in the output of every run, not in a comment.

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

use tepsim_core::math;
use tepsim_oracle::Oracle;
use tepsim_oracle::tier5::battery::{
    ALPHA, CALIBRATION_SPLITS, MEAN_MARGIN_FRACTION, Report, VARIANCE_MARGIN_LOG, compare,
};
use tepsim_oracle::tier5::{Battery, Scenario, VARIABLES, run_fortran, run_port, seed, start};

/// Run one scenario on both sources and compare.
fn scenario_report(oracle: &mut Oracle, battery: Battery, scenario: Scenario) -> Report {
    let start = start(oracle);
    let reference: Vec<_> = (0..battery.seeds)
        .map(|s| run_fortran(oracle, &start, scenario, seed(s), battery.hours))
        .collect();
    let candidate: Vec<_> = (0..battery.seeds)
        .map(|s| run_port(&start, scenario, seed(s), battery.hours))
        .collect();
    compare(scenario, &reference, &candidate)
}

fn print_report(report: &Report) {
    println!(
        "\n=== {} ({} seeds, {}) ===",
        report.scenario.label(),
        report.seeds,
        if report.gated {
            "gating"
        } else {
            "reporting only: too few seeds to calibrate"
        }
    );
    println!(
        "  {:<4} {:>11} {:>11} {:>10} {:>10} {:>10} {:>10}",
        "var", "paired p", "var p", "KS", "energy", "ACF", "spectrum"
    );
    // The five worst variables by TOST p-value, so the table stays readable
    // while still naming whatever is closest to failing.
    let mut order: Vec<usize> = (0..report.variables.len()).collect();
    order.sort_by(|a, b| {
        report.variables[*b]
            .paired
            .p
            .total_cmp(&report.variables[*a].paired.p)
    });
    for index in order.iter().take(5) {
        let v = &report.variables[*index];
        println!(
            "  {:<4} {:>11.4} {:>11.4} {:>10.4} {:>10.4} {:>10.4} {:>10.4}",
            v.variable + 1,
            v.paired.p,
            v.variance.p,
            v.ks.p_value(),
            v.energy.p_value(),
            v.autocorrelation.p_value(),
            v.spectrum.p_value()
        );
    }
    println!(
        "  mean-test power: unpaired {:.2}, paired {:.4} x the margin",
        report.worst_mean_power(),
        report.worst_paired_power()
    );
    println!(
        "  gate uses the paired test ({})",
        if report.worst_paired_power() < 1.0 {
            "enough seeds to declare equivalence".to_string()
        } else {
            format!(
                "the margin needs about {} seeds; {} were run",
                report.seeds_for_power(),
                report.seeds
            )
        }
    );
    let undecided = report.undecided();
    if !undecided.is_empty() {
        println!(
            "  {} variable(s) undecided on the mean. Worst five, as (variable, \
             power, |difference| / margin):",
            undecided.len()
        );
        for (variable, power, gap) in undecided.iter().take(5) {
            println!("    XMV/XMEAS {variable:<3} power {power:.2}  gap {gap:.3}");
        }
    }
    println!(
        "  frobenius: cross {:.6e}, within-source max {:.6e}, p {:.4}",
        report.structure.frobenius.cross,
        report.structure.frobenius.within_max(),
        report.structure.frobenius.p_value()
    );
    if let Some((i, j, x, y)) = report.structure.worst_pair {
        println!("  worst correlation pair: ({i}, {j}) {x:.6} against {y:.6}");
    }
}

/// The battery over whatever scenarios the selected size covers.
#[test]
fn the_battery_finds_the_two_sources_equivalent() {
    let battery = Battery::selected();
    println!(
        "Tier 5 battery: {} scenarios x {} seeds x {} h, {} samples per run, \
         {} libm",
        battery.scenarios,
        battery.seeds,
        battery.hours,
        battery.samples(),
        if math::USES_SYSTEM_LIBM {
            "platform"
        } else {
            "vendored"
        }
    );
    println!(
        "margins: mean within {MEAN_MARGIN_FRACTION} of the reference sd; \
         variance within {:.1}%; everything else calibrated against {} \
         half-splits at alpha {ALPHA}",
        VARIANCE_MARGIN_LOG.exp_m1() * 100.0,
        CALIBRATION_SPLITS
    );

    let mut oracle = Oracle::lock();
    let mut all_failures = Vec::new();
    let mut all_undecided = Vec::new();
    let mut worst_gap = (0.0_f64, String::new(), 0usize);
    let mut gated = false;

    for scenario in Scenario::all().take(battery.scenarios) {
        let report = scenario_report(&mut oracle, battery, scenario);
        print_report(&report);
        gated |= report.gated;
        for (variable, what) in report.failures() {
            all_failures.push((scenario.label(), variable, what));
        }
        for (variable, power, gap) in report.undecided() {
            if gap > worst_gap.0 {
                worst_gap = (gap, scenario.label(), variable);
            }
            all_undecided.push((scenario.label(), variable, power, gap));
        }
    }

    if !all_undecided.is_empty() {
        println!(
            "\n{} scenario-variable mean tests could not decide: the interval \
             is wider than the margin, so TOST could not have declared \
             equivalence at any difference. Worst measured gap among them is \
             {:.3} of the margin, on {} variable {}.",
            all_undecided.len(),
            worst_gap.0,
            worst_gap.1,
            worst_gap.2
        );
        // An undecided test whose measured difference is *inside* the margin
        // is uninformative. One whose difference is well outside it is a
        // finding even though the test could not formally reject.
        assert!(
            worst_gap.0 < 1.0,
            "an undecided mean test has a measured difference of {:.3} margins \
             on {} variable {}. That is outside the margin even though the \
             test lacked the power to say so formally, and per CLAUDE.md it is \
             a BLOCKED item with the numbers rather than a widened margin.",
            worst_gap.0,
            worst_gap.1,
            worst_gap.2
        );
    }

    if !gated {
        println!(
            "\nreporting only: {} seeds give too few half-splits for the \
             permutation tests to reject at {ALPHA}. The paired mean test \
             does decide at this size, and its verdicts above are real. Run \
             `cargo xtask validate --tiers 5` for the whole gate.",
            battery.seeds
        );
        println!(
            "{} statistic(s) are outside their margins at this size, which at \
             this size is a statement about the sample and not about the port.",
            all_failures.len()
        );
        return;
    }

    assert!(
        all_failures.is_empty(),
        "Tier 5 found {} statistic(s) outside their margin: {all_failures:?}. \
         Per CLAUDE.md this is a BLOCKED item with the numbers, never a \
         widened margin.",
        all_failures.len()
    );
}

/// The battery machinery, checked against a case whose answer is known.
///
/// Comparing the Fortran with *itself* must find equivalence on every
/// statistic: the two sides are literally the same runs. If that fails, the
/// battery is broken rather than the port.
#[test]
fn the_battery_finds_a_source_equivalent_to_itself() {
    let mut oracle = Oracle::lock();
    let battery = Battery::SMOKE;
    let start = start(&mut oracle);
    let runs: Vec<_> = (0..battery.seeds)
        .map(|s| {
            run_fortran(
                &mut oracle,
                &start,
                Scenario::NOMINAL,
                seed(s),
                battery.hours,
            )
        })
        .collect();

    let report = compare(Scenario::NOMINAL, &runs, &runs);

    // Every distance statistic is exactly zero against itself.
    for variable in &report.variables {
        let v = variable.variable + 1;
        // KS and the spectrum are exactly zero: one is a count difference and
        // the other a ratio of identical sums. The energy distance is not,
        // because it is `2A - B - C` with the three terms of order one and the
        // answer of order zero, so its floor is an ulp of the terms.
        assert_eq!(variable.ks.cross, 0.0, "KS on variable {v}");
        assert_eq!(variable.spectrum.cross, 0.0, "spectrum on variable {v}");
        // The bound scales with the data, because the energy distance has the
        // units of the data. XMEAS(2) sits near 3664 and varies by a few, so
        // an absolute bound would be a statement about that variable's units
        // rather than about the arithmetic.
        let spread = tepsim_stats::Summary::of(
            &runs
                .iter()
                .flat_map(|r| r.series(variable.variable))
                .collect::<Vec<f64>>(),
        )
        .sd()
        .max(1.0);
        assert!(
            variable.energy.cross.abs() < 1e-12 * spread,
            "energy on variable {v} is {}, past the cancellation floor for a \
             variable of spread {spread:.3}",
            variable.energy.cross
        );
        if !variable.constant {
            // A constant series has no autocorrelation and no spectrum: the
            // ratio defining them is 0/0. XMV(12), the agitator, is exactly
            // that, and it is variable 53.
            assert_eq!(variable.autocorrelation.cross, 0.0, "ACF on variable {v}");
        }
        // And every permutation test passes, because zero cannot be the
        // strict maximum of a set that contains it. A constant variable's
        // serial statistics are undefined and have no null to compare against,
        // which `Calibrated::passes` reports as `None` rather than as a pass.
        for statistic in variable.calibrated() {
            if statistic.within.is_empty() {
                assert!(
                    variable.constant,
                    "{} on variable {v} has no null",
                    statistic.name
                );
                continue;
            }
            assert_eq!(
                statistic.p_value(),
                1.0,
                "{} on variable {v}",
                statistic.name
            );
        }
        // The moment tests are *centred* on zero, whatever their verdict.
        if !variable.constant {
            assert_eq!(variable.mean.welch.difference, 0.0, "mean on variable {v}");
            assert_eq!(variable.paired.welch.difference, 0.0, "paired on {v}");
            assert_eq!(variable.variance.welch.difference, 0.0, "variance on {v}");
        }
    }
    assert_eq!(report.structure.frobenius.cross, 0.0);
    assert_eq!(report.structure.frobenius.p_value(), 1.0);

    // The TOST *verdicts* are deliberately not asserted. At four seeds the
    // confidence interval is wider than the margin however small the
    // difference, so TOST cannot declare equivalence even here. That is the
    // sample size speaking, not the data.
    println!(
        "a source is identical to itself on all {VARIABLES} variables; TOST \
         power at {} seeds is {:.2} x the margin, so its verdicts are not \
         asserted",
        battery.seeds,
        report.worst_mean_power()
    );
    assert!(
        report.worst_mean_power() > 1.0,
        "four seeds now have enough power to declare equivalence, so this \
         test's reasoning about why the verdicts are unasserted is stale"
    );

    // The one gate that does apply at any size: a variable the reference never
    // moved must not move in the candidate.
    assert!(
        report.failures().is_empty(),
        "identical runs produced failures: {:?}",
        report.failures()
    );
}

/// And it must *reject* a source that is genuinely different.
///
/// This is the teeth check, and without it the previous test would be
/// satisfied by a battery that says yes to everything. The candidate here is
/// the port with its measurements shifted by a fifth of a standard deviation,
/// which is twice the mean margin and should fail on the mean at least.
#[test]
fn the_battery_rejects_a_source_that_is_actually_different() {
    let mut oracle = Oracle::lock();
    let battery = Battery::SMOKE;
    let start = start(&mut oracle);

    let reference: Vec<_> = (0..battery.seeds)
        .map(|s| {
            run_fortran(
                &mut oracle,
                &start,
                Scenario::NOMINAL,
                seed(s),
                battery.hours,
            )
        })
        .collect();
    let mut shifted = reference.clone();

    // Variable 7, reactor pressure: shift it by a fifth of its own spread.
    let variable = 6;
    let spread = tepsim_stats::Summary::of(
        &reference
            .iter()
            .flat_map(|r| r.series(variable))
            .collect::<Vec<f64>>(),
    )
    .sd();
    let shift = 0.2 * spread;
    for run in &mut shifted {
        for sample in &mut run.samples {
            sample[variable] += shift;
        }
    }
    println!("reactor pressure sd {spread:.4}, shifted by {shift:.4}");

    let report = compare(Scenario::NOMINAL, &reference, &shifted);

    // The verdicts have no power at four seeds, so this asserts on the
    // measured quantities instead: the shifted variable's mean difference is
    // the shift, and every other variable's is exactly zero.
    let shifted_report = &report.variables[variable];
    println!(
        "XMEAS({}): mean difference {:.6}, expected {shift:.6}; KS {:.4}, \
         energy {:.6}",
        variable + 1,
        shifted_report.mean.welch.difference,
        shifted_report.ks.cross,
        shifted_report.energy.cross
    );
    assert!(
        (shifted_report.mean.welch.difference - shift).abs() < 1e-9,
        "the mean difference is {}, not the {shift} that was injected",
        shifted_report.mean.welch.difference
    );
    // And the paired test sees it exactly, because every seed's difference is
    // the shift.
    assert!(
        (shifted_report.paired.welch.difference - shift).abs() < 1e-9,
        "the paired difference is {}",
        shifted_report.paired.welch.difference
    );
    assert!(
        !shifted_report.paired.equivalent,
        "a shift of twice the margin was declared equivalent by the paired test"
    );
    assert!(
        shifted_report.ks.cross > 0.0,
        "the distributions are shifted apart and KS is {}",
        shifted_report.ks.cross
    );
    assert!(
        shifted_report.energy.cross > 0.0,
        "the distributions are shifted apart and the energy distance is {}",
        shifted_report.energy.cross
    );
    // The shift is a translation, so serial structure is untouched. That the
    // battery reports it as untouched is itself worth checking: a statistic
    // that moved here would be responding to the mean rather than to the
    // dynamics.
    // Not exactly zero: the autocorrelation centres its series, and centring a
    // translated series subtracts a different constant, so the last bits move.
    assert!(
        shifted_report.autocorrelation.cross < 1e-13,
        "a pure translation moved the autocorrelation by {}, far past the \
         re-centring's last bits",
        shifted_report.autocorrelation.cross
    );
    // A log ratio, so dimensionless, and its floor is the last bits of the
    // re-centring rather than anything to do with the variable's units.
    assert!(
        shifted_report.spectrum.cross < 1e-12,
        "a pure translation moved the spectrum by {}, so the segment mean is \
         not being removed",
        shifted_report.spectrum.cross
    );

    // Every other variable is untouched, so the battery localises the change
    // rather than smearing it.
    for other in &report.variables {
        if other.variable == variable || other.constant {
            continue;
        }
        assert_eq!(
            other.mean.welch.difference,
            0.0,
            "variable {} moved although only {} was shifted",
            other.variable + 1,
            variable + 1
        );
        assert_eq!(other.ks.cross, 0.0, "KS on variable {}", other.variable + 1);
        let spread = tepsim_stats::Summary::of(
            &reference
                .iter()
                .flat_map(|r| r.series(other.variable))
                .collect::<Vec<f64>>(),
        )
        .sd()
        .max(1.0);
        assert!(
            other.energy.cross.abs() < 1e-12 * spread,
            "energy on variable {}",
            other.variable + 1
        );
    }
}
