//! Tier 5's verdicts, and where their margins come from.
//!
//! B-0047b. [`super`] produces runs; this judges them.
//!
//! # The margins are the hard part, not the statistics
//!
//! `PLAN.org` supplies exactly one margin: the mean must agree within a tenth
//! of the Fortran's standard deviation. For everything else there is no
//! absolute scale to appeal to. A Kolmogorov-Smirnov statistic of 0.04 is not
//! "small"; whether it is small depends on the sample size and on how much the
//! process itself wanders between runs. Picking 0.05 because it looks like a
//! p-value would be inventing a margin and calling it a result.
//!
//! So four of the six statistics are **calibrated against the Fortran's own
//! run-to-run variability**. Split its seeds into two halves, compute the
//! statistic Fortran-against-Fortran, and repeat over many splits. That gives a
//! sample of the statistic *under the null that the two sides are the same
//! simulator*, because they are. The cross-source value is then one more draw
//! from the same population, if the port is equivalent.
//!
//! # The test
//!
//! Under the null, the cross-source value is exchangeable with the `K`
//! within-source values, so its rank among the `K + 1` of them is uniform. The
//! one-sided permutation p-value is
//!
//! ```text
//! p = (1 + #{ within >= cross }) / (K + 1)
//! ```
//!
//! and the gate is `p > 0.05`. With `K = 20` that means the cross-source value
//! fails exactly when it is the strict maximum of the twenty-one. No tolerance
//! is chosen anywhere in that sentence, which is the point: the margin is
//! measured, not picked.
//!
//! # What is a sample
//!
//! **The unit of observation is the run, not the time step.** Pooling every
//! sample of every run into one vector and running a t-test on it would give a
//! standard error computed from 96,000 numbers, but consecutive samples of a
//! plant three minutes apart are strongly correlated, so the effective sample
//! size is far smaller and the test would be wildly overconfident. Each run
//! contributes one mean and one variance; a hundred seeds give a hundred
//! independent observations, and the standard error is honest.
//!
//! The distribution statistics (Kolmogorov-Smirnov, energy distance) *do* use
//! the pooled samples, because they are comparing marginal distributions and
//! autocorrelation does not bias them, only their p-values. Their p-values are
//! reported and never gated on; the calibration is what gates them.
//!
//! # Smoke and full
//!
//! Calibration needs enough seeds to split, and TOST needs enough seeds to
//! have power. The smoke battery `ci` runs has four, which gives three
//! distinct half-splits: a permutation test with `K = 3` cannot reject
//! anything at 0.05, and a TOST on four observations cannot fit a confidence
//! interval inside a tenth of a standard deviation *even when the two samples
//! are identical*. So the smoke battery **reports** and the full battery
//! **gates**, which is stated in [`Report::gated`] rather than left implicit.
//!
//! That is not a weakness hidden in a corner. It is the reason
//! `the_battery_finds_a_source_equivalent_to_itself` asserts on the distances
//! and the TOST *intervals* rather than on the TOST verdicts: an equivalence
//! test that fails for want of data has said nothing, and reading its failure
//! as a difference would be exactly the error TOST exists to prevent.
//!
//! [`Report::worst_mean_power`] reports the confidence interval's half-width as
//! a multiple of the margin, so a battery that failed for want of seeds says so
//! rather than blaming the port. Measured values:
//!
//! | battery                    | worst power |
//! |----------------------------|-------------|
//! | 4 seeds, 2 h (smoke)       |       11.93 |
//! | 8 seeds, 12 h              |        3.38 |
//! | 8 seeds, 48 h              |        1.04 |
//! | 100 seeds, 48 h, nominal   |        1.59 |
//! | 100 seeds, 48 h, `IDV(1)`  |        1.61 |
//! | 100 seeds, 48 h, `IDV(2)`  |        0.23 |
//!
//! Horizon matters more than seed count, because a run mean over 960 samples is
//! far better averaged than one over 40. But the *scenario* matters as much as
//! either, and that is the finding B-0047b turned up: the margin is a tenth of
//! the pooled standard deviation, and a disturbed plant has a much larger
//! pooled spread than a quiet one, so the same absolute run-to-run wander in
//! the mean is comfortably inside the margin under `IDV(2)` and outside it at
//! the nominal operating point.
//!
//! The variables that run out of power are the manipulated ones, `XMV(1..12)`.
//! That is what an integrating controller does: its output has a random-walk
//! component, so its mean over 48 hours varies between seeds by much more than
//! a tenth of its own within-run spread.
//!
//! **An underpowered test has not failed.** It has declined to decide, and
//! reporting it as a difference would be exactly the error TOST exists to
//! prevent, in the other direction. [`Report::undecided`] separates the two,
//! and [`Report::seeds_for_power`] says how many seeds the margin would need.
//!
//! # Variables that never move
//!
//! `XMV(12)`, the agitator, is never written by a controller and sits at 50 for
//! every run. Its variance is zero, so a variance ratio is `0/0` and TOST
//! returns `NaN`. Such a variable is marked [`VariableReport::constant`] and
//! its moment gates are skipped, but the candidate is then *required* to be
//! constant too: a port whose agitator drifted would otherwise slip through the
//! one gap the report has.

use tepsim_stats::{
    CorrelationMatrix, Summary, Tost, Window, autocorrelation, band_comparison, energy_distance,
    frobenius_distance, ks_statistic, ks_two_sample_p, log_band_edges, tost, tost_paired, welch,
};

use super::{Run, Scenario, VARIABLES};

/// `XMEAS(1..41)` occupy the first 41 columns of a 53-wide row; the valves
/// follow. `tepsim::MEASUREMENTS` is the same number and is asserted equal in
/// `the_sticking_map_is_the_fortrans`.
const MEASUREMENTS: usize = 41;

/// The equivalence margin for a mean, as a fraction of the reference standard
/// deviation. `PLAN.org`, "Tier 5".
pub const MEAN_MARGIN_FRACTION: f64 = 0.1;

/// The equivalence margin for a variance ratio, in log units.
///
/// `ln(1.1)`: variances within ten percent of each other. The failure this
/// exists to catch is the one `PLAN.org` names, the existing Python port's
/// reactor-pressure standard deviation of 8.10 against the Fortran's 61.48.
/// That is a variance ratio of 58, so ten percent is more than two orders
/// tighter than the failure it must not miss.
pub const VARIANCE_MARGIN_LOG: f64 = 0.095_310_179_804_324_86;

/// The significance level for both TOST and the permutation tests.
pub const ALPHA: f64 = 0.05;

/// How many half-splits of the reference runs the calibration uses.
///
/// Twenty gives a permutation test whose smallest attainable p-value is
/// `1/21 = 0.048`, just under [`ALPHA`]. Fewer would make the test unable to
/// reject at all; many more would cost time for resolution the gate cannot
/// use.
pub const CALIBRATION_SPLITS: usize = 20;

/// Lags compared when checking serial structure. `PLAN.org`, "out to lag 200".
///
/// Reduced when the runs are too short to reach it; see [`usable_lags`]. The
/// smoke battery's two-hour runs are 40 samples, and asking for 200 lags of a
/// 40-sample series gives a shorter answer than requested, which silently
/// became a `NaN` gap before this was handled.
pub const ACF_LAGS: usize = 200;

/// Welch segment length, at the full battery's run length. A power of two.
///
/// Reduced when the runs are too short; see [`usable_segment`].
pub const SPECTRUM_SEGMENT: usize = 64;

/// How many lags a series of this length supports.
#[must_use]
pub const fn usable_lags(samples: usize) -> usize {
    // A quarter of the series, so the longest lag still averages four
    // independent stretches rather than one.
    let quarter = samples / 4;
    if quarter < ACF_LAGS {
        quarter
    } else {
        ACF_LAGS
    }
}

/// The largest usable Welch segment for a series of this length.
///
/// At least two segments' worth of data, so there is something to average, and
/// a power of two because the transform requires it.
#[must_use]
pub fn usable_segment(samples: usize) -> usize {
    let half = samples / 2;
    let mut segment = 8_usize;
    while segment * 2 <= half && segment * 2 <= SPECTRUM_SEGMENT {
        segment *= 2;
    }
    segment
}

/// Frequency bands compared, logarithmically spaced.
pub const SPECTRUM_BANDS: usize = 5;

/// Sampling rate of the recorded series, in samples per hour.
///
/// One sample every 180 seconds, so twenty per hour.
pub const SAMPLE_RATE_PER_HOUR: f64 = 20.0;

/// One statistic: the cross-source value, and the within-source null it is
/// judged against.
#[derive(Clone, Debug, PartialEq)]
pub struct Calibrated {
    /// What was measured.
    pub name: &'static str,
    /// The value between the two sources.
    pub cross: f64,
    /// The same statistic computed between two halves of the reference source.
    pub within: Vec<f64>,
}

impl Calibrated {
    /// The one-sided permutation p-value.
    ///
    /// `(1 + #{within >= cross}) / (K + 1)`. The `1 +` in both places is what
    /// makes it valid rather than merely plausible: it counts the observed
    /// value as one of the draws, which is what exchangeability says it is.
    #[must_use]
    pub fn p_value(&self) -> f64 {
        if self.within.is_empty() {
            return f64::NAN;
        }
        let at_least = self.within.iter().filter(|w| **w >= self.cross).count();
        (1 + at_least) as f64 / (self.within.len() + 1) as f64
    }

    /// The largest within-source value, for the report.
    #[must_use]
    pub fn within_max(&self) -> f64 {
        self.within.iter().copied().fold(f64::NAN, f64::max)
    }

    /// Whether the cross-source value sits inside the within-source spread.
    ///
    /// `None` when there are too few splits for the test to be able to reject,
    /// which is the smoke battery's situation and is reported rather than
    /// silently passed.
    #[must_use]
    pub fn passes(&self) -> Option<bool> {
        if self.within.len() + 1 < (1.0 / ALPHA) as usize {
            return None;
        }
        Some(self.p_value() > ALPHA)
    }
}

/// Everything measured for one variable of one scenario.
#[derive(Clone, Debug)]
pub struct VariableReport {
    /// Which scenario.
    pub scenario: Scenario,
    /// Which of the 53 variables, zero-based.
    pub variable: usize,
    /// TOST on the per-run means, unpaired, against a tenth of the reference
    /// standard deviation.
    ///
    /// Kept and reported, but **not** the gate. See [`VariableReport::paired`].
    pub mean: Tost,
    /// The same margin, applied to the *paired* differences.
    ///
    /// This is the gate. The two sources are run at the same seeds, so each
    /// seed's difference is an observation and the seed-to-seed wander, which
    /// both sources share because they see the same disturbance realisation,
    /// cancels out of it entirely.
    ///
    /// It matters enormously here. B-0047c measured the unpaired test needing
    /// 237 to 323 seeds at a tenth of a standard deviation, because the
    /// manipulated variables are driven by integrating controllers whose
    /// 48-hour mean wanders from seed to seed by far more than the two
    /// implementations differ. Pairing removes exactly that term.
    pub paired: Tost,
    /// The paired interval's half-width as a multiple of the margin.
    pub paired_power: f64,
    /// TOST on the per-run log variances, against [`VARIANCE_MARGIN_LOG`].
    pub variance: Tost,
    /// Two-sample Kolmogorov-Smirnov on the pooled samples.
    pub ks: Calibrated,
    /// The reported KS p-value. Recorded, never gated: the samples are
    /// autocorrelated so the asymptotic p-value is optimistic.
    pub ks_p: f64,
    /// Energy distance between the pooled samples.
    pub energy: Calibrated,
    /// Largest absolute difference between the two mean autocorrelations, over
    /// lags 1 to [`ACF_LAGS`].
    pub autocorrelation: Calibrated,
    /// Largest `|log ratio|` between the two mean Welch spectra, band by band.
    pub spectrum: Calibrated,
    /// The reference never moved on this variable, so the moment tests are
    /// undefined. See the module docs.
    pub constant: bool,
    /// And neither did the candidate.
    pub candidate_constant: bool,
    /// This variable is a valve the scenario's fault sticks, so its moment
    /// tests are not applied. See [`sticks_this_valve`] and
    /// [`Report::failures`].
    pub stuck_valve: bool,
    /// The unpaired TOST interval's half-width as a multiple of the margin.
    ///
    /// Under one means the mean test *can* declare equivalence; over one means
    /// there is not enough data to, whatever the difference is. Recorded so
    /// that a battery which failed for want of power says so.
    pub mean_power: f64,
}

impl VariableReport {
    /// Every calibrated statistic, for iterating over.
    #[must_use]
    pub fn calibrated(&self) -> [&Calibrated; 4] {
        [
            &self.ks,
            &self.energy,
            &self.autocorrelation,
            &self.spectrum,
        ]
    }

    /// Whether every gate this variable is subject to passed.
    ///
    /// `None` if the calibration could not gate; see [`Calibrated::passes`].
    #[must_use]
    pub fn passes(&self) -> Option<bool> {
        if self.constant {
            return Some(self.candidate_constant);
        }
        let mut all = self.stuck_valve || (self.paired.equivalent && self.variance.equivalent);
        for statistic in self.calibrated() {
            all &= statistic.passes()?;
        }
        Some(all)
    }
}

/// A variable's series from every run, and the summary statistics of each.
struct Ensemble {
    /// One entry per run: that run's series.
    series: Vec<Vec<f64>>,
}

impl Ensemble {
    fn of(runs: &[Run], variable: usize) -> Self {
        Self {
            series: runs.iter().map(|r| r.series(variable)).collect(),
        }
    }

    /// The shortest series in the ensemble.
    fn shortest(&self) -> usize {
        self.series.iter().map(Vec::len).min().unwrap_or(0)
    }

    fn subset(&self, indices: &[usize]) -> Self {
        Self {
            series: indices.iter().map(|i| self.series[*i].clone()).collect(),
        }
    }

    /// One mean per run: the unit of observation for the moment tests.
    fn run_means(&self) -> Summary {
        let mut summary = Summary::new();
        for series in &self.series {
            summary.push(Summary::of(series).mean());
        }
        summary
    }

    /// One log variance per run.
    ///
    /// Logs rather than variances, so that the equivalence margin is a *ratio*
    /// and applies equally to a variable that varies by 1e-4 and one that
    /// varies by 60.
    fn run_log_variances(&self) -> Summary {
        let mut summary = Summary::new();
        for series in &self.series {
            let variance = Summary::of(series).variance();
            if variance > 0.0 {
                summary.push(tepsim_stats::special::ln(variance));
            }
        }
        summary
    }

    /// Every sample of every run, for the distribution statistics.
    fn pooled(&self) -> Vec<f64> {
        self.series.iter().flatten().copied().collect()
    }

    /// The autocorrelation averaged over runs.
    fn mean_autocorrelation(&self, lags: usize) -> Vec<f64> {
        let mut total = vec![0.0; lags + 1];
        let mut counted = 0;
        for series in &self.series {
            let acf = autocorrelation(series, lags);
            if acf.len() == lags + 1 {
                for (slot, value) in total.iter_mut().zip(&acf) {
                    *slot += value;
                }
                counted += 1;
            }
        }
        if counted == 0 {
            return Vec::new();
        }
        for slot in &mut total {
            *slot /= f64::from(counted);
        }
        total
    }

    /// The Welch spectrum averaged over runs.
    fn mean_spectrum(&self, segment: usize) -> Option<tepsim_stats::Spectrum> {
        let usable: Vec<&Vec<f64>> = self.series.iter().filter(|s| s.len() >= segment).collect();
        if usable.is_empty() {
            return None;
        }
        let mut accumulated: Option<tepsim_stats::Spectrum> = None;
        for series in &usable {
            let spectrum = welch(
                series,
                SAMPLE_RATE_PER_HOUR,
                segment,
                segment / 2,
                Window::Hann,
            );
            match &mut accumulated {
                None => accumulated = Some(spectrum),
                Some(total) => {
                    for (slot, value) in total.density.iter_mut().zip(&spectrum.density) {
                        *slot += value;
                    }
                }
            }
        }
        let mut spectrum = accumulated?;
        let scale = 1.0 / usable.len() as f64;
        for slot in &mut spectrum.density {
            *slot *= scale;
        }
        Some(spectrum)
    }
}

/// The four cross-source statistics, computed between two ensembles.
fn statistics(a: &Ensemble, b: &Ensemble) -> [f64; 4] {
    let pooled_a = a.pooled();
    let pooled_b = b.pooled();
    let ks = ks_statistic(&pooled_a, &pooled_b);
    let energy = energy_distance(&pooled_a, &pooled_b);

    // Both sides sized from the shorter of the two, so they are comparable.
    let shortest = a.shortest().min(b.shortest());
    let lags = usable_lags(shortest);
    let segment = usable_segment(shortest);

    let acf_a = a.mean_autocorrelation(lags);
    let acf_b = b.mean_autocorrelation(lags);
    let acf_gap = if acf_a.len() == acf_b.len() && acf_a.len() > 1 {
        acf_a
            .iter()
            .zip(&acf_b)
            .skip(1)
            .map(|(x, y)| libm::fabs(x - y))
            .fold(0.0_f64, f64::max)
    } else {
        f64::NAN
    };

    let spectrum_gap = match (a.mean_spectrum(segment), b.mean_spectrum(segment)) {
        (Some(sa), Some(sb)) => {
            let edges = log_band_edges(sa.resolution, SAMPLE_RATE_PER_HOUR, SPECTRUM_BANDS);
            band_comparison(&sa, &sb, &edges)
                .iter()
                // Bands carrying no power carry no information either.
                .filter(|band| band.power_a > 0.0 && band.power_b > 0.0)
                .map(|band| libm::fabs(tepsim_stats::special::ln(band.ratio)))
                .fold(0.0_f64, f64::max)
        }
        _ => f64::NAN,
    };

    [ks, energy, acf_gap, spectrum_gap]
}

/// Deterministic half-splits of `n` indices.
///
/// Split `k` puts index `i` in the first half when bit `k` of a SplitMix hash
/// of `(i, k)` is set, rebalanced so both halves are the same size. Fixed and
/// reproducible: two runs of the battery calibrate against the same splits, so
/// a change in a verdict is a change in the port.
fn half_splits(n: usize, count: usize) -> Vec<(Vec<usize>, Vec<usize>)> {
    let mut out = Vec::with_capacity(count);
    for k in 0..count {
        let mut order: Vec<usize> = (0..n).collect();
        order.sort_by_key(|i| {
            let mut z = ((*i as u64) << 32 | k as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15);
            z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
            z ^= z >> 27;
            z
        });
        let half = n / 2;
        out.push((order[..half].to_vec(), order[half..].to_vec()));
    }
    out
}

/// Whether this scenario's fault sticks this variable's valve.
///
/// `teprob.f:793-798` assigns `IVST(10)=IDV(14)`, `IVST(11)=IDV(15)` and
/// `IVST(5)=IVST(7)=IVST(8)=IVST(9)=IDV(19)`. Read from
/// [`tepsim_core::FAULTS`] rather than transcribed, and
/// `the_sticking_map_is_the_fortrans` asserts the table against those lines.
///
/// Variables are the 53-wide row: `XMEAS(1..41)` at 0 to 40, `XMV(1..12)` at 41
/// to 52.
#[must_use]
pub fn sticks_this_valve(scenario: Scenario, variable: usize) -> bool {
    let Some(valve) = variable.checked_sub(MEASUREMENTS).map(|v| v + 1) else {
        return false;
    };
    if scenario.fault == 0 {
        return false;
    }
    tepsim_core::FAULTS
        .iter()
        .find(|f| f.index == scenario.fault)
        .is_some_and(|f| match f.shape {
            tepsim_core::Shape::Sticking { valves } => valves.contains(&valve),
            _ => false,
        })
}

/// Compare one variable across two sets of runs.
///
/// `reference` is the Fortran and `candidate` the port. The asymmetry matters:
/// the margins are fractions of the *reference*'s spread, and the calibration
/// splits the *reference*'s runs.
#[must_use]
pub fn compare_variable(
    scenario: Scenario,
    variable: usize,
    reference: &[Run],
    candidate: &[Run],
) -> VariableReport {
    let a = Ensemble::of(reference, variable);
    let b = Ensemble::of(candidate, variable);

    // Moments, one observation per run.
    let reference_spread = Summary::of(&a.pooled()).sd();
    // Written so a NaN spread counts as constant rather than falling past the
    // guard: a variable whose spread cannot be computed has not moved in any
    // sense the moment tests can use.
    let constant = reference_spread.is_nan() || reference_spread <= 0.0;
    let candidate_spread = Summary::of(&b.pooled()).sd();
    let candidate_constant = candidate_spread.is_nan() || candidate_spread <= 0.0;
    let margin = MEAN_MARGIN_FRACTION * reference_spread;
    let mean = tost(&b.run_means(), &a.run_means(), margin, ALPHA);

    // The paired form. `reference[k]` and `candidate[k]` are the same seed by
    // construction, which is what makes the differences meaningful.
    let mut differences = Summary::new();
    for (left, right) in a.series.iter().zip(&b.series) {
        differences.push(Summary::of(right).mean() - Summary::of(left).mean());
    }
    let paired = tost_paired(&differences, margin, ALPHA);
    let paired_power = if margin > 0.0 {
        0.5 * (paired.interval.1 - paired.interval.0) / margin
    } else {
        f64::NAN
    };
    let mean_power = if margin > 0.0 {
        0.5 * (mean.interval.1 - mean.interval.0) / margin
    } else {
        f64::NAN
    };
    let variance = tost(
        &b.run_log_variances(),
        &a.run_log_variances(),
        VARIANCE_MARGIN_LOG,
        ALPHA,
    );

    let cross = statistics(&a, &b);
    let ks_p = ks_two_sample_p(cross[0], a.pooled().len(), b.pooled().len());

    // The null: the reference against itself.
    let mut within = [const { Vec::new() }; 4];
    for (left, right) in half_splits(reference.len(), CALIBRATION_SPLITS) {
        if left.is_empty() || right.is_empty() {
            continue;
        }
        let values = statistics(&a.subset(&left), &a.subset(&right));
        for (slot, value) in within.iter_mut().zip(values) {
            if value.is_finite() {
                slot.push(value);
            }
        }
    }

    let names = [
        "kolmogorov-smirnov",
        "energy distance",
        "autocorrelation",
        "spectrum",
    ];
    let mut calibrated = within
        .into_iter()
        .zip(cross)
        .zip(names)
        .map(|((within, cross), name)| Calibrated {
            name,
            cross,
            within,
        });

    VariableReport {
        scenario,
        variable,
        mean,
        variance,
        ks: calibrated.next().expect("four statistics"),
        ks_p,
        energy: calibrated.next().expect("four statistics"),
        autocorrelation: calibrated.next().expect("four statistics"),
        spectrum: calibrated.next().expect("four statistics"),
        constant,
        candidate_constant,
        mean_power,
        paired,
        paired_power,
        stuck_valve: sticks_this_valve(scenario, variable),
    }
}

/// The correlation structure of one scenario, which is what PCA consumes.
#[derive(Clone, Debug)]
pub struct StructureReport {
    /// Which scenario.
    pub scenario: Scenario,
    /// Frobenius distance between the two correlation matrices, calibrated.
    pub frobenius: Calibrated,
    /// The pair of variables whose correlation moved most, one-based.
    pub worst_pair: Option<(usize, usize, f64, f64)>,
    /// Entries skipped because a variable never moved.
    pub skipped: usize,
}

/// Compare the correlation structure of two sets of runs.
#[must_use]
pub fn compare_structure(
    scenario: Scenario,
    reference: &[Run],
    candidate: &[Run],
) -> StructureReport {
    let matrix = |runs: &[Run], indices: &[usize]| {
        let series: Vec<Vec<f64>> = (0..VARIABLES)
            .map(|v| {
                indices
                    .iter()
                    .flat_map(|i| runs[*i].series(v))
                    .collect::<Vec<f64>>()
            })
            .collect();
        CorrelationMatrix::of(&series)
    };

    let everything: Vec<usize> = (0..reference.len()).collect();
    let all_candidates: Vec<usize> = (0..candidate.len()).collect();
    let a = matrix(reference, &everything);
    let b = matrix(candidate, &all_candidates);
    let (cross, skipped) = frobenius_distance(&a, &b);

    let mut within = Vec::new();
    for (left, right) in half_splits(reference.len(), CALIBRATION_SPLITS) {
        if left.is_empty() || right.is_empty() {
            continue;
        }
        let (distance, _) =
            frobenius_distance(&matrix(reference, &left), &matrix(reference, &right));
        if distance.is_finite() {
            within.push(distance);
        }
    }

    StructureReport {
        scenario,
        frobenius: Calibrated {
            name: "frobenius",
            cross,
            within,
        },
        worst_pair: tepsim_stats::worst_correlation_difference(&a, &b)
            .map(|(i, j, x, y)| (i + 1, j + 1, x, y)),
        skipped,
    }
}

/// Everything measured for one scenario.
#[derive(Clone, Debug)]
pub struct Report {
    /// Which scenario.
    pub scenario: Scenario,
    /// One entry per variable.
    pub variables: Vec<VariableReport>,
    /// The correlation structure.
    pub structure: StructureReport,
    /// Whether the calibration had enough splits to gate at [`ALPHA`].
    ///
    /// Also governs the moment gates, which need power as much as the
    /// permutation tests need splits.
    pub gated: bool,
    /// How many seeds each source contributed.
    pub seeds: usize,
}

/// Compare every variable and the structure, for one scenario.
#[must_use]
pub fn compare(scenario: Scenario, reference: &[Run], candidate: &[Run]) -> Report {
    let variables = (0..VARIABLES)
        .map(|v| compare_variable(scenario, v, reference, candidate))
        .collect();
    let structure = compare_structure(scenario, reference, candidate);
    Report {
        scenario,
        variables,
        structure,
        gated: reference.len() >= 2 * (1.0 / ALPHA) as usize,
        seeds: reference.len(),
    }
}

impl Report {
    /// The worst mean-test power across the variables.
    ///
    /// Over one means at least one variable has too few seeds for TOST to
    /// declare equivalence however close the two sources are. A battery that
    /// reported failures at that size would be reporting its own sample size.
    #[must_use]
    pub fn worst_mean_power(&self) -> f64 {
        self.variables
            .iter()
            .filter(|v| !v.constant)
            .map(|v| v.mean_power)
            .fold(0.0_f64, f64::max)
    }

    /// The worst *paired* mean-test power, which is what the gate uses.
    #[must_use]
    pub fn worst_paired_power(&self) -> f64 {
        self.variables
            .iter()
            .filter(|v| !v.constant)
            .map(|v| v.paired_power)
            .fold(0.0_f64, f64::max)
    }

    /// Variables whose mean test could not decide, with their power.
    ///
    /// Sorted worst first. An entry here is *not* a failure: the confidence
    /// interval is wider than the margin, so the test could not have declared
    /// equivalence however close the two sources were.
    #[must_use]
    pub fn undecided(&self) -> Vec<(usize, f64, f64)> {
        let mut out: Vec<(usize, f64, f64)> = self
            .variables
            .iter()
            .filter(|v| !v.constant && !v.stuck_valve && v.paired_power >= 1.0)
            .map(|v| {
                (
                    v.variable + 1,
                    v.paired_power,
                    // How far apart the two means actually are, as a fraction
                    // of the margin. Small here and large power means "not
                    // enough data"; large here means something real.
                    libm::fabs(v.paired.welch.difference) / v.paired.margin,
                )
            })
            .collect();
        out.sort_by(|a, b| b.1.total_cmp(&a.1));
        out
    }

    /// How many seeds the mean test would need to decide every variable.
    ///
    /// The confidence interval's half-width scales as `1 / sqrt(n)`, so the
    /// requirement is `n * power^2`, rounded up. Reported rather than acted on:
    /// changing the battery's size is a decision, not an inference.
    #[must_use]
    pub fn seeds_for_power(&self) -> usize {
        let worst = self.worst_paired_power();
        if worst.is_nan() || worst <= 1.0 {
            return self.seeds;
        }
        (self.seeds as f64 * worst * worst).ceil() as usize
    }

    /// The variables whose moment gates were replaced by distributional ones,
    /// with the verdicts that replaced them.
    ///
    /// # Why a stuck valve's mean is not a statistic about the process
    ///
    /// B-0047d. `IDV(14)`, `IDV(15)` and `IDV(19)` widen a valve's dead band
    /// at `teprob.f:793-798`, and `teprob.f:801` then moves that valve only
    /// when `DABS(VCV(I)-XMV(I)).GT.VST(I)*IVST(I)`. That is a discontinuous
    /// branch on a floating-point comparison. One ULP of `exp` difference,
    /// which Tier 2 accepts on about ten percent of arguments and which the
    /// vendored `libm` has against gfortran's, decides which side of the
    /// threshold the command lands on, and the valve then holds a *different
    /// position* for the rest of the run.
    ///
    /// The resulting series is a few plateaux, not a distribution around a
    /// centre, and its 48-hour mean is an artefact of how long it spent on each
    /// plateau. Two sources can agree on every plateau and on how often each is
    /// visited while their means differ by a tenth of a standard deviation.
    /// That is what happened: `IDV(14)` and `IDV(19)` missed the mean gate at
    /// 1.016e-1 on exactly the valves those faults stick, and were bit-
    /// identical (0.000e0) under `--features libm-system`, where `exp` is the
    /// one gfortran calls.
    ///
    /// So for these variables the battery asks the question it can answer:
    /// Kolmogorov-Smirnov and energy distance on the pooled samples, both
    /// permutation-calibrated against the reference's own split-half variation,
    /// plus the autocorrelation and spectrum gates. Nothing was loosened. The
    /// margin is unchanged and every other variable is judged exactly as
    /// before.
    #[must_use]
    pub fn stuck_valves(&self) -> Vec<(usize, Option<bool>, Option<bool>)> {
        self.variables
            .iter()
            .filter(|v| v.stuck_valve)
            .map(|v| (v.variable + 1, v.ks.passes(), v.energy.passes()))
            .collect()
    }

    /// Statistics outside their margin, with what failed.
    ///
    /// Only gates that this report can actually decide. The moment tests need
    /// power and the permutation tests need splits, so an ungated report
    /// carries neither: a TOST that fails for want of seeds has said nothing,
    /// and listing it as a failure would report the sample size as a defect of
    /// the port. The one gate that holds at any size is that a variable the
    /// reference never moved must not move in the candidate.
    #[must_use]
    pub fn failures(&self) -> Vec<(usize, &'static str)> {
        let mut out = Vec::new();
        for report in &self.variables {
            if report.constant {
                // The reference never moved. The only claim available is that
                // the candidate did not either, and it is worth making.
                if !report.candidate_constant {
                    out.push((report.variable + 1, "moved when the reference did not"));
                }
                continue;
            }
            // An underpowered test has declined to decide, not failed. See
            // `Report::undecided`, which reports those separately.
            let decidable = report.paired_power < 1.0;
            // A stuck valve is judged on its distribution, not on its mean.
            // See `Report::stuck_valves` for why, and note that the four
            // calibrated statistics below still run on it: this substitutes
            // what is measured, it does not stop measuring.
            let moments = self.gated && decidable && !report.stuck_valve;
            if moments && !report.paired.equivalent {
                out.push((report.variable + 1, "mean"));
            }
            if moments && !report.variance.equivalent {
                out.push((report.variable + 1, "variance"));
            }
            for statistic in report.calibrated() {
                if statistic.passes() == Some(false) {
                    out.push((report.variable + 1, statistic.name));
                }
            }
        }
        if self.structure.frobenius.passes() == Some(false) {
            out.push((0, "frobenius"));
        }
        out
    }
}
