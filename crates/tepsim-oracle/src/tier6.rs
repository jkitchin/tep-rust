//! Tier 6: can a fault detector tell which simulator generated its data?
//!
//! B-0050. `PLAN.org` calls this "the operational definition of practically
//! equivalent for this project's research audience", and the reason is that it
//! is the only tier whose units are the ones a reader of the TEP literature
//! already has intuitions about: detection rate, false alarm rate, detection
//! delay. A Frobenius distance of 0.03 between two correlation matrices means
//! nothing to anyone. "The detector caught 97.4% of `IDV(5)` on the Fortran
//! and 97.5% on the port, and two halves of the Fortran differ by 0.4%" means
//! something to everyone.
//!
//! # The experiment
//!
//! Four combinations, per fault:
//!
//! | trained on | evaluated on | what it is |
//! |---|---|---|
//! | Fortran | Fortran | within-source |
//! | Fortran | port    | **cross-source** |
//! | port    | port    | within-source |
//! | port    | Fortran | **cross-source** |
//!
//! The claim is not that the cross-source numbers equal the within-source
//! numbers. Nothing equals anything in a stochastic simulation. The claim is
//! that the cross-source difference is **not larger than the difference the
//! reference simulator shows against itself**, and that requires estimating
//! the second quantity rather than asserting it.
//!
//! # Where the null comes from
//!
//! The same construction [`crate::tier5::battery`] uses for its margins: split
//! the reference's test seeds into two halves and compute the statistic
//! Fortran-against-Fortran. That is a draw of the statistic *under the null
//! that the two sides are the same simulator*, because they are. Repeat over
//! [`CALIBRATION_SPLITS`] distinct splits and the cross-source value is one
//! more draw from the same population, if the port is equivalent.
//! [`Calibrated`] then supplies the one-sided permutation p-value
//! `(1 + #{within >= cross}) / (K + 1)`, and the gate is `p > ALPHA`.
//!
//! No margin is chosen anywhere in that sentence. That is the point of the
//! design, and it is why a detection rate that differs by 0.2 between the two
//! sources can be a pass while one that differs by 0.02 can be a failure: what
//! matters is how much the reference's own seeds move it.
//!
//! # Two corrections the naive version needs
//!
//! **Group size.** A half-split compares `n/2` runs against `n/2`, while the
//! cross-source comparison has `n` on each side. A difference of means between
//! smaller groups is larger simply because it is noisier, so comparing a raw
//! `|mean - mean|` across the two would make the null about `sqrt(2)` too wide
//! and bias every verdict toward passing. [`separation`] removes that by
//! dividing by `sqrt(1/|A| + 1/|B|)`, which makes the statistic's null
//! expectation depend only on the per-run variance and not on the group sizes.
//! The raw difference is still reported, because that is the number a reader
//! wants; the *test* runs on the normalised one.
//!
//! **Pairing.** Both sources are run at the same generator words, so run `s` of
//! one is the partner of run `s` of the other. If those partners stayed
//! correlated, the cross-source difference would be a paired difference with
//! less variance than the unpaired null implies, and the test would again be
//! biased toward passing. Whether they do is a fact about the process and not
//! something to assume, so [`MetricComparison::pairing`] measures it: the
//! Pearson correlation between the two sources' per-run values across seeds. A
//! value near zero says the trajectories have separated by the time the metric
//! is taken and the unpaired null is the right one. A large value is a caveat
//! that has to be reported with the verdict.
//!
//! # The detector, fixed in advance
//!
//! Every parameter below is stated here, with its reason, and none of them was
//! chosen after seeing a number. That ordering is the whole methodological
//! content of a detector comparison: a confidence level picked once the
//! detection rates are on the screen is not a parameter, it is a result.
//!
//! - **Variables**: all 53 recorded, `XMEAS(1..41)` then `XMV(1..12)`. The TEP
//!   literature uses 52, excluding `XMV(12)`, the agitator, because the base
//!   control scheme never moves it. [`tepsim_stats::Pca`] standardises a
//!   constant column to exactly zero, so it contributes nothing to either
//!   statistic and the fitted model is numerically the literature's 52-variable
//!   one. Keeping the column makes the exclusion visible in
//!   `Pca::constant_columns` rather than hidden in an index list.
//! - **Detector**: PCA on the correlation matrix, with Hotelling's T-squared
//!   and the squared prediction error. `PLAN.org` names it first, and it is the
//!   baseline of the TEP fault-detection literature.
//! - **Lags**: two detectors, at [`LAG_COUNTS`] = 0 and 2. Zero is static PCA;
//!   two is dynamic PCA at the lag count Ku et al. and essentially all the TEP
//!   literature use. Two detectors rather than one because a claim that rests
//!   on a single detector is a claim about that detector.
//! - **Retention**: [`RETENTION`], the smallest number of components explaining
//!   90% of the variance. A *rule*, applied identically to whichever training
//!   set the detector is handed, so the component count is an output of the
//!   experiment. If the two sources produce different component counts, that is
//!   itself a finding and [`DetectorReport::components`] reports it.
//! - **Confidence**: [`CONFIDENCE`] = 0.99 for both limits, the level the TEP
//!   literature quotes.
//! - **Persistence**: [`PERSISTENCE`] = 3 consecutive alarms for a detection
//!   delay. At a 1% nominal false alarm rate three in a row happen by chance
//!   once in a million samples, so a delay found this way is not luck. One
//!   alarm would mostly measure luck; six is also used in the literature and
//!   would give different numbers. The detection *rate* uses no persistence, so
//!   that `1 - FDR` stays the missed detection rate the published tables report.
//! - **Training seeds are disjoint from test seeds.** The T-squared limit in
//!   [`tepsim_stats::t_squared_limit`] is the one for an observation
//!   independent of the training set, so the evaluation has to actually be
//!   independent, or the false alarm rate is optimistic by construction.
//!
//! # What the runs are, and what they are not
//!
//! [`crate::tier5::run_fortran`] and [`crate::tier5::run_port`] switch the
//! requested `IDV` on at `t = 0`, so a faulted run is faulted throughout and
//! its whole record is post-fault: the detection rate is taken with
//! `onset = 0`. There is no fault-free prologue to measure false alarms on, so
//! the false alarm rate is measured on **separate fault-free runs** at the test
//! seeds, which is the call `detection::false_alarm_rate` documents for
//! `d00_te`.
//!
//! The driver also forces `IDV(12)` on at eight hours, in both sources, exactly
//! as `temain_mod.f:366-368` does. So a 48-hour "fault-free" run carries a
//! random-variation disturbance for forty of its forty-eight hours. That is not
//! a defect of this experiment: it is how the published `d00`-`d21` files were
//! generated, since they came from this driver. It does mean the absolute
//! detection rates here are not comparable with published tables that assume a
//! quiet training set; comparing against published rates is B-0049a's job, and
//! nothing in the cross-source claim depends on it, because both sources carry
//! the same `IDV(12)`.
//!
//! # Battery size
//!
//! Selected by `TEP_TIER6`, the way `TEP_TIER5` and `TEP_TIER4_HOURS` are. The
//! smoke battery runs in seconds and **cannot conclude anything**: with six
//! test seeds there are only ten distinct half-splits, so the permutation test
//! cannot produce a p-value below 0.05 whatever the data say, and
//! [`Calibrated::passes`] returns `None` rather than a verdict. That is not a
//! weakness hidden in a corner; it falls out of the arithmetic and is printed
//! on every run. See [`Battery`].
//!
//! # Why this tier also measures how close the alarms came to flipping
//!
//! A permutation test on a statistic that is exactly zero passes, and passing
//! that way is uninformative on its own. It happens here, and the reason is a
//! result rather than an accident: B-0044 measured the closed loop's worst
//! `XMEAS` disagreement at 8.7e-14 after four hours on the vendored `libm`,
//! because the controllers pull the plant back to a setpoint and stop a
//! one-ULP `exp` difference from compounding. Tier 4's *open*-loop divergence
//! is the number people remember, and it does not apply to a closed-loop run.
//!
//! So "the detection rates were identical" is not a claim this tier is willing
//! to leave at that. [`Agreement`] measures the two quantities that decide
//! whether an alarm *could* have differed:
//!
//! - the largest gap between the two sources' monitoring statistics over every
//!   scored sample, and
//! - the smallest distance from any of those statistics to its control limit.
//!
//! When the first is many orders below the second, no thresholding of those
//! statistics could have produced a different alarm on any sample, and the
//! identical detection rates follow from arithmetic rather than from luck.
//! [`Trajectory`] does the same one level down, on the raw samples, in units of
//! each variable's own standard deviation.
//!
//! The gates are still the permutation tests. These two exist so that a run of
//! zeros can be read.

use tepsim_stats::{
    ControlLimits, CorrelationMatrix, Pca, Retention, Summary, alarms_above, augment_with_lags,
    detection_report,
};

use crate::Oracle;
use crate::tier5::battery::Calibrated;
use crate::tier5::{Run, Scenario, VARIABLES, run_fortran, run_port, seed, start};

/// Confidence level both control limits are drawn at.
///
/// 0.99, the level the TEP fault-detection literature quotes, giving a 1%
/// nominal false alarm rate. The *realised* rate differs from the nominal one
/// because plant data are not the independent normal draws the limits assume,
/// and the size of that gap is one of the numbers this tier reports.
pub const CONFIDENCE: f64 = 0.99;

/// How many principal components to keep.
///
/// The smallest number explaining 90% of the variance. See the module docs for
/// why a rule rather than a count.
pub const RETENTION: Retention = Retention::CumulativeVariance(0.90);

/// The lag counts the experiment fits a detector at.
///
/// Zero is static PCA. Two is dynamic PCA at the lag count Ku, Storer and
/// Georgakis and the TEP literature after them use.
pub const LAG_COUNTS: [usize; 2] = [0, 2];

/// How many consecutive alarms make a detection, for the delay only.
pub const PERSISTENCE: usize = 3;

/// How many distinct half-splits of the reference's test seeds the null uses.
///
/// Twenty, for [`crate::tier5::battery::CALIBRATION_SPLITS`]'s reason: the
/// smallest attainable permutation p-value is `1/21 = 0.048`, just under
/// [`ALPHA`], so the test can reject; fewer and it could not.
pub const CALIBRATION_SPLITS: usize = 20;

/// The significance level of the permutation test.
pub const ALPHA: f64 = 0.05;

/// Which simulator produced a set of runs.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Source {
    /// The vendored Fortran, driven through the FFI.
    Fortran,
    /// The Rust port, driven through [`tepsim::Simulation`].
    Port,
}

impl Source {
    /// Both, reference first.
    pub const ALL: [Self; 2] = [Self::Fortran, Self::Port];

    /// A short label for a report row.
    #[must_use]
    pub const fn label(self) -> &'static str {
        match self {
            Self::Fortran => "fortran",
            Self::Port => "rust",
        }
    }
}

/// How big an experiment to run, from `TEP_TIER6`.
///
/// `full` is the gate. Anything else, the variable being absent included, is
/// the smoke experiment, which exercises every code path and concludes
/// nothing: see [`Experiment::gated`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Battery {
    /// How many faults, `IDV(1)` through `IDV(faults)`.
    pub faults: usize,
    /// Fault-free runs per source used to fit the detectors.
    pub train_seeds: usize,
    /// Runs per source per scenario used to evaluate them.
    pub test_seeds: usize,
    /// Simulated hours per run.
    pub hours: usize,
}

impl Battery {
    /// The environment variable that selects the size.
    pub const ENV: &'static str = "TEP_TIER6";

    /// What `cargo test` runs: two faults, four hours, six test seeds.
    ///
    /// Six test seeds admit only `C(6,3)/2 = 10` distinct half-splits, so the
    /// permutation test cannot reject at 0.05 and every verdict comes back
    /// undecided. It covers the plumbing: both sources run, four models fit,
    /// every metric is finite, the two sources trip in the same places. It
    /// cannot support any statement about equivalence.
    pub const SMOKE: Self = Self {
        faults: 2,
        train_seeds: 6,
        test_seeds: 6,
        hours: 4,
    };

    /// What `TEP_TIER6=full` runs: every fault, 48 hours, twenty test seeds.
    ///
    /// Twenty test seeds give `C(20,10)/2 = 92,378` distinct half-splits, far
    /// more than the twenty the null draws, so the draws are a genuine sample
    /// of the null rather than the same few splits repeated.
    ///
    /// 430 runs per source, 860 in total, at 48 hours each.
    pub const FULL: Self = Self {
        faults: crate::tier5::FAULTS,
        train_seeds: 10,
        test_seeds: 20,
        hours: 48,
    };

    /// Which experiment to run.
    #[must_use]
    pub fn selected() -> Self {
        match std::env::var(Self::ENV).as_deref() {
            Ok("full") => Self::FULL,
            _ => Self::SMOKE,
        }
    }

    /// How many runs this experiment is, per source.
    #[must_use]
    pub const fn runs_per_source(&self) -> usize {
        self.train_seeds + self.test_seeds + self.faults * self.test_seeds
    }

    /// The seed indices the detectors are fitted to.
    pub fn training_seeds(&self) -> impl Iterator<Item = usize> + use<> {
        0..self.train_seeds
    }

    /// The seed indices the detectors are evaluated on, disjoint from the
    /// training ones.
    pub fn evaluation_seeds(&self) -> impl Iterator<Item = usize> + use<> {
        self.train_seeds..self.train_seeds + self.test_seeds
    }
}

/// A run's samples as one row-major matrix, `samples` by [`VARIABLES`].
fn flatten(run: &Run) -> Vec<f64> {
    let mut out = Vec::with_capacity(run.samples.len() * VARIABLES);
    for row in &run.samples {
        out.extend_from_slice(row);
    }
    out
}

/// Fit one detector to pooled fault-free runs.
///
/// Each run is lag-augmented **on its own** and the augmented rows are then
/// concatenated. Concatenating first and augmenting the result would
/// manufacture `lags` rows per boundary that pair the end of one run with the
/// start of the next, which are not observations of anything. There are only
/// twenty such rows in the full battery, so they would never show up as a
/// failure; they would show up as a slightly wrong model, forever.
///
/// # Panics
///
/// If `training` is empty, or if a run is shorter than `lags + 2` samples.
#[must_use]
pub fn fit(training: &[Run], lags: usize) -> Pca {
    let columns = VARIABLES * (lags + 1);
    let mut data = Vec::new();
    let mut rows = 0;
    for run in training {
        let flat = flatten(run);
        let (augmented, r, c) = augment_with_lags(&flat, run.samples.len(), VARIABLES, lags);
        assert_eq!(c, columns, "augmentation produced the wrong width");
        data.extend_from_slice(&augmented);
        rows += r;
    }
    Pca::fit(&data, rows, columns, RETENTION)
}

/// Both monitoring statistic series for every run in a set.
///
/// Each series is `samples - lags` long: the first `lags` samples of a run have
/// no complete history and cannot be scored by a dynamic model. That offset is
/// added back in [`Performance`] so that a delay is always in samples since the
/// fault, whatever the lag count.
///
/// Kept rather than thresholded straight away, because [`agreement`] needs the
/// statistics themselves and scoring every run twice to get them would double
/// the most expensive part of the evaluation.
#[derive(Clone, Debug)]
pub struct Scores {
    /// One Hotelling T-squared series per run.
    pub t_squared: Vec<Vec<f64>>,
    /// One squared-prediction-error series per run.
    pub spe: Vec<Vec<f64>>,
    /// How many samples of each run were unscorable, which is the lag count.
    pub offset: usize,
}

/// Score every run in a set.
#[must_use]
pub fn score(model: &Pca, runs: &[Run], lags: usize) -> Scores {
    let mut t_squared = Vec::with_capacity(runs.len());
    let mut spe = Vec::with_capacity(runs.len());
    for run in runs {
        let flat = flatten(run);
        let samples = run.samples.len();
        let scorable = samples.saturating_sub(lags);
        let mut t2 = Vec::with_capacity(scorable);
        let mut q = Vec::with_capacity(scorable);
        for t in lags..samples {
            // The window is `lags + 1` samples in time order, oldest first,
            // which is exactly what `augment_with_lags` turns into one
            // present-first row.
            let window = &flat[(t - lags) * VARIABLES..(t + 1) * VARIABLES];
            let (row, _, _) = augment_with_lags(window, lags + 1, VARIABLES, lags);
            let statistics = model.statistics(&row);
            t2.push(statistics.t_squared);
            q.push(statistics.spe);
        }
        t_squared.push(t2);
        spe.push(q);
    }
    Scores {
        t_squared,
        spe,
        offset: lags,
    }
}

/// How close the two sources' scored statistics came to each other, and how
/// close either came to its own control limit.
///
/// The pair of numbers that decides whether an alarm could have differed at
/// all. See the module docs.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Agreement {
    /// Largest `|T2_fortran - T2_port|` over every scored sample.
    pub t_squared_gap: f64,
    /// Largest `|SPE_fortran - SPE_port|`.
    pub spe_gap: f64,
    /// Smallest `|T2 - limit|` over every scored sample of either source.
    pub t_squared_approach: f64,
    /// Smallest `|SPE - limit|`.
    pub spe_approach: f64,
    /// How many samples went into both.
    pub samples: usize,
}

impl Agreement {
    /// The largest statistic gap divided by the closest approach to a limit.
    ///
    /// Below one, no alarm on any sample could have been decided differently by
    /// the two sources, because the two statistics never got closer to their
    /// limit than they are to each other. Above one is not a failure; it means
    /// at least one alarm *could* have flipped and the permutation tests are
    /// then doing real work rather than confirming an identity.
    #[must_use]
    pub fn decisive(&self) -> f64 {
        let t = self.t_squared_gap / self.t_squared_approach;
        let q = self.spe_gap / self.spe_approach;
        t.max(q)
    }
}

/// Compare two sources' scores sample by sample.
///
/// The two sets must be the same shape: the same runs, at the same seeds, in
/// the same order.
#[must_use]
pub fn agreement(a: &Scores, b: &Scores, limits: &ControlLimits) -> Agreement {
    let mut t_squared_gap = 0.0_f64;
    let mut spe_gap = 0.0_f64;
    let mut t_squared_approach = f64::INFINITY;
    let mut spe_approach = f64::INFINITY;
    let mut samples = 0;
    let pairs = [
        (
            &a.t_squared,
            &b.t_squared,
            limits.t_squared,
            &mut t_squared_gap,
            &mut t_squared_approach,
        ),
        (&a.spe, &b.spe, limits.spe, &mut spe_gap, &mut spe_approach),
    ];
    for (left, right, limit, gap, approach) in pairs {
        for (x, y) in left.iter().zip(right) {
            for (p, q) in x.iter().zip(y) {
                *gap = gap.max(libm::fabs(p - q));
                *approach = approach
                    .min(libm::fabs(p - limit))
                    .min(libm::fabs(q - limit));
            }
        }
    }
    for series in &a.t_squared {
        samples += series.len();
    }
    Agreement {
        t_squared_gap,
        spe_gap,
        t_squared_approach,
        spe_approach,
        samples,
    }
}

/// How far apart the two sources' raw samples are, in units of each variable's
/// own standard deviation.
///
/// Detector-independent, so it is measured once per scenario. Standard
/// deviations rather than absolute units because the 53 variables span kPa near
/// 2705 and mole fractions near 0.01, and one absolute number over the lot
/// would only ever report the pressure.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Trajectory {
    /// The largest absolute difference between paired samples.
    pub max_absolute: f64,
    /// The same difference divided by that variable's standard deviation over
    /// the reference's runs.
    pub max_scaled: f64,
    /// Which variable the scaled maximum was on, zero-based.
    pub variable: usize,
    /// Which sample of which run, zero-based.
    pub sample: usize,
    /// How many variables never moved in the reference and so were skipped.
    pub skipped: usize,
}

/// Compare two sources' raw samples, run for run.
///
/// # Panics
///
/// If the two sets are not the same shape.
#[must_use]
pub fn trajectory(fortran: &[Run], port: &[Run]) -> Trajectory {
    assert_eq!(fortran.len(), port.len(), "the two sources must be paired");
    // The reference's spread per variable, pooled over its runs, which is the
    // scale a difference is only meaningful against.
    let mut spread = [0.0_f64; VARIABLES];
    for (v, slot) in spread.iter_mut().enumerate() {
        let mut summary = Summary::new();
        for run in fortran {
            for row in &run.samples {
                summary.push(row[v]);
            }
        }
        *slot = summary.sd();
    }
    // Spelled out rather than `!(sd > 0.0)` so that the `NaN` case is visible
    // rather than riding on the negation.
    let skipped = spread
        .iter()
        .filter(|s| !s.is_finite() || **s <= 0.0)
        .count();

    let mut max_absolute = 0.0_f64;
    let mut max_scaled = 0.0_f64;
    let mut variable = 0;
    let mut sample = 0;
    for (run, (a, b)) in fortran.iter().zip(port).enumerate() {
        for (t, (x, y)) in a.samples.iter().zip(&b.samples).enumerate() {
            for v in 0..VARIABLES {
                let difference = libm::fabs(x[v] - y[v]);
                max_absolute = max_absolute.max(difference);
                if spread[v] > 0.0 {
                    let scaled = difference / spread[v];
                    if scaled > max_scaled {
                        max_scaled = scaled;
                        variable = v;
                        sample = run * a.samples.len() + t;
                    }
                }
            }
        }
    }
    Trajectory {
        max_absolute,
        max_scaled,
        variable,
        sample,
        skipped,
    }
}

/// One monitoring statistic's per-run performance over a set of runs.
#[derive(Clone, Debug)]
pub struct Performance {
    /// One alarm rate per run. The fault detection rate on faulted runs, the
    /// false alarm rate on fault-free ones; the arithmetic is the same and the
    /// reading differs.
    pub rate: Vec<f64>,
    /// One detection delay per run, in samples since the fault, censored at the
    /// record length when the persistence rule never fired.
    ///
    /// Censoring rather than dropping. A run that never detects is evidence
    /// about the detector, and omitting it would let a detector that fails on
    /// half its runs report a *better* mean delay than one that detects late on
    /// all of them.
    pub delay: Vec<f64>,
    /// How many runs the persistence rule fired on.
    pub detected: usize,
}

/// One source's runs for one scenario, under one detector.
#[derive(Clone, Debug)]
pub struct SourcePerformance {
    /// Hotelling's T-squared.
    pub t_squared: Performance,
    /// The squared prediction error.
    pub spe: Performance,
}

/// Threshold a set of scores and reduce them to per-run metrics.
///
/// `faulted` says whether the fault is present, which is what decides whether
/// the whole record is post-fault (detection rate, delay) or pre-fault (false
/// alarm rate). Both are read out of the same [`detection_report`] call, so the
/// definitions are the crate's and not restated here.
#[must_use]
pub fn evaluate(scores: &Scores, limits: &ControlLimits, faulted: bool) -> SourcePerformance {
    SourcePerformance {
        t_squared: reduce(&scores.t_squared, limits.t_squared, scores.offset, faulted),
        spe: reduce(&scores.spe, limits.spe, scores.offset, faulted),
    }
}

/// One statistic's series, thresholded and reduced to per-run metrics.
fn reduce(series: &[Vec<f64>], limit: f64, offset: usize, faulted: bool) -> Performance {
    let mut performance = Performance {
        rate: Vec::with_capacity(series.len()),
        delay: Vec::with_capacity(series.len()),
        detected: 0,
    };
    for statistic in series {
        let alarms = alarms_above(statistic, limit);
        // A faulted run is faulted from its first sample, so the whole record
        // is post-fault; a fault-free one has no post-fault part at all.
        let onset = if faulted { 0 } else { alarms.len() };
        let report = detection_report(&alarms, onset, PERSISTENCE);
        performance.rate.push(if faulted {
            report.fault_detection_rate
        } else {
            report.false_alarm_rate
        });
        // `offset` samples of the record were unscorable, so a delay measured
        // on the alarm series is that many samples late in plant time.
        match report.detection_delay {
            Some(delay) => {
                performance.detected += 1;
                performance.delay.push((delay + offset) as f64);
            }
            None => performance.delay.push((alarms.len() + offset) as f64),
        }
    }
    performance
}

/// The size-normalised distance between two groups of per-run values.
///
/// ```text
/// d(A, B) = |mean(A) - mean(B)| / sqrt(1/|A| + 1/|B|)
/// ```
///
/// The division is what makes a half-split of `n/2` against `n/2` comparable
/// with a cross-source comparison of `n` against `n`: under a common per-run
/// variance the expectation of `d^2` is that variance, whatever the group
/// sizes. Without it the null would be about `sqrt(2)` too wide and every
/// verdict would be biased toward passing.
///
/// `NaN` if either group is empty, which [`calibrate`] drops rather than
/// counting as a draw.
#[must_use]
pub fn separation(a: &[f64], b: &[f64]) -> f64 {
    if a.is_empty() || b.is_empty() {
        return f64::NAN;
    }
    let difference = Summary::of(a).mean() - Summary::of(b).mean();
    let scale = tepsim_stats::special::sqrt(1.0 / a.len() as f64 + 1.0 / b.len() as f64);
    libm::fabs(difference) / scale
}

/// SplitMix64's finaliser, for the split generator.
const fn mix(z: u64) -> u64 {
    let z = z.wrapping_mul(0x9E37_79B9_7F4A_7C15);
    let z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z ^ (z >> 27)
}

/// Up to `count` **distinct** balanced half-splits of `n` indices.
///
/// Deterministic: split `k` sorts the indices by a SplitMix hash of
/// `(index, k)` and takes the first and second halves. Two runs of the
/// experiment calibrate against the same splits, so a change in a verdict is a
/// change in the port.
///
/// Distinct is the part that matters, and it is why this is not
/// `tier5::battery`'s generator with a different name. With six test seeds
/// there are only ten balanced half-splits in existence; asking for twenty
/// would hand back the same ones repeatedly and the permutation p-value would
/// have a denominator of 21 with ten real draws behind it. Returning fewer is
/// the honest answer, and [`Calibrated::passes`] turns "fewer than nineteen"
/// into `None` rather than a verdict.
///
/// A split and its mirror image are the same split, so they are canonicalised
/// to the lexicographically smaller half first before the duplicate check. When
/// `n` is odd one index is left out of both halves, so that the two sides stay
/// the same size.
#[must_use]
pub fn half_splits(n: usize, count: usize) -> Vec<(Vec<usize>, Vec<usize>)> {
    let mut out: Vec<(Vec<usize>, Vec<usize>)> = Vec::new();
    if n < 2 || count == 0 {
        return out;
    }
    let half = n / 2;
    // A bounded search: the hashed order can repeat a split, and without a
    // ceiling a small `n` would spin forever looking for a split that does not
    // exist.
    let attempts = (count as u64).saturating_mul(64);
    for k in 0..attempts {
        if out.len() == count {
            break;
        }
        let mut order: Vec<usize> = (0..n).collect();
        order.sort_by_key(|i| mix(((*i as u64) << 32) | k));
        let mut left = order[..half].to_vec();
        let mut right = order[half..half * 2].to_vec();
        left.sort_unstable();
        right.sort_unstable();
        if right < left {
            core::mem::swap(&mut left, &mut right);
        }
        if !out.iter().any(|(l, r)| *l == left && *r == right) {
            out.push((left, right));
        }
    }
    out
}

/// Calibrate one cross-source statistic against the reference's own spread.
///
/// `reference` and `candidate` are per-run values at the same seeds, one from
/// each source. The null is built from `reference` alone, matching
/// [`crate::tier5::battery`]'s convention that margins come from the reference:
/// it is the stricter choice, because a port with a *wider* run-to-run spread
/// than the Fortran would be judged against the Fortran's narrower one.
#[must_use]
pub fn calibrate(
    name: &'static str,
    reference: &[f64],
    candidate: &[f64],
    splits: &[(Vec<usize>, Vec<usize>)],
) -> Calibrated {
    let take = |group: &[usize]| -> Vec<f64> { group.iter().map(|i| reference[*i]).collect() };
    let within = splits
        .iter()
        .map(|(left, right)| separation(&take(left), &take(right)))
        .filter(|value| value.is_finite())
        .collect();
    Calibrated {
        name,
        cross: separation(reference, candidate),
        within,
    }
}

/// The Pearson correlation between two sources' per-run values across seeds.
///
/// The pairing diagnostic described in the module docs. `NaN` when either side
/// never moves, which is the common and harmless case of a detection rate
/// pinned at one on every seed.
#[must_use]
pub fn pairing(a: &[f64], b: &[f64]) -> f64 {
    if a.len() != b.len() || a.len() < 2 {
        return f64::NAN;
    }
    CorrelationMatrix::of(&[a.to_vec(), b.to_vec()]).get(0, 1)
}

/// One metric, compared between the two evaluation sources for a fixed
/// detector.
#[derive(Clone, Debug)]
pub struct MetricComparison {
    /// What was measured.
    pub metric: &'static str,
    /// The per-run values on the Fortran's runs.
    pub fortran: Summary,
    /// The per-run values on the port's runs.
    pub port: Summary,
    /// `mean(port) - mean(fortran)`, in the metric's own units. Signed, and
    /// reported because it is the number a reader wants; the test runs on
    /// [`separation`] instead, for the reason in the module docs.
    pub difference: f64,
    /// [`pairing`], the seed-to-seed correlation between the two sources.
    pub pairing: f64,
    /// The largest within-source draw taken from the **port**'s seeds rather
    /// than the Fortran's. A diagnostic only: if it is much larger than
    /// `calibrated.within_max()`, the port's run-to-run spread is wider than
    /// the reference's and the verdict should be read with that in mind.
    pub port_within_max: f64,
    /// The cross-source value against the Fortran's within-source null.
    pub calibrated: Calibrated,
}

impl MetricComparison {
    /// Build one from two sources' per-run values.
    #[must_use]
    pub fn of(
        metric: &'static str,
        fortran: &[f64],
        port: &[f64],
        splits: &[(Vec<usize>, Vec<usize>)],
    ) -> Self {
        let port_null = calibrate(metric, port, fortran, splits);
        Self {
            metric,
            fortran: Summary::of(fortran),
            port: Summary::of(port),
            difference: Summary::of(port).mean() - Summary::of(fortran).mean(),
            pairing: pairing(fortran, port),
            port_within_max: port_null.within_max(),
            calibrated: calibrate(metric, fortran, port, splits),
        }
    }

    /// Whether this metric's gate passed, or `None` if it could not be gated.
    #[must_use]
    pub fn passes(&self) -> Option<bool> {
        self.calibrated.passes()
    }
}

/// Everything measured for one fault, for one trained detector.
#[derive(Clone, Debug)]
pub struct FaultReport {
    /// Which scenario.
    pub scenario: Scenario,
    /// How many of the Fortran's evaluation runs tripped.
    pub fortran_tripped: usize,
    /// How many of the port's did.
    pub port_tripped: usize,
    /// Detection rate and delay, for T-squared and for SPE.
    pub metrics: Vec<MetricComparison>,
    /// How close the two sources' statistics came to each other and to the
    /// limits, for this detector on this scenario.
    pub agreement: Agreement,
}

/// Everything measured for one detector: one training source, one lag count.
#[derive(Clone, Debug)]
pub struct DetectorReport {
    /// Which simulator's fault-free runs it was fitted to.
    pub trained_on: Source,
    /// The lag count. Zero is static PCA.
    pub lags: usize,
    /// How many components the retention rule kept.
    pub components: usize,
    /// What fraction of the variance they explain.
    pub explained_variance: f64,
    /// How many augmented rows it was fitted to.
    pub training_rows: usize,
    /// How many augmented columns, `53 * (lags + 1)`.
    pub columns: usize,
    /// Both control limits at [`CONFIDENCE`].
    pub limits: ControlLimits,
    /// The false alarm rate, on fault-free evaluation runs, for T-squared and
    /// for SPE.
    pub false_alarm: Vec<MetricComparison>,
    /// The same agreement diagnostic, on the fault-free evaluation runs.
    pub nominal_agreement: Agreement,
    /// One entry per fault.
    pub faults: Vec<FaultReport>,
}

impl DetectorReport {
    /// A label for a report row.
    #[must_use]
    pub fn label(&self) -> String {
        format!(
            "{} / {}",
            self.trained_on.label(),
            if self.lags == 0 {
                "PCA".to_string()
            } else {
                format!("DPCA({})", self.lags)
            }
        )
    }

    /// Every comparison this detector produced, false alarms first.
    pub fn comparisons(&self) -> impl Iterator<Item = (Scenario, &MetricComparison)> {
        self.false_alarm
            .iter()
            .map(|m| (Scenario::NOMINAL, m))
            .chain(
                self.faults
                    .iter()
                    .flat_map(|f| f.metrics.iter().map(|m| (f.scenario, m))),
            )
    }

    /// Every agreement diagnostic this detector produced, nominal first.
    pub fn agreements(&self) -> impl Iterator<Item = (Scenario, Agreement)> {
        core::iter::once((Scenario::NOMINAL, self.nominal_agreement))
            .chain(self.faults.iter().map(|f| (f.scenario, f.agreement)))
    }
}

/// The whole experiment.
#[derive(Clone, Debug)]
pub struct Experiment {
    /// The size it was run at.
    pub battery: Battery,
    /// Four detectors: two training sources by two lag counts.
    pub detectors: Vec<DetectorReport>,
    /// How many distinct half-splits the null actually had.
    pub splits: usize,
    /// The raw sample agreement, one entry per scenario, nominal first.
    ///
    /// Detector-independent, so it is measured once rather than four times.
    pub trajectories: Vec<(Scenario, Trajectory)>,
}

impl Experiment {
    /// Whether the permutation test can reject at [`ALPHA`] at this size.
    ///
    /// It can only when there are at least `1/ALPHA - 1` distinct half-splits,
    /// because the smallest attainable p-value is `1 / (splits + 1)`. This is
    /// the same arithmetic [`Calibrated::passes`] does per statistic, surfaced
    /// once so a report can say up front whether it is gating or reporting.
    #[must_use]
    pub fn gated(&self) -> bool {
        self.splits + 1 >= (1.0 / ALPHA) as usize
    }

    /// Every gate that this experiment could decide and that came out false.
    ///
    /// Returns `(detector label, scenario, metric)`. Empty at a size that
    /// cannot gate, which is why [`Experiment::gated`] has to be read beside
    /// it: an empty list from an ungated run is not a pass.
    #[must_use]
    pub fn failures(&self) -> Vec<(String, Scenario, &'static str)> {
        let mut out = Vec::new();
        for detector in &self.detectors {
            for (scenario, comparison) in detector.comparisons() {
                if comparison.passes() == Some(false) {
                    out.push((detector.label(), scenario, comparison.metric));
                }
            }
        }
        out
    }

    /// The smallest permutation p-value anywhere in the experiment, and where.
    ///
    /// The number to record: a run whose worst p-value moved from 0.62 to 0.10
    /// has changed even though both pass, and that is exactly the drift
    /// `CLAUDE.md` asks to be visible across iterations.
    #[must_use]
    pub fn worst(&self) -> Option<(String, Scenario, &'static str, f64)> {
        let mut worst: Option<(String, Scenario, &'static str, f64)> = None;
        for detector in &self.detectors {
            for (scenario, comparison) in detector.comparisons() {
                let p = comparison.calibrated.p_value();
                if !p.is_finite() {
                    continue;
                }
                if worst.as_ref().is_none_or(|(_, _, _, best)| p < *best) {
                    worst = Some((detector.label(), scenario, comparison.metric, p));
                }
            }
        }
        worst
    }

    /// How many comparisons had a cross-source statistic of exactly zero.
    ///
    /// Reported as a count because it is the difference between "the test
    /// passed" and "there was nothing to test". See the module docs.
    #[must_use]
    pub fn identical(&self) -> (usize, usize) {
        let mut zero = 0;
        let mut total = 0;
        for detector in &self.detectors {
            for (_, comparison) in detector.comparisons() {
                total += 1;
                if comparison.calibrated.cross == 0.0 {
                    zero += 1;
                }
            }
        }
        (zero, total)
    }

    /// The worst agreement anywhere: the largest statistic gap, the closest
    /// approach to a limit, and the largest ratio of the two.
    ///
    /// The third is the one that matters. Below one, no alarm decision anywhere
    /// in the experiment could have gone the other way.
    #[must_use]
    pub fn worst_agreement(&self) -> (f64, f64, f64) {
        let mut gap = 0.0_f64;
        let mut approach = f64::INFINITY;
        let mut decisive = 0.0_f64;
        for detector in &self.detectors {
            for (_, a) in detector.agreements() {
                gap = gap.max(a.t_squared_gap).max(a.spe_gap);
                approach = approach.min(a.t_squared_approach).min(a.spe_approach);
                let d = a.decisive();
                if d.is_finite() {
                    decisive = decisive.max(d);
                }
            }
        }
        (gap, approach, decisive)
    }

    /// The worst raw-sample disagreement anywhere, in units of the variable's
    /// own standard deviation, with the scenario it happened in.
    #[must_use]
    pub fn worst_trajectory(&self) -> Option<(Scenario, Trajectory)> {
        self.trajectories
            .iter()
            .copied()
            .max_by(|a, b| a.1.max_scaled.total_cmp(&b.1.max_scaled))
    }

    /// The largest seed-to-seed pairing correlation seen anywhere.
    ///
    /// Near zero means the two sources' runs at a shared seed have separated by
    /// the time the metric is taken, so treating them as unpaired is right. See
    /// the module docs.
    #[must_use]
    pub fn worst_pairing(&self) -> f64 {
        self.detectors
            .iter()
            .flat_map(DetectorReport::comparisons)
            .map(|(_, c)| libm::fabs(c.pairing))
            .filter(|v| v.is_finite())
            .fold(0.0_f64, f64::max)
    }
}

/// The names of the four metrics compared per fault, in order.
const FAULT_METRICS: [&str; 4] = [
    "T2 detection rate",
    "T2 detection delay",
    "SPE detection rate",
    "SPE detection delay",
];

/// The names of the two compared on the fault-free runs.
const NOMINAL_METRICS: [&str; 2] = ["T2 false alarm rate", "SPE false alarm rate"];

/// Compare one scenario across the two evaluation sources, for one detector.
fn compare(
    fortran: &SourcePerformance,
    port: &SourcePerformance,
    faulted: bool,
    splits: &[(Vec<usize>, Vec<usize>)],
) -> Vec<MetricComparison> {
    if faulted {
        let sources = [
            (
                FAULT_METRICS[0],
                &fortran.t_squared.rate,
                &port.t_squared.rate,
            ),
            (
                FAULT_METRICS[1],
                &fortran.t_squared.delay,
                &port.t_squared.delay,
            ),
            (FAULT_METRICS[2], &fortran.spe.rate, &port.spe.rate),
            (FAULT_METRICS[3], &fortran.spe.delay, &port.spe.delay),
        ];
        sources
            .into_iter()
            .map(|(name, a, b)| MetricComparison::of(name, a, b, splits))
            .collect()
    } else {
        let sources = [
            (
                NOMINAL_METRICS[0],
                &fortran.t_squared.rate,
                &port.t_squared.rate,
            ),
            (NOMINAL_METRICS[1], &fortran.spe.rate, &port.spe.rate),
        ];
        sources
            .into_iter()
            .map(|(name, a, b)| MetricComparison::of(name, a, b, splits))
            .collect()
    }
}

/// Run the whole cross-source experiment.
///
/// Runs are generated one scenario at a time and dropped once every detector
/// has scored them, so the peak memory is two sources' worth of one scenario
/// rather than the whole battery: 430 runs per source at 48 hours would
/// otherwise be about 350 MB of samples.
///
/// # Panics
///
/// If the battery asks for fewer than two training or test seeds, or for a
/// horizon too short to lag-augment.
#[must_use]
pub fn run(oracle: &mut Oracle, battery: Battery) -> Experiment {
    assert!(
        battery.train_seeds >= 2 && battery.test_seeds >= 2,
        "a source needs at least two runs to have a run-to-run spread"
    );
    let start = start(oracle);

    // Fault-free training runs. Generated once, used to fit every detector,
    // then dropped: they are the largest thing held at any moment.
    let mut training = Vec::new();
    for source in Source::ALL {
        let runs: Vec<Run> = battery
            .training_seeds()
            .map(|s| match source {
                Source::Fortran => {
                    run_fortran(oracle, &start, Scenario::NOMINAL, seed(s), battery.hours)
                }
                Source::Port => run_port(&start, Scenario::NOMINAL, seed(s), battery.hours),
            })
            .collect();
        training.push((source, runs));
    }

    let mut detectors = Vec::new();
    let mut models = Vec::new();
    for lags in LAG_COUNTS {
        for (source, runs) in &training {
            let model = fit(runs, lags);
            let limits = model.limits(CONFIDENCE);
            detectors.push(DetectorReport {
                trained_on: *source,
                lags,
                components: model.retained(),
                explained_variance: model.explained_variance(),
                training_rows: model.samples(),
                columns: model.variables(),
                limits,
                false_alarm: Vec::new(),
                nominal_agreement: Agreement {
                    t_squared_gap: f64::NAN,
                    spe_gap: f64::NAN,
                    t_squared_approach: f64::NAN,
                    spe_approach: f64::NAN,
                    samples: 0,
                },
                faults: Vec::new(),
            });
            models.push((model, limits, lags));
        }
    }
    drop(training);

    let splits = half_splits(battery.test_seeds, CALIBRATION_SPLITS);

    // The evaluation runs for one scenario, at seeds the models never saw.
    let evaluation = |oracle: &mut Oracle, scenario: Scenario| -> [Vec<Run>; 2] {
        let fortran = battery
            .evaluation_seeds()
            .map(|s| run_fortran(oracle, &start, scenario, seed(s), battery.hours))
            .collect();
        let port = battery
            .evaluation_seeds()
            .map(|s| run_port(&start, scenario, seed(s), battery.hours))
            .collect();
        [fortran, port]
    };
    let tripped = |runs: &[Run]| runs.iter().filter(|r| r.tripped.is_some()).count();

    let mut trajectories = Vec::with_capacity(battery.faults + 1);

    let [fortran, port] = evaluation(oracle, Scenario::NOMINAL);
    trajectories.push((Scenario::NOMINAL, trajectory(&fortran, &port)));
    for (report, (model, limits, lags)) in detectors.iter_mut().zip(&models) {
        let a = score(model, &fortran, *lags);
        let b = score(model, &port, *lags);
        report.nominal_agreement = agreement(&a, &b, limits);
        report.false_alarm = compare(
            &evaluate(&a, limits, false),
            &evaluate(&b, limits, false),
            false,
            &splits,
        );
    }
    drop((fortran, port));

    for fault in 1..=battery.faults {
        let scenario = Scenario::fault(fault);
        let [fortran, port] = evaluation(oracle, scenario);
        trajectories.push((scenario, trajectory(&fortran, &port)));
        let (fortran_tripped, port_tripped) = (tripped(&fortran), tripped(&port));
        for (report, (model, limits, lags)) in detectors.iter_mut().zip(&models) {
            let a = score(model, &fortran, *lags);
            let b = score(model, &port, *lags);
            report.faults.push(FaultReport {
                scenario,
                fortran_tripped,
                port_tripped,
                metrics: compare(
                    &evaluate(&a, limits, true),
                    &evaluate(&b, limits, true),
                    true,
                    &splits,
                ),
                agreement: agreement(&a, &b, limits),
            });
        }
    }

    Experiment {
        battery,
        detectors,
        splits: splits.len(),
        trajectories,
    }
}
