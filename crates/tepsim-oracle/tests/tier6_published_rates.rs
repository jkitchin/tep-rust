//! Tier 6, second half: do the detectors reproduce the fault detection rates
//! the TEP literature reports?
//!
//! B-0049a. Tier 6 proper asks whether a detector can tell the Fortran from
//! the port. It answers that without ever asking whether the detector is any
//! good, and it does not have to: a broken detector that is broken identically
//! on both sources still passes a cross-source test. This file asks the other
//! question. Train the Tier 6 PCA detector on the published `d00`, evaluate it
//! on the published `d01_te` through `d21_te`, and compare the detection rates
//! against what thirty years of TEP fault-detection papers say they should be.
//!
//! # What is asserted, and why it is not a table of published numbers
//!
//! **No published number is quoted here as a target.** The papers that carry
//! the per-fault tables were not available while this was written, and
//! transcribing remembered numbers into a test would produce an assertion that
//! looks like a citation and is not one. That is worse than no test: a future
//! reader cannot tell a checked number from an invented one, and the whole
//! point of this project's validation ladder is that they can.
//!
//! What *is* asserted is the part of the literature that is qualitative,
//! unanimous, and reproducible from first principles, namely which faults a
//! PCA detector can see at all:
//!
//! - **Faults 3, 9 and 15 are near-undetectable.** They are step and
//!   random-variation changes in the D feed temperature (3), the C feed
//!   temperature (9) and the condenser cooling water valve (15) whose effect
//!   on the 52 measured variables is inside the ordinary operating spread. The
//!   TEP literature reports missed detection rates near 90% and above for all
//!   three, which is to say detection rates near the false alarm rate. This is
//!   the most-repeated fact about the benchmark; every survey of it says so,
//!   and papers routinely exclude the three from their averages for that
//!   reason.
//! - **Faults 1, 2, 6, 7, 8, 12, 13, 14 and 18 are detected essentially
//!   always.** Step changes in feed composition, feed ratio, feed loss and
//!   header pressure move the plant far enough that both statistics leave
//!   their limits and stay outside them.
//!
//! Those two facts are the assertions. They are checked as a **separation**
//! rather than against a threshold: the worst easy fault has to score above
//! the best hard fault by a wide margin, and the hard faults have to sit near
//! the detector's own false alarm rate. A separation claim cannot be met by a
//! detector that alarms on everything or on nothing, which a pair of absolute
//! thresholds could be.
//!
//! If a fault in the easy list came out hard, that would be a finding about
//! this implementation. If 3, 9 or 15 came out easy, that would be a finding
//! too, and a more alarming one: it would mean the detector is firing on
//! something other than the fault.
//!
//! # The detector is Tier 6's, unchanged
//!
//! [`tepsim_oracle::tier6`] fixes every parameter with its reason, before any
//! number was on the screen: PCA on the correlation matrix, the smallest
//! number of components explaining 90% of the variance, 99% control limits,
//! Hotelling's T-squared and SPE. This file re-uses them rather than choosing
//! its own, so the rates below are the rates of the detector the cross-source
//! experiment used, and `the_default_detector_is_tier_sixs` asserts that the
//! local fit really is the same fit.
//!
//! # What the published files bring that Tier 6's own runs do not
//!
//! Tier 6 switches a fault on at `t = 0`, so a run has no fault-free prologue
//! and the false alarm rate has to come from separate nominal runs. The
//! published `_te` files have 160 fault-free samples before the fault, which
//! is what makes them comparable with the literature: the same record supplies
//! the detection rate and a false alarm rate under the same conditions.
//!
//! They also carry the `IDV(12)` question. Tier 6's runs force `IDV(12)` on at
//! eight hours because `temain_mod.f:366-368` does; Tier 7 measured the
//! published files and concluded they were generated with that line
//! **replaced** rather than kept. Nothing here depends on which is right,
//! because nothing here regenerates anything: the published bytes are the
//! input.
//!
//! # Size
//!
//! `TEP_TIER6_RATES` selects it, the way `TEP_TIER5`, `TEP_TIER6` and
//! `TEP_TIER7` do. The default runs the one fixed detector over all 22 files
//! and takes well under a second, because it reads 44,000 rows off disk and
//! diagonalises one 53-by-53 matrix. `full` adds dynamic PCA at two lags and a
//! sweep over retention rules, which is the sensitivity analysis: if the
//! easy/hard separation survives every retention rule and both lag counts, it
//! is a property of the data and not of a parameter choice.

#![cfg(feature = "oracle")]
#![allow(
    clippy::needless_range_loop,
    reason = "matrix code indexes row-major storage by (i, j); enumerate() obscures it"
)]

use tepsim_oracle::tier5::{Run, VARIABLES};
use tepsim_oracle::tier6::{CONFIDENCE, PERSISTENCE, RETENTION, score};
use tepsim_oracle::tier7::{Published, Split, TESTING_ROWS};
use tepsim_stats::{Pca, Retention, alarms_above, augment_with_lags, detection_report};

/// The sample the fault arrives at in every `dNN_te.dat`.
///
/// `temain_mod.f:226` sets `SSPTS = 3600 * 8` and the driver samples every 180
/// steps, so eight hours is row 160. Tier 7 confirms it against the bytes:
/// every `_te` file sits at the nominal operating point for its first 160 rows
/// and departs in row 160.
const ONSET: usize = 160;

/// The faults the TEP literature agrees a PCA detector essentially always
/// catches.
///
/// Step changes in the A/C feed ratio (1), the B composition (2), the reactor
/// cooling water inlet temperature (6, a feed loss), the condenser cooling
/// water inlet temperature (7, a header pressure loss), the A/B/C feed
/// composition (8), and the unknown faults 12, 13, 14 and 18.
const EASY: [usize; 9] = [1, 2, 6, 7, 8, 12, 13, 14, 18];

/// The faults the TEP literature agrees a PCA detector cannot see.
///
/// The D feed temperature step (3), the C feed temperature random variation
/// (9), and the condenser cooling water valve (15). Their effect on the 52
/// measured variables is inside the ordinary operating spread, so a detector
/// scoring them well above its own false alarm rate is detecting something
/// that is not the fault.
const HARD: [usize; 3] = [3, 9, 15];

/// How big an analysis to run, from `TEP_TIER6_RATES`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct Size {
    /// Whether to add DPCA at two lags and the retention sweep.
    sweep: bool,
}

impl Size {
    const ENV: &'static str = "TEP_TIER6_RATES";

    fn selected() -> Self {
        Self {
            sweep: std::env::var(Self::ENV).as_deref() == Ok("full"),
        }
    }
}

/// One detector: a lag count and a retention rule.
#[derive(Clone, Copy, Debug)]
struct Detector {
    lags: usize,
    retention: Retention,
    label: &'static str,
}

/// Fit a PCA model to a set of runs, at a chosen lag count and retention rule.
///
/// The body of [`tepsim_oracle::tier6::fit`] with the retention rule made a
/// parameter, so the sweep below can vary it. Each run is lag-augmented on its
/// own before the rows are concatenated, for the reason that function gives:
/// augmenting a concatenation would manufacture rows pairing the end of one
/// run with the start of the next.
///
/// `the_default_detector_is_tier_sixs` asserts that this reduces to
/// `tier6::fit` at `tier6::RETENTION`, so the numbers reported here belong to
/// the detector Tier 6 fixed rather than to a variant of it.
fn fit_with(training: &[Run], lags: usize, retention: Retention) -> Pca {
    let columns = VARIABLES * (lags + 1);
    let mut data = Vec::new();
    let mut rows = 0;
    for run in training {
        let mut flat = Vec::with_capacity(run.samples.len() * VARIABLES);
        for row in &run.samples {
            flat.extend_from_slice(row);
        }
        let (augmented, r, c) = augment_with_lags(&flat, run.samples.len(), VARIABLES, lags);
        assert_eq!(c, columns, "augmentation produced the wrong width");
        data.extend_from_slice(&augmented);
        rows += r;
    }
    Pca::fit(&data, rows, columns, retention)
}

/// One fault's numbers under one detector.
#[derive(Clone, Copy, Debug)]
struct Rates {
    fault: usize,
    /// Detection rate over the post-fault samples, T-squared.
    t_squared: f64,
    /// The same, SPE.
    spe: f64,
    /// False alarms among the file's own 160 fault-free samples, T-squared.
    t_squared_false: f64,
    /// The same, SPE.
    spe_false: f64,
    /// Samples after the onset before three consecutive alarms, T-squared.
    t_squared_delay: Option<usize>,
    /// The same, SPE.
    spe_delay: Option<usize>,
    /// Where the onset landed in the statistic series, which is `ONSET` less
    /// the lag count.
    ///
    /// Carried so a test can assert it corresponds to plant sample [`ONSET`].
    /// The lag shift is the one piece of index arithmetic in this file that a
    /// detection rate cannot check: moving a boundary two samples out of 960
    /// changes a rate over 800 samples in the fourth decimal, and a mutation
    /// run confirmed that dropping the shift entirely fails nothing.
    onset_index: usize,
    /// How long the statistic series was.
    series: usize,
}

/// Score one published testing file against a model and reduce it to rates.
///
/// The lag offset is what makes this fiddly and it is handled once, here. A
/// model at `lags` lags cannot score the first `lags` samples of a record, so
/// the statistic series is that much shorter and its index `i` is sample
/// `i + lags`. The fault onset therefore sits at index `ONSET - lags` in the
/// series, and a delay measured on the series is already in samples since the
/// fault because both ends shifted by the same amount.
fn rates_for(
    model: &Pca,
    limits: &tepsim_stats::ControlLimits,
    file: &Published,
    lags: usize,
) -> Rates {
    let run = file.run();
    assert_eq!(
        run.samples.len(),
        TESTING_ROWS,
        "{} is not a 960-row testing file",
        file.name()
    );
    let scores = score(model, core::slice::from_ref(&run), lags);
    let onset = ONSET - lags;

    let t2 = alarms_above(&scores.t_squared[0], limits.t_squared);
    let q = alarms_above(&scores.spe[0], limits.spe);
    let t2_report = detection_report(&t2, onset, PERSISTENCE);
    let q_report = detection_report(&q, onset, PERSISTENCE);

    Rates {
        fault: file.fault,
        t_squared: t2_report.fault_detection_rate,
        spe: q_report.fault_detection_rate,
        t_squared_false: t2_report.false_alarm_rate,
        spe_false: q_report.false_alarm_rate,
        t_squared_delay: t2_report.detection_delay,
        spe_delay: q_report.detection_delay,
        onset_index: onset,
        series: scores.t_squared[0].len(),
    }
}

/// The alarm rate on a fault-free record, for each statistic.
///
/// `onset` is the record length, so every sample counts as pre-fault and the
/// answer is a pure false alarm rate. This is the number the hard faults'
/// detection rates are read against: a fault whose detection rate equals this
/// is invisible to the detector, whatever the absolute value happens to be.
fn false_alarms_on(
    model: &Pca,
    limits: &tepsim_stats::ControlLimits,
    run: &Run,
    lags: usize,
) -> (f64, f64) {
    let scores = score(model, core::slice::from_ref(run), lags);
    let length = scores.t_squared[0].len();
    let t2 = alarms_above(&scores.t_squared[0], limits.t_squared);
    let q = alarms_above(&scores.spe[0], limits.spe);
    (
        detection_report(&t2, length, PERSISTENCE).false_alarm_rate,
        detection_report(&q, length, PERSISTENCE).false_alarm_rate,
    )
}

/// Everything one detector produced.
struct Evaluation {
    model: Pca,
    limits: tepsim_stats::ControlLimits,
    /// `(T-squared, SPE)` alarm rates on `d00_te`, which has no fault.
    nominal: (f64, f64),
    /// The same on `d00.dat`, the file the model was fitted to.
    ///
    /// The phase-I rate. Read beside `nominal` it separates "the control limit
    /// does not describe this data" from "the test file is not like the
    /// training file": if both are far above `1 - CONFIDENCE` the limit's
    /// distributional assumption is what fails, and if only the second is, the
    /// two published files differ.
    training: (f64, f64),
    rates: Vec<Rates>,
}

/// Fit one detector to the published `d00` training file and evaluate it on
/// every published testing file.
fn evaluate(detector: Detector) -> Evaluation {
    let training_run = Published {
        fault: 0,
        split: Split::Training,
    }
    .run();
    let model = fit_with(
        core::slice::from_ref(&training_run),
        detector.lags,
        detector.retention,
    );
    let limits = model.limits(CONFIDENCE);
    let nominal_run = Published {
        fault: 0,
        split: Split::Testing,
    }
    .run();
    Evaluation {
        nominal: false_alarms_on(&model, &limits, &nominal_run, detector.lags),
        training: false_alarms_on(&model, &limits, &training_run, detector.lags),
        rates: (1..=21)
            .map(|fault| {
                rates_for(
                    &model,
                    &limits,
                    &Published {
                        fault,
                        split: Split::Testing,
                    },
                    detector.lags,
                )
            })
            .collect(),
        model,
        limits,
    }
}

/// The class a fault is in, for the report.
fn class(fault: usize) -> &'static str {
    if EASY.contains(&fault) {
        "easy"
    } else if HARD.contains(&fault) {
        "HARD"
    } else {
        ""
    }
}

fn print_table(detector: Detector, evaluation: &Evaluation) {
    let Evaluation {
        model,
        limits,
        nominal,
        training,
        rates,
    } = evaluation;
    println!(
        "\n-- {} --\n  {} components explaining {:.4} of the variance, from {} \
         rows of {} columns; limits T2 {:.3}, SPE {:.3} at {CONFIDENCE} confidence\n  \
         alarm rate on d00_te   (960 fault-free samples): T2 {:.4}, SPE {:.4}\n  \
         alarm rate on d00.dat  (the training file):      T2 {:.4}, SPE {:.4}\n  \
         nominal, from the confidence level:              {:.4}",
        detector.label,
        model.retained(),
        model.explained_variance(),
        model.samples(),
        model.variables(),
        limits.t_squared,
        limits.spe,
        nominal.0,
        nominal.1,
        training.0,
        training.1,
        1.0 - CONFIDENCE
    );
    println!(
        "  {:<6} {:<6} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9}",
        "fault", "class", "T2 rate", "SPE rate", "T2 FAR", "SPE FAR", "T2 delay", "SPE delay"
    );
    let delay = |d: Option<usize>| match d {
        Some(d) => format!("{d}"),
        None => "none".to_string(),
    };
    for r in rates {
        println!(
            "  {:<6} {:<6} {:>9.4} {:>9.4} {:>9.4} {:>9.4} {:>9} {:>9}",
            r.fault,
            class(r.fault),
            r.t_squared,
            r.spe,
            r.t_squared_false,
            r.spe_false,
            delay(r.t_squared_delay),
            delay(r.spe_delay)
        );
    }
}

/// The largest fault-free alarm rate a statistic can have and still be called
/// a detector.
///
/// A statistic that alarms on a quarter of fault-free samples carries almost
/// no information about whether a fault is present, and a detection rate read
/// off it is not a detection rate. The threshold is deliberately far above any
/// nominal level and far below the one case that trips it: dynamic PCA's SPE
/// on these files alarms on **52.7%** of `d00_te`, which is a coin flip. See
/// [`Statistic::usable`].
const UNUSABLE_ABOVE: f64 = 0.25;

/// Which monitoring statistic a row of the report is about.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Statistic {
    TSquared,
    Spe,
}

impl Statistic {
    const BOTH: [Self; 2] = [Self::TSquared, Self::Spe];

    const fn label(self) -> &'static str {
        match self {
            Self::TSquared => "T2",
            Self::Spe => "SPE",
        }
    }

    const fn rate(self, r: &Rates) -> f64 {
        match self {
            Self::TSquared => r.t_squared,
            Self::Spe => r.spe,
        }
    }

    /// This statistic's alarm rate on the fault-free `d00_te`.
    const fn false_alarms(self, nominal: (f64, f64)) -> f64 {
        match self {
            Self::TSquared => nominal.0,
            Self::Spe => nominal.1,
        }
    }

    /// Whether the statistic is discriminating at all on this data.
    fn usable(self, nominal: (f64, f64)) -> bool {
        self.false_alarms(nominal) <= UNUSABLE_ABOVE
    }
}

/// How one statistic separates the two fault classes.
#[derive(Clone, Copy, Debug)]
struct Separation {
    statistic: Statistic,
    /// The lowest detection rate among [`EASY`].
    worst_easy: f64,
    /// The highest among [`HARD`].
    best_hard: f64,
    /// This statistic's own false alarm rate on `d00_te`.
    false_alarms: f64,
    usable: bool,
}

impl Separation {
    fn gap(&self) -> f64 {
        self.worst_easy - self.best_hard
    }

    /// Whether this statistic on its own reproduces the literature's split:
    /// it discriminates, it catches every easy fault, and the classes do not
    /// overlap.
    fn carries_the_claim(&self) -> bool {
        self.usable && self.worst_easy >= EASY_FLOOR && self.gap() > MINIMUM_GAP
    }
}

/// The separation each statistic achieves.
///
/// Per statistic, not on the better of the two, and that is the correction the
/// dynamic-PCA row forced. A `max(T2, SPE)` summary is dominated by whichever
/// statistic is *louder*, which is the wrong one when one of them is alarming
/// on half the fault-free data: it reported fault 3 at 0.62 for DPCA(2) and
/// called the classes overlapping, when what had actually happened was that
/// SPE had stopped being a detector and T-squared was separating the classes
/// by 0.86.
fn separations(rates: &[Rates], nominal: (f64, f64)) -> Vec<Separation> {
    Statistic::BOTH
        .into_iter()
        .map(|statistic| Separation {
            statistic,
            worst_easy: rates
                .iter()
                .filter(|r| EASY.contains(&r.fault))
                .map(|r| statistic.rate(r))
                .fold(f64::INFINITY, f64::min),
            best_hard: rates
                .iter()
                .filter(|r| HARD.contains(&r.fault))
                .map(|r| statistic.rate(r))
                .fold(0.0_f64, f64::max),
            false_alarms: statistic.false_alarms(nominal),
            usable: statistic.usable(nominal),
        })
        .collect()
}

/// How many times its own false alarm rate a hard fault is allowed to reach.
///
/// The literature's claim about faults 3, 9 and 15 is that they are
/// *indistinguishable from fault-free operation*, so the claim has to be made
/// against the detector's own alarm rate on fault-free data and not against an
/// absolute number.
///
/// **That is not the form this test was first written in, and the difference
/// is a finding.** The first version asserted "a hard fault scores at most
/// 0.20", which is what an absolute reading of "near-undetectable" suggests.
/// It failed, on fault 15's SPE at 0.2112. Fault 15 is not detectable; the
/// threshold was wrong, and the table this test prints says why.
///
/// The SPE limit alarms on 15.0% of `d00_te`, against the 1% its confidence
/// level nominally buys. The obvious reading is that the Jackson-Mudholkar
/// limit does not describe autocorrelated plant data, which is the caveat
/// [`tepsim_stats::detection::false_alarm_rate`] exists to raise. The numbers
/// say otherwise: on `d00.dat`, the file the model was **fitted** to, the same
/// limit alarms on 0.2% of samples and T-squared on none at all, both *under*
/// nominal. So the limit describes the training data perfectly well, and the
/// 15% is not a calibration failure. It is `d00_te` differing from `d00`.
///
/// That is consistent with what Tier 7 found independently and could not
/// resolve: [`tepsim_oracle::tier7::Unknown::TrainingProtocolIsUndocumented`].
/// The 22 training files are 480 or 500 rows against the driver's 960 and
/// nothing in `temain_mod.f` produces them, so they were made by editing
/// constants nobody recorded. A residual structure that shifts between the two
/// halves of the published pair is exactly what an unrecorded protocol change
/// would look like, and it is why an absolute alarm-rate threshold cannot be
/// the form of this claim.
///
/// Three, then, with the largest observed ratio well under it. The value is a
/// round number and not a fitted one; what makes the test meaningful is the
/// gap between this and [`EASY_MULTIPLE`], which is a factor of sixteen in the
/// measured data.
const HARD_MULTIPLE: f64 = 3.0;

/// How many times its own false alarm rate an easy fault must reach, on at
/// least one statistic.
const EASY_MULTIPLE: f64 = 5.0;

/// The detection rate every easy fault must reach on the statistic that
/// carries the claim.
///
/// The literature says "essentially always" for these nine. Four fifths is a
/// weak reading of that and is the point: the test is not trying to pin a
/// published number, it is trying to fail loudly if a fault the literature
/// calls trivial comes out hard.
const EASY_FLOOR: f64 = 0.80;

/// How far apart the two classes must be on the statistic that carries the
/// claim.
const MINIMUM_GAP: f64 = 0.5;

/// The two facts the TEP literature agrees on, asserted.
///
/// Applied to every detector in the sweep, not only the headline one, so the
/// claim is about the benchmark rather than about a parameter choice.
fn check(detector: Detector, evaluation: &Evaluation) {
    let Evaluation { nominal, rates, .. } = evaluation;
    let who = detector.label;

    // Every rate is a fraction of a real denominator, and the window the
    // denominators come from is the one this file claims: a series `lags`
    // shorter than the file, split at the score index that corresponds to
    // plant sample `ONSET`.
    for r in rates {
        for (name, value) in [
            ("T2 rate", r.t_squared),
            ("SPE rate", r.spe),
            ("T2 FAR", r.t_squared_false),
            ("SPE FAR", r.spe_false),
        ] {
            assert!(
                (0.0..=1.0).contains(&value),
                "{who}: fault {}: {name} is {value}, which is not a rate",
                r.fault
            );
        }
        assert_eq!(
            r.series,
            TESTING_ROWS - detector.lags,
            "{who}: fault {}: a {}-lag model scored {} of {TESTING_ROWS} samples",
            r.fault,
            detector.lags,
            r.series
        );
        assert_eq!(
            r.onset_index + detector.lags,
            ONSET,
            "{who}: fault {}: score index {} is plant sample {}, not {ONSET}",
            r.fault,
            r.onset_index,
            r.onset_index + detector.lags
        );
    }

    // Fact one: every easy fault alarms far more often than fault-free data
    // does. Stated as a ratio, so it holds for a statistic whose absolute
    // false alarm rate is high as well as for one whose is low.
    for r in rates.iter().filter(|r| EASY.contains(&r.fault)) {
        let ratio = (r.t_squared / nominal.0).max(r.spe / nominal.1);
        assert!(
            ratio >= EASY_MULTIPLE,
            "{who}: fault {} is one the TEP literature reports as detected \
             essentially always, and it alarms only {ratio:.2} times as often \
             as fault-free data (T2 {:.4} against {:.4}, SPE {:.4} against \
             {:.4}). That is a finding about the implementation, not about the \
             benchmark.",
            r.fault,
            r.t_squared,
            nominal.0,
            r.spe,
            nominal.1
        );
    }

    // Fact two: the famous three sit within a small multiple of the
    // detector's own false alarm rate, on **both** statistics. The check is
    // weak for a statistic whose false alarm rate is already large, which is
    // why `usable` is reported beside it rather than folded into it.
    for r in rates.iter().filter(|r| HARD.contains(&r.fault)) {
        for statistic in Statistic::BOTH {
            let (rate, far) = (statistic.rate(r), statistic.false_alarms(*nominal));
            assert!(
                rate <= HARD_MULTIPLE * far,
                "{who}: fault {} is one of the three the TEP literature reports \
                 as near-undetectable by PCA, and {} caught {rate:.4} of its \
                 post-fault samples against a false alarm rate of {far:.4} on \
                 fault-free data, a ratio of {:.2}. A detector that finds fault \
                 {} easy is firing on something that is not the fault.",
                r.fault,
                statistic.label(),
                rate / far,
                r.fault
            );
        }
    }

    // Fact three, the one a pair of thresholds could not express: some
    // monitoring statistic in this scheme separates the two classes cleanly.
    let separations = separations(rates, *nominal);
    println!(
        "  {:<5} {:>10} {:>11} {:>8} {:>10} {:>8}",
        "stat", "worst easy", "best hard", "gap", "d00_te FAR", "usable"
    );
    for s in &separations {
        println!(
            "  {:<5} {:>10.4} {:>11.4} {:>8.4} {:>10.4} {:>8}",
            s.statistic.label(),
            s.worst_easy,
            s.best_hard,
            s.gap(),
            s.false_alarms,
            if s.usable { "yes" } else { "NO" }
        );
    }
    assert!(
        separations.iter().any(Separation::carries_the_claim),
        "{who}: no monitoring statistic reproduces the literature's split. \
         A statistic carries it when it alarms on at most {UNUSABLE_ABOVE} of \
         fault-free samples, catches every easy fault at {EASY_FLOOR} or \
         better, and leaves a gap above {MINIMUM_GAP} to the best of faults \
         {HARD:?}. Measured: {separations:?}"
    );

    // Fact four: the fault really is at sample `ONSET`. On a statistic that
    // carries the claim, an easy fault's own 160 pre-fault samples have to be
    // quiet compared with its post-fault samples. Without this the onset could
    // be anywhere: moving it to 320 sweeps 160 detected samples into the
    // pre-fault window and out of the post-fault one, which changes a rate over
    // 800 samples in the second decimal and nothing else notices. A mutation
    // run confirmed that, and this is the assertion that answers it.
    for s in separations.iter().filter(|s| s.carries_the_claim()) {
        for r in rates.iter().filter(|r| EASY.contains(&r.fault)) {
            let (post, pre) = match s.statistic {
                Statistic::TSquared => (r.t_squared, r.t_squared_false),
                Statistic::Spe => (r.spe, r.spe_false),
            };
            assert!(
                pre < 0.25 * post,
                "{who}: {} alarms on {pre:.4} of fault {}'s 160 pre-fault \
                 samples against {post:.4} of its post-fault samples. Either \
                 the fault is not at sample {ONSET} or the detector is alarming \
                 on fault-free data.",
                s.statistic.label(),
                r.fault
            );
        }
    }
}

/// The local fit is Tier 6's fit, so the rates below belong to Tier 6's
/// detector.
///
/// `fit_with` exists only so the sensitivity sweep can vary the retention
/// rule. If it drifted from `tier6::fit` the headline table would be of some
/// other detector while claiming to be of that one, and nothing else in this
/// file would notice.
#[test]
fn the_default_detector_is_tier_sixs() {
    let training = Published {
        fault: 0,
        split: Split::Training,
    }
    .run();
    for lags in [0_usize, 2] {
        let local = fit_with(core::slice::from_ref(&training), lags, RETENTION);
        let theirs = tepsim_oracle::tier6::fit(core::slice::from_ref(&training), lags);
        assert_eq!(
            local, theirs,
            "the local fit diverged from tier6::fit at {lags} lags"
        );
    }
}

/// `d00.dat` is 500 rows and every `dNN_te.dat` is 960, with `XMV(12)`
/// supplied as a constant.
///
/// Asserted because every rate in this file is a fraction whose denominator is
/// one of those two numbers. A file that silently loaded short would move
/// every rate without failing anything else.
#[test]
fn the_published_files_have_the_shape_the_rates_assume() {
    let training = Published {
        fault: 0,
        split: Split::Training,
    }
    .run();
    assert_eq!(training.samples.len(), 500);
    for fault in 0..=21 {
        let run = Published {
            fault,
            split: Split::Testing,
        }
        .run();
        assert_eq!(
            run.samples.len(),
            TESTING_ROWS,
            "d{fault:02}_te is not 960 rows"
        );
    }
    // The 53rd channel is the unrecorded agitator, constant in every published
    // file, so PCA standardises it to zero and the fitted model is the
    // literature's 52-variable one. Checked rather than assumed, because if it
    // were not constant it would be a 53rd variable nobody asked for.
    let model = fit_with(core::slice::from_ref(&training), 0, RETENTION);
    assert_eq!(
        model.constant_columns(),
        &[VARIABLES - 1],
        "exactly the agitator column should be constant in d00"
    );
}

/// The headline: train on `d00`, evaluate on `d01_te` through `d21_te`, and
/// check the two qualitative facts the TEP literature agrees on.
#[test]
fn the_published_detection_rates_separate_the_easy_faults_from_the_famous_three() {
    let size = Size::selected();
    println!(
        "\nTier 6: detection rates on the published d00-d21\n\
         ================================================\n\
         trained on d00.dat (500 samples), evaluated on dNN_te.dat (960 samples,\n\
         fault at row {ONSET}); detection rate over the 800 post-fault samples,\n\
         false alarm rate over the 160 before it; {PERSISTENCE} consecutive\n\
         alarms for a delay. No published per-fault number is quoted as a\n\
         target: see the module docs.\n\
         size: {} ({}=full adds DPCA and the retention sweep)",
        if size.sweep { "full" } else { "default" },
        Size::ENV
    );

    let detector = Detector {
        lags: 0,
        retention: RETENTION,
        label: "PCA, 90% variance, 99% limits (the Tier 6 detector)",
    };
    let evaluation = evaluate(detector);
    print_table(detector, &evaluation);
    check(detector, &evaluation);

    println!(
        "\n  easy faults {EASY:?}\n  hard faults {HARD:?}\n  \
         every hard fault alarms at most {:.2}x as often as fault-free data; \
         every easy fault at least {:.2}x",
        evaluation
            .rates
            .iter()
            .filter(|r| HARD.contains(&r.fault))
            .map(|r| (r.t_squared / evaluation.nominal.0).max(r.spe / evaluation.nominal.1))
            .fold(0.0_f64, f64::max),
        evaluation
            .rates
            .iter()
            .filter(|r| EASY.contains(&r.fault))
            .map(|r| (r.t_squared / evaluation.nominal.0).max(r.spe / evaluation.nominal.1))
            .fold(f64::INFINITY, f64::min),
    );

    if !size.sweep {
        println!(
            "\n  Set {}=full for dynamic PCA and the retention sweep, which is \
             the sensitivity analysis.",
            Size::ENV
        );
        return;
    }

    // The sensitivity analysis. If the separation survives every retention
    // rule and both lag counts, it is a property of the data.
    println!("\n-- sensitivity: does the separation depend on the parameters? --");
    let sweep = [
        Detector {
            lags: 0,
            retention: Retention::CumulativeVariance(0.85),
            label: "PCA, 85% variance",
        },
        Detector {
            lags: 0,
            retention: Retention::CumulativeVariance(0.95),
            label: "PCA, 95% variance",
        },
        Detector {
            lags: 0,
            retention: Retention::CumulativeVariance(0.99),
            label: "PCA, 99% variance",
        },
        Detector {
            lags: 0,
            retention: Retention::Kaiser,
            label: "PCA, Kaiser's rule",
        },
        Detector {
            lags: 2,
            retention: RETENTION,
            label: "DPCA(2), 90% variance",
        },
        Detector {
            lags: 2,
            retention: Retention::CumulativeVariance(0.99),
            label: "DPCA(2), 99% variance",
        },
    ];
    let mut summary = Vec::new();
    for detector in sweep {
        let evaluation = evaluate(detector);
        print_table(detector, &evaluation);
        // The same two facts, at every parameter setting. The literature's
        // split is a property of the benchmark, not of a retention rule, so a
        // detector that loses it is wrong.
        check(detector, &evaluation);
        for s in separations(&evaluation.rates, evaluation.nominal) {
            summary.push((detector.label, evaluation.model.retained(), s));
        }
    }

    println!(
        "\n-- sensitivity summary --\n  {:<24} {:>5} {:>5} {:>10} {:>11} {:>8} {:>10} {:>8}",
        "detector", "kept", "stat", "worst easy", "best hard", "gap", "d00_te FAR", "carries"
    );
    for (label, kept, s) in &summary {
        println!(
            "  {:<24} {:>5} {:>5} {:>10.4} {:>11.4} {:>8.4} {:>10.4} {:>8}",
            label,
            kept,
            s.statistic.label(),
            s.worst_easy,
            s.best_hard,
            s.gap(),
            s.false_alarms,
            if s.carries_the_claim() { "yes" } else { "no" }
        );
    }
}
