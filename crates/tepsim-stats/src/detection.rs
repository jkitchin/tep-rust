//! Fault detection rate, false alarm rate, and detection delay.
//!
//! The three numbers Tier 6 compares between simulators. `PLAN.org` states the
//! claim: a detector trained on one simulator and tested on the other should
//! perform the way it does trained and tested on the same one. That comparison
//! is only meaningful if the three metrics are defined the same way on both
//! sides, and the literature does not define them the same way as each other,
//! so they are pinned here.
//!
//! # The convention these follow
//!
//! Every metric is over **samples**, not over runs. One fault file yields one
//! detection rate, computed as a fraction of the observations in it, exactly as
//! Chiang, Russell and Braatz report for TEP: their tables of missed detection
//! rates are `1 - FDR` in this module's terms, over the 800 post-fault
//! observations of each test file.
//!
//! Chiang, L. H., Russell, E. L. and Braatz, R. D. (2001), *Fault Detection and
//! Diagnosis in Industrial Systems*, Springer, chapter 10.
//!
//! # The benchmark's shape
//!
//! The published TEP test files `d01_te` through `d21_te` are 960 samples at a
//! three-minute interval: 160 fault-free samples, then the fault is introduced
//! and 800 more follow. `d00_te` is 960 fault-free samples throughout. So the
//! usual call is `onset = 160`, and a false alarm rate measured on `d00_te` is
//! the call with `onset = alarms.len()`.
//!
//! # Alarms, not statistics
//!
//! These take a boolean alarm series rather than a statistic and a limit, so
//! that a detector combining T-squared and SPE, or one with a persistence rule,
//! feeds the same functions. [`alarms_above`] builds the simple series.

use alloc::vec::Vec;

/// Threshold a statistic into an alarm series: `statistic > limit`.
///
/// **Strictly** greater. A statistic landing exactly on its control limit is
/// not an exceedance. The case is unreachable on real data and reachable on
/// constructed data, and picking the side deliberately is cheaper than
/// discovering later which side was picked.
#[must_use]
pub fn alarms_above(statistic: &[f64], limit: f64) -> Vec<bool> {
    // A `NaN` statistic compares false and so raises no alarm. That is the
    // conservative direction and it is not silent: a detector producing `NaN`
    // will show a detection rate of zero, which is not a result anyone reads
    // past.
    statistic.iter().map(|&s| s > limit).collect()
}

/// Everything the three metrics report, computed in one pass over the series.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DetectionReport {
    /// The index the fault was introduced at.
    pub onset: usize,
    /// How many samples in the series.
    pub samples: usize,
    /// Samples before the onset.
    pub pre_fault: usize,
    /// Samples at or after the onset.
    pub post_fault: usize,
    /// Alarms among the pre-fault samples.
    pub false_alarms: usize,
    /// Alarms among the post-fault samples.
    pub detections: usize,
    /// [`fault_detection_rate`].
    pub fault_detection_rate: f64,
    /// [`false_alarm_rate`].
    pub false_alarm_rate: f64,
    /// [`detection_delay`], in samples.
    pub detection_delay: Option<usize>,
    /// The persistence the delay was measured with.
    pub consecutive: usize,
}

/// The fraction of post-fault samples that raised an alarm.
///
/// ```text
/// FDR = |{ i >= onset : alarms[i] }| / |{ i >= onset }|
/// ```
///
/// The sample at index `onset` is the first one **with** the fault present, so
/// it counts as post-fault. Its complement `1 - FDR` is the missed detection
/// rate the TEP tables report.
///
/// This is a rate over samples and not a per-fault yes-or-no. A detector that
/// catches a fault immediately and then loses it scores badly here, which is
/// the intent: a monitoring statistic that drops back inside its limit while
/// the fault is still running is not detecting the fault, it is flickering.
///
/// Returns `NaN` when there are no post-fault samples, rather than zero. Zero
/// is a detector that failed; no data is not.
#[must_use]
pub fn fault_detection_rate(alarms: &[bool], onset: usize) -> f64 {
    if onset >= alarms.len() {
        return f64::NAN;
    }
    let post = &alarms[onset..];
    post.iter().filter(|&&a| a).count() as f64 / post.len() as f64
}

/// The fraction of pre-fault samples that raised an alarm.
///
/// ```text
/// FAR = |{ i < onset : alarms[i] }| / |{ i < onset }|
/// ```
///
/// Measured on the fault-free part of the same record, so it is the detector's
/// type-I error rate under the operating conditions it is actually being asked
/// about, rather than the nominal `1 - confidence` its control limit was drawn
/// at. The two differ, and the size of the gap is the interesting number: a
/// limit set at 99% confidence that produces 4% false alarms is telling you the
/// data are not what the limit's distributional assumption says they are.
///
/// Returns `NaN` when `onset` is zero, because a record with no fault-free part
/// carries no information about false alarms.
#[must_use]
pub fn false_alarm_rate(alarms: &[bool], onset: usize) -> f64 {
    let boundary = onset.min(alarms.len());
    if boundary == 0 {
        return f64::NAN;
    }
    let pre = &alarms[..boundary];
    pre.iter().filter(|&&a| a).count() as f64 / boundary as f64
}

/// How many samples after the onset the fault was first detected.
///
/// The smallest `d >= 0` such that samples `onset + d` through
/// `onset + d + consecutive - 1` all raised an alarm, or `None` if no such run
/// exists inside the record. In **samples**; multiply by the sampling interval
/// for a time. The TEP files are three minutes apart, so a delay of 5 is
/// fifteen minutes.
///
/// A delay of zero means the very first post-fault sample started the run.
///
/// # Why a persistence requirement
///
/// With `consecutive = 1` this is the first alarm after the onset, and on a
/// detector with any false alarm rate at all that is mostly a measurement of
/// luck: at a 1% false alarm rate the first post-fault sample alarms by chance
/// one time in a hundred, and reporting a delay of zero for it flatters the
/// detector. Requiring a run of alarms is the standard remedy, and the run
/// length is a parameter because the literature does not agree on it: three and
/// six both appear, and the choice moves the answer.
///
/// The run must fit inside the record. A fault detected on the last two samples
/// with `consecutive = 3` is not detected, which is the honest reading: the
/// evidence for it does not exist.
///
/// Returns `None` when there are no post-fault samples.
///
/// # Panics
///
/// If `consecutive` is zero. There is no sensible reading of "detected after a
/// run of no alarms".
#[must_use]
pub fn detection_delay(alarms: &[bool], onset: usize, consecutive: usize) -> Option<usize> {
    assert!(consecutive > 0, "a detection needs at least one alarm");
    if onset >= alarms.len() {
        return None;
    }
    let post = &alarms[onset..];
    if post.len() < consecutive {
        return None;
    }
    (0..=(post.len() - consecutive))
        .find(|&start| post[start..start + consecutive].iter().all(|&a| a))
}

/// All three metrics and the counts behind them.
///
/// The counts are reported alongside the rates because a rate without its
/// denominator cannot be compared against the next run, which is `CLAUDE.md`'s
/// rule about numbers rather than verdicts applied to a detector.
///
/// # Panics
///
/// If `consecutive` is zero.
#[must_use]
pub fn detection_report(alarms: &[bool], onset: usize, consecutive: usize) -> DetectionReport {
    let boundary = onset.min(alarms.len());
    DetectionReport {
        onset,
        samples: alarms.len(),
        pre_fault: boundary,
        post_fault: alarms.len() - boundary,
        false_alarms: alarms[..boundary].iter().filter(|&&a| a).count(),
        detections: alarms[boundary..].iter().filter(|&&a| a).count(),
        fault_detection_rate: fault_detection_rate(alarms, onset),
        false_alarm_rate: false_alarm_rate(alarms, onset),
        detection_delay: detection_delay(alarms, onset, consecutive),
        consecutive,
    }
}
