# Building a detector, and measuring it

`tepsim-stats` carries the statistics the validation ladder runs on, and they
are equally usable as a monitoring scheme. This tutorial fits a PCA model to a
day of the fault-free plant, turns Hotelling's T-squared and the squared
prediction error into alarms, and then measures what those alarms are worth.

`tepsim-stats` is a development-only crate. Nothing shipped depends on it and
`cargo xtask ci` asserts that, so add it as a `dev-dependency` or use it from a
separate crate of your own.

The measurement is the point of the tutorial. Any detector can be made to look
good by reporting only the faults it catches, so the program below reports the
false alarm rate on held-out fault-free data alongside every detection rate, and
one of the two statistics comes out badly.

```rust,ignore
use tepsim::run::CHANNELS;
use tepsim::{Run, Scenario, Simulation};
use tepsim_stats::{Pca, Retention, alarms_above, detection_report};

/// A run flattened into the row-major matrix `Pca::fit` wants: `samples` rows
/// of `CHANNELS` values, which is exactly `Sample::row` stacked.
fn matrix(run: &Run) -> Vec<f64> {
    let mut out = Vec::with_capacity(run.samples.len() * CHANNELS);
    for sample in &run.samples {
        out.extend_from_slice(&sample.row());
    }
    out
}

/// The fault-free plant, with the driver's hour-eight IDV(12) switched off so
/// that a long fault-free record really is fault free.
fn clean(hours: f64, seed: f64) -> Scenario {
    Scenario {
        driver_forces_idv12: false,
        ..Scenario::baseline()
    }
    .with_hours(hours)
    .with_seed(seed)
}

/// Both monitoring statistics for every sample of a run.
fn monitor(model: &Pca, run: &Run) -> (Vec<f64>, Vec<f64>) {
    let mut t2 = Vec::with_capacity(run.samples.len());
    let mut spe = Vec::with_capacity(run.samples.len());
    for sample in &run.samples {
        let s = model.statistics(&sample.row());
        t2.push(s.t_squared);
        spe.push(s.spe);
    }
    (t2, spe)
}

fn main() {
    let hours = 24.0;

    let training = Simulation::new(clean(hours, 4_651_207_995.0)).run();
    let model = Pca::fit(
        &matrix(&training),
        training.samples.len(),
        CHANNELS,
        Retention::CumulativeVariance(0.9),
    );
    let limits = model.limits(0.99);

    println!("--- the model ---");
    println!("training samples:   {}", training.samples.len());
    println!("variables:          {CHANNELS}");
    println!("constant columns:   {:?}", model.constant_columns());
    println!("components kept:    {}", limits.components);
    println!("variance explained: {:.4}", model.explained_variance());
    println!("T-squared limit:    {:.3}", limits.t_squared);
    println!("SPE limit:          {:.3}", limits.spe);

    // A second fault-free record, from a different seed, is the pre-fault half
    // of every test below. It has to be data the model did not see, or the
    // false alarm rate is a measurement of the fit and not of the detector.
    let free = Simulation::new(clean(hours, 1_234_567_891.0)).run();
    let (free_t2, free_spe) = monitor(&model, &free);
    let onset = free.samples.len();

    println!();
    println!("--- false alarms, on {onset} fault-free samples ---");
    let free_t2_alarms = alarms_above(&free_t2, limits.t_squared);
    let free_spe_alarms = alarms_above(&free_spe, limits.spe);
    println!(
        "  T-squared: {} alarms, FAR {:.4}",
        free_t2_alarms.iter().filter(|a| **a).count(),
        free_t2_alarms.iter().filter(|a| **a).count() as f64 / onset as f64,
    );
    println!(
        "  SPE:       {} alarms, FAR {:.4}",
        free_spe_alarms.iter().filter(|a| **a).count(),
        free_spe_alarms.iter().filter(|a| **a).count() as f64 / onset as f64,
    );

    println!();
    println!("--- detection, four disturbances ---");
    println!("  fault     T2 FDR  T2 delay      Q FDR   Q delay");
    for fault in [1, 3, 4, 11] {
        let faulted = Simulation::new(clean(hours, 1_234_567_891.0).with_fault(fault)).run();
        let (fault_t2, fault_spe) = monitor(&model, &faulted);

        // The record the literature builds: fault-free samples, then faulted
        // ones, with the onset at the join.
        let mut t2 = free_t2.clone();
        t2.extend_from_slice(&fault_t2);
        let mut spe = free_spe.clone();
        spe.extend_from_slice(&fault_spe);

        let t2_report = detection_report(&alarms_above(&t2, limits.t_squared), onset, 3);
        let spe_report = detection_report(&alarms_above(&spe, limits.spe), onset, 3);
        let delay = |d: Option<usize>| d.map_or("never".into(), |d| format!("{d}"));
        println!(
            "  IDV({fault:>2})   {:>7.3}  {:>8}   {:>8.3}  {:>8}",
            t2_report.fault_detection_rate,
            delay(t2_report.detection_delay),
            spe_report.fault_detection_rate,
            delay(spe_report.detection_delay),
        );
    }
}
```

```text
--- the model ---
training samples:   480
variables:          53
constant columns:   [52]
components kept:    33
variance explained: 0.9110
T-squared limit:    60.386
SPE limit:          10.481

--- false alarms, on 480 fault-free samples ---
  T-squared: 10 alarms, FAR 0.0208
  SPE:       98 alarms, FAR 0.2042

--- detection, four disturbances ---
  fault     T2 FDR  T2 delay      Q FDR   Q delay
  IDV( 1)     0.994         3      0.998         1
  IDV( 3)     0.025     never      0.210       182
  IDV( 4)     0.502        29      0.998         1
  IDV(11)     0.492        18      0.731         8
```

## What the model did

`Pca::fit` standardises each column to zero mean and unit variance, forms the
correlation matrix, and diagonalises it with a cyclic Jacobi sweep, which is
deterministic to the last bit on every architecture. `Retention` names the rule
for how many components to keep, and it is a named enum rather than a bare `k`
because two detectors that retain different numbers of components are different
detectors: a report that does not say which rule produced it cannot be
reproduced. Here `CumulativeVariance(0.9)` kept 33 of 53 components, which
together explain 91.1% of the variance.

`constant_columns()` reports `[52]`, the zero-based position of `XMV(12)`. The
agitator speed never moves, so its standard deviation is zero and it cannot be
standardised. Rather than divide by zero or quietly drop the column,
`tepsim-stats` records it and excludes it from the model, and says so when
asked. A silent drop here is how a 53-variable model becomes a 52-variable model
that nobody can reproduce.

## What the alarms are worth

`alarms_above` thresholds a statistic into a boolean series, strictly greater
than the limit. `detection_report` then computes three numbers from that series
and an onset index: the fault detection rate over the post-onset samples, the
false alarm rate over the pre-onset ones, and the detection delay, which is the
first run of `consecutive` alarms after the onset. Three was used above.

The persistence requirement matters. With `consecutive = 1` the delay is the
first alarm after the onset, and on a detector with any false alarm rate at all
that is mostly luck: at a 2% false alarm rate the first post-onset sample
alarms by chance one time in fifty, and calling that a delay of zero flatters
the detector. The literature uses three and six and does not agree, so the run
length is a parameter and the report carries the value it was measured with.

`IDV(1)`, the A/C feed ratio step, is caught by both statistics almost
immediately and stays caught: a detection rate of 0.994 means the alarm was up
for essentially the whole record rather than flickering. `IDV(3)`, the D feed
temperature step, is caught by neither, which is the expected answer: it is one
of the three faults the TEP literature reports as effectively undetectable by
PCA monitoring, and a detection rate of 0.025 against a false alarm rate of
0.021 is a detector producing noise. `IDV(4)` is the case the previous tutorial
looked at by hand: the residual subspace sees it at once, and T-squared, which
lives in the subspace the training data actually spans, only catches it half the
time.

Delays are in samples. The cadence is 180 steps, so one sample is three minutes
and the T-squared delay of 3 on `IDV(1)` is nine minutes. The 182 on `IDV(3)`
is nine hours, which is not a detection.

## The number the detector would rather you did not see

The SPE limit was drawn at 99% confidence and produced alarms on 20.4% of
held-out fault-free samples. That is a factor of twenty, and it is not a bug in
`spe_limit`, which computes exactly the Jackson and Mudholkar expression it
claims to.

It is the assumption underneath that expression failing. The limit is derived
for residuals that are normal and independent, and the plant's are neither. The
feed compositions are driven by slow random walks, so a 24-hour record wanders
somewhere a different 24-hour record did not, and every sample in the excursion
alarms together. The T-squared limit, drawn from an F distribution over the
retained subspace, holds up much better at 2.1% against a nominal 1%.

This is the reason `tepsim-stats` reports counts alongside rates everywhere: 98
alarms out of 480 can be compared against the next run and "the detector had a
high false alarm rate" cannot. The fix, if you want one, is more training data,
a limit estimated empirically from the training residuals rather than
analytically, or `dpca`, which augments the matrix with lagged copies and models
the serial correlation instead of assuming it away. That is a different tutorial
and a real research question, which is rather the point of the simulator.
