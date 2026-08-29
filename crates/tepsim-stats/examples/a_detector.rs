//! The worked example behind `book/src/tutorials/a-detector.md`.
//!
//! The book shows this file's body verbatim, and
//! `crates/tepsim/tests/book_examples.rs` asserts the two are the same bytes
//! and that the transcript beside it is the output this produces.

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
