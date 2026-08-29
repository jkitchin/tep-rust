//! The worked example behind `book/src/tutorials/injecting-a-fault.md`.
//!
//! The book shows this file's body verbatim, and
//! `crates/tepsim/tests/book_examples.rs` asserts the two are the same bytes
//! and that the transcript beside it is the output this produces.

use tepsim::{Scenario, Simulation};

/// Mean and the largest absolute departure from a reference series.
fn compare(name: &str, base: &[f64], faulted: &[f64]) {
    let mean = |x: &[f64]| x.iter().sum::<f64>() / x.len() as f64;
    let worst = base
        .iter()
        .zip(faulted)
        .map(|(b, f)| (f - b).abs())
        .fold(0.0_f64, f64::max);
    println!(
        "  {name:<28} base {:>9.3}   fault {:>9.3}   worst gap {:>8.3}",
        mean(base),
        mean(faulted),
        worst
    );
}

fn main() {
    let hours = 8.0;
    let base = Simulation::new(Scenario::baseline().with_hours(hours)).run();
    let faulted = Simulation::new(Scenario::fault(4).with_hours(hours)).run();

    println!("--- IDV(4): reactor cooling water inlet temperature, step ---");
    println!("samples: {} each", base.samples.len());
    compare(
        "XMEAS(9) reactor temp",
        &base.measurement(9),
        &faulted.measurement(9),
    );
    compare(
        "XMEAS(21) reactor cw out",
        &base.measurement(21),
        &faulted.measurement(21),
    );
    compare(
        "XMV(10) reactor cw valve",
        &base.manipulated(10),
        &faulted.manipulated(10),
    );

    println!();
    println!("--- the first hour, sample by sample ---");
    let temp = faulted.measurement(9);
    let valve = faulted.manipulated(10);
    for i in [0, 1, 2, 3, 4, 6, 9, 14, 19] {
        println!(
            "  {:>6.3} h   XMEAS(9) {:>8.4}   XMV(10) {:>7.3}",
            faulted.samples[i].hours, temp[i], valve[i]
        );
    }

    println!();
    println!("--- ground truth ---");
    let first = &faulted.samples[0];
    println!("  faulted:      {}", first.labels.faulted());
    println!(
        "  faults:       {:?}",
        first.labels.faults().collect::<Vec<_>>()
    );
    println!("  since onset:  {:?} h", first.labels.since_onset[3]);

    println!();
    println!("--- the driver's IDV(12), which is no longer the default ---");
    let forced = Simulation::new(Scenario::faithful().with_hours(9.0).with_fault(4)).run();
    for sample in &forced.samples {
        let faults: Vec<usize> = sample.labels.faults().collect();
        if faults.len() > 1 {
            println!(
                "  faithful: first sample with two faults at {:.3} h {faults:?}",
                sample.hours
            );
            break;
        }
    }
    let plain = Simulation::new(Scenario::fault(4).with_hours(9.0)).run();
    let ever: Vec<usize> = plain
        .samples
        .last()
        .expect("a run has samples")
        .labels
        .faults()
        .collect();
    println!("  default:  at 9 h, still only {ever:?}");
}
