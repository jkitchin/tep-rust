//! The worked example behind `book/src/tutorials/scheduling-a-fault.md`.
//!
//! The book shows this file's body verbatim, and
//! `crates/tepsim/tests/book_examples.rs` asserts the two are the same bytes,
//! so a listing in the book cannot drift away from code that compiles. Its
//! printed output is the transcript the tutorial quotes.

use tepsim::tepsim_core::Extensions;
use tepsim::{Action, Event, Run, Scenario, Simulation};

/// The value of one manipulated variable at a few chosen hours.
fn trace(run: &Run, mv: usize, at: &[f64]) {
    let series = run.manipulated(mv);
    for hour in at {
        let index = run
            .samples
            .iter()
            .position(|s| s.hours >= *hour)
            .unwrap_or(run.samples.len() - 1);
        let faults: Vec<usize> = run.samples[index].labels.faults().collect();
        println!(
            "  {:>6.2} h   XMV({mv}) {:>7.3}   active {faults:?}",
            run.samples[index].hours, series[index]
        );
    }
}

/// Runs the four demonstrations the tutorial walks through.
fn main() {
    println!("--- a fault that arrives at hour 6 and clears at hour 12 ---");
    let scenario = Scenario::baseline()
        .with_hours(18.0)
        .with_event(Event::start(6.0, 4))
        .with_event(Event::stop(12.0, 4));
    println!("events: {}", scenario.schedule.len());
    for event in scenario.schedule.events() {
        println!("  {:>5.1} h  {:?}", event.at_hours, event.action);
    }
    let run = Simulation::new(scenario).run();
    trace(&run, 10, &[4.0, 5.9, 6.1, 8.0, 11.9, 12.1, 14.0, 17.0]);

    println!();
    println!("--- ground truth follows the schedule ---");
    for hour in [5.9, 6.1, 11.9, 12.1] {
        let sample = run
            .samples
            .iter()
            .find(|s| s.hours >= hour)
            .expect("the run is long enough");
        println!(
            "  {:>5.2} h   faulted {}   IDV(4) since onset {:?}",
            sample.hours,
            sample.labels.faulted(),
            sample.labels.since_onset[3].map(|h| (h * 100.0).round() / 100.0)
        );
    }

    println!();
    println!("--- half a fault, which the original cannot express ---");
    let partial = Scenario::baseline()
        .with_hours(18.0)
        .with_continuous_disturbances()
        .with_event(Event::new(
            6.0,
            Action::SetMagnitude {
                fault: 4,
                magnitude: 0.5,
            },
        ))
        .with_event(Event::new(
            12.0,
            Action::SetMagnitude {
                fault: 4,
                magnitude: 0.0,
            },
        ));
    let half = Simulation::new(partial).run();
    trace(&half, 10, &[5.9, 8.0, 11.9, 14.0]);

    println!();
    println!("--- without the extension it is refused, not rounded ---");
    let refused = Scenario {
        extensions: Extensions::none(),
        ..partial
    };
    println!("  {:?}", refused.validate());

    println!();
    println!("--- the schedule travels with the scenario ---");
    let text = partial.to_text();
    println!("  {text}");
    let back = Scenario::from_text(&text).expect("its own text parses");
    println!("  digest before {}", hex(&partial));
    println!("  digest after  {}", hex(&back));
    println!("  equal: {}", back == partial);
}

/// A scenario's digest, as the sixteen hex characters a URL fragment carries.
fn hex(scenario: &Scenario) -> String {
    String::from_utf8(scenario.digest_hex().to_vec()).expect("hex is ascii")
}
