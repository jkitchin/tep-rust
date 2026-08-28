# A fault that arrives, and clears

Every published Tennessee Eastman dataset is the same shape of experiment: set
some `IDV` flags before the run, then leave them. That is all the original
admits, which is why the literature's fault onsets are all at the same place and
why almost nobody studies a fault that arrives, persists, and then goes away.

A `Scenario` carries a `Schedule` of `Event`s, each at a time and each doing one
thing. This tutorial starts `IDV(4)` at hour six, stops it at hour twelve, and
then does something the original cannot express at all: runs it at half
strength.

```rust,ignore
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

fn hex(scenario: &Scenario) -> String {
    String::from_utf8(scenario.digest_hex().to_vec()).expect("hex is ascii")
}
```

```text
--- a fault that arrives at hour 6 and clears at hour 12 ---
events: 2
    6.0 h  Start { fault: 4 }
   12.0 h  Stop { fault: 4 }
    4.05 h   XMV(10)  41.063   active []
    5.95 h   XMV(10)  40.582   active []
    6.15 h   XMV(10)  45.254   active [4]
    8.05 h   XMV(10)  44.594   active [4, 12]
   11.95 h   XMV(10)  43.941   active [4, 12]
   12.15 h   XMV(10)  40.828   active [12]
   14.05 h   XMV(10)  39.513   active [12]
   17.05 h   XMV(10)  40.314   active [12]

--- ground truth follows the schedule ---
   5.95 h   faulted false   IDV(4) since onset None
   6.15 h   faulted true   IDV(4) since onset Some(0.15)
  11.95 h   faulted true   IDV(4) since onset Some(5.95)
  12.15 h   faulted true   IDV(4) since onset None

--- half a fault, which the original cannot express ---
    5.95 h   XMV(10)  40.582   active []
    8.05 h   XMV(10)  42.671   active [4, 12]
   11.95 h   XMV(10)  41.937   active [4, 12]
   14.05 h   XMV(10)  39.506   active [12]

--- without the extension it is refused, not rounded ---
  Err(ContinuousDisturbancesNotEnabled)

--- the schedule travels with the scenario ---
  tepsim.scenario.v1;seed=4651207995;hours=18;step=2.777777777777778e-4;every=180;faults=;controlled=1;idv12=1;trip=0;continuous=1;integrator=euler;events=6:magnitude:4:0.5,12:magnitude:4:0
  digest before 246b561b64789ada
  digest after  246b561b64789ada
  equal: true
```

## Reading the trace

The cooling water valve sits near 41% open, jumps to 45% within a sample of the
fault arriving, holds there for the six hours the fault is live, and drops back
to 41% within a sample of it clearing. The controller does not know a fault
happened in either direction; it is rejecting a disturbance both times.

An event is applied on the first step whose time is at or past its own, so an
event at hour six lands between the samples at 5.95 and 6.15. The window is
half-open at the start and closed at the end, so an event is applied exactly
once however the step size divides its time, and an event at time zero is
applied on the first step rather than never.

Two events at the same instant keep the order they were added, because `Stop`
then `Start` on one fault is a different scenario from the reverse, and the
digest has to be able to tell them apart.

## The labels are honest, including where that is inconvenient

At 12.15 hours `faulted()` still reports `true` while `since_onset[3]` has gone
back to `None`. Both are correct. `IDV(4)` really has cleared, and `IDV(12)` is
still on because the driver switched it on at hour eight, as `temain_mod.f`
does. `faulted()` asks whether anything is wrong and `since_onset[i]` asks about
one specific disturbance, and on a scheduled run those questions genuinely come
apart. Set `driver_forces_idv12` to `false` to study a schedule without the
driver's contribution.

## Half a fault

`teprob.f:341-346` opens `TEFUNC` by forcing every `IDV` to zero or one, so the
original has exactly two states per disturbance. The disturbances are then used
multiplicatively, `XST(1,4) = TESUB8(1,TIME) - IDV(1)*0.03` at `teprob.f:407`
being the pattern, so a magnitude between zero and one scales the fault smoothly
and needs no other change to the model.

That is the `continuous_disturbances` extension, and it is off by default. With
it on, `IDV(4)` at magnitude 0.5 moves the valve to 42.7% where the full fault
moves it to 44.6% and the fault-free plant sits at 41%, which is about what
halfway should look like.

Two things are worth saying plainly about it. A run with the extension on is not
comparable to any published dataset and none of the validation ladder applies to
it. And a fractional magnitude without the extension is *refused* rather than
rounded: `Scenario::validate` returns
`Err(Invalid::ContinuousDisturbancesNotEnabled)`. Silently turning a request for
half a fault into a whole one would produce a run that does not match its own
description, which is exactly what the content hash exists to prevent. A
magnitude of exactly 1.0 is bit-identical to the faithful path, which is
asserted, so the extension can be left on for a study that mixes full and
partial faults.

## The schedule is part of the description

The last block is the reason all of this is worth doing. The scenario, schedule
included, writes out as one line of text, that line parses back to an equal
scenario, and the content hash is unchanged: `246b561b64789ada` before and
after. A dataset generated from this scenario can carry the line in its header,
and a reader can reconstruct the exact run rather than the description of a run.

The same line is what TEP Studio puts in a URL fragment, so a scheduled
experiment can be handed to somebody as a link. The format is versioned and its
parser is strict: a missing field, an unknown field, a value out of range or a
version this build does not know is an error that says what was wrong, rather
than a default quietly substituted. See `tepsim::text`.

A schedule holds up to 32 events, which is far more than any experiment in the
literature uses. The published datasets have exactly one event each: the fault,
switched on before the run.
