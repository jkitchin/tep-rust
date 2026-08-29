//! Scheduled events, continuous magnitudes, and the content hash.
//!
//! The original admits one shape of experiment: set some flags before the run
//! and leave them. Every published Tennessee Eastman dataset is that. These
//! check the things it cannot express.

#![allow(
    clippy::float_cmp,
    reason = "the run is deterministic, so equality is exact or the test fails"
)]

use tepsim::{Action, Event, Invalid, Outcome, Scenario, Schedule, Simulation};

const SHORT: f64 = 1.0;

// ---------------------------------------------------------------------------
// The schedule itself
// ---------------------------------------------------------------------------

#[test]
fn events_are_kept_in_time_order() {
    let schedule = Schedule::new()
        .with(Event::start(3.0, 1))
        .with(Event::stop(1.0, 2))
        .with(Event::start(2.0, 3));

    let times: Vec<f64> = schedule.events().map(|e| e.at_hours).collect();
    assert_eq!(times, vec![1.0, 2.0, 3.0]);
    assert_eq!(schedule.len(), 3);
}

/// Two events at the same instant keep the order they were added.
///
/// `Stop` then `Start` on one fault at one time is a different experiment from
/// `Start` then `Stop`, and a schedule that reordered them would silently turn
/// one into the other.
#[test]
fn events_at_the_same_instant_keep_their_insertion_order() {
    let forward = Schedule::new()
        .with(Event::stop(5.0, 4))
        .with(Event::start(5.0, 4));
    let backward = Schedule::new()
        .with(Event::start(5.0, 4))
        .with(Event::stop(5.0, 4));

    let first = |s: &Schedule| s.events().next().expect("an event").action;
    assert_eq!(first(&forward), Action::Stop { fault: 4 });
    assert_eq!(first(&backward), Action::Start { fault: 4 });
    assert_ne!(forward, backward);
}

/// An event is applied exactly once, whatever the step size.
#[test]
fn each_event_fires_exactly_once() {
    let scenario = Scenario::baseline()
        .with_hours(SHORT)
        .with_event(Event::start(0.25, 4))
        .with_event(Event::stop(0.5, 4));
    let run = Simulation::new(scenario).run();

    // Before 0.25 h: off. Between 0.25 and 0.5: on. After 0.5: off again.
    let state_at = |hours: f64| {
        run.samples
            .iter()
            .rfind(|s| s.hours <= hours)
            .map(|s| s.labels.active[3])
    };
    assert_eq!(state_at(0.2), Some(false), "fired early");
    assert_eq!(state_at(0.4), Some(true), "did not fire");
    assert_eq!(state_at(0.9), Some(false), "did not clear");
}

/// A fault that stops and starts again has a *new* onset, not a continuing one.
#[test]
fn a_restarted_fault_gets_a_fresh_onset() {
    let scenario = Scenario::baseline()
        .with_hours(SHORT)
        .with_event(Event::start(0.1, 6))
        .with_event(Event::stop(0.3, 6))
        .with_event(Event::start(0.6, 6));
    let run = Simulation::new(scenario).run();

    let last = run.samples.last().expect("samples");
    assert!(last.labels.active[5]);
    let since = last.labels.since_onset[5].expect("an onset");

    // Measured from the second onset at 0.6 h, not from the first at 0.1 h.
    //
    // Within one step, not exactly: an event at 0.6 h fires on the first step
    // whose time has reached 0.6, and the onset is stamped with that step's
    // time. So the recorded onset is in `[0.6, 0.6 + dt)` and the elapsed time
    // is short by up to one step. Asking for exact equality would be asking
    // the scheduler to interpolate within a step, which it deliberately does
    // not: an event that took effect part way through a step would make the
    // run depend on the step size in a way nothing else does.
    let expected = last.hours - 0.6;
    let step = run.scenario.step_hours;
    assert!(
        (since - expected).abs() <= step,
        "since_onset is {since:.6} at {:.6} h; measured from the second onset \
         it should be about {expected:.6}, and from the first it would be \
         {:.6}",
        last.hours,
        last.hours - 0.1
    );
    // And it is emphatically not measured from the first onset.
    assert!(
        (since - (last.hours - 0.1)).abs() > 10.0 * step,
        "since_onset appears to be measured from the first onset"
    );
}

/// A scheduled fault actually changes the plant.
#[test]
fn a_scheduled_fault_moves_the_plant() {
    let clean = Simulation::new(Scenario::baseline().with_hours(SHORT)).run();
    let scheduled = Simulation::new(
        Scenario::baseline()
            .with_hours(SHORT)
            .with_event(Event::start(0.2, 1)),
    )
    .run();

    // Identical before the event, different after.
    let before = |run: &tepsim::Run| {
        run.samples
            .iter()
            .filter(|s| s.hours < 0.2)
            .map(tepsim::Sample::row)
            .collect::<Vec<_>>()
    };
    assert_eq!(
        before(&clean),
        before(&scheduled),
        "differed before the event"
    );

    let worst = clean
        .samples
        .iter()
        .zip(&scheduled.samples)
        .filter(|(a, _)| a.hours > 0.5)
        .flat_map(|(a, b)| a.row().into_iter().zip(b.row()))
        .filter(|(x, _)| *x != 0.0)
        .map(|(x, y)| (x - y).abs() / x.abs())
        .fold(0.0_f64, f64::max);
    println!("IDV(1) from 0.2 h: worst departure {worst:.3e}");
    assert!(worst > 1e-4, "the scheduled fault did nothing: {worst:.3e}");
}

// ---------------------------------------------------------------------------
// Continuous magnitudes
// ---------------------------------------------------------------------------

/// A fractional magnitude is refused unless the extension is on.
///
/// Refused, not rounded. Silently turning a request for half a fault into a
/// whole one would produce a run that does not match its own description,
/// which is exactly what the content hash exists to prevent.
#[test]
fn a_fractional_magnitude_needs_the_extension() {
    let half = Event::new(
        1.0,
        Action::SetMagnitude {
            fault: 4,
            magnitude: 0.5,
        },
    );

    let refused = Scenario::baseline().with_event(half);
    assert_eq!(
        refused.validate(),
        Err(Invalid::ContinuousDisturbancesNotEnabled)
    );

    let allowed = refused.with_continuous_disturbances();
    assert_eq!(allowed.validate(), Ok(()));

    // Exactly 0 and exactly 1 are what the original admits, so they need
    // nothing.
    for magnitude in [0.0, 1.0] {
        let scenario = Scenario::baseline().with_event(Event::new(
            1.0,
            Action::SetMagnitude {
                fault: 4,
                magnitude,
            },
        ));
        assert_eq!(scenario.validate(), Ok(()), "magnitude {magnitude}");
    }
}

/// A magnitude of exactly one is bit-identical to the faithful path.
///
/// This is what lets the extension stay on for a study that mixes full and
/// partial faults: turning it on does not by itself change anything.
#[test]
fn magnitude_one_is_identical_to_the_faithful_path() {
    let faithful = Simulation::new(Scenario::fault(4).with_hours(SHORT)).run();
    let extended = Simulation::new(
        Scenario::baseline()
            .with_hours(SHORT)
            .with_continuous_disturbances()
            .with_event(Event::new(
                0.0,
                Action::SetMagnitude {
                    fault: 4,
                    magnitude: 1.0,
                },
            )),
    )
    .run();

    assert_eq!(faithful.samples.len(), extended.samples.len());
    for (index, (a, b)) in faithful.samples.iter().zip(&extended.samples).enumerate() {
        assert_eq!(a.row(), b.row(), "sample {index}");
    }
}

/// Half a fault sits between no fault and a whole one.
#[test]
fn a_half_magnitude_fault_sits_between_none_and_all() {
    let at = |magnitude: f64| {
        let mut scenario = Scenario::baseline()
            .with_hours(SHORT)
            .with_continuous_disturbances();
        if magnitude > 0.0 {
            scenario = scenario.with_event(Event::new(
                0.0,
                Action::SetMagnitude {
                    fault: 1,
                    magnitude,
                },
            ));
        }
        // XMEAS(1) is the A feed, which IDV(1) acts on directly.
        let run = Simulation::new(scenario).run();
        let series = run.measurement(1);
        series.iter().sum::<f64>() / series.len() as f64
    };

    let none = at(0.0);
    let half = at(0.5);
    let full = at(1.0);
    println!("mean XMEAS(1): none {none:.6}, half {half:.6}, full {full:.6}");

    assert!(
        (half - none).abs() > 1e-6,
        "half a fault is indistinguishable from none"
    );
    assert!(
        (half - full).abs() > 1e-6,
        "half a fault is indistinguishable from a whole one"
    );
    let between = (half - none) / (full - none);
    println!("half sits at {between:.3} of the way from none to full");
    assert!(
        (0.2..=0.8).contains(&between),
        "half a fault sits at {between:.3} of the way, which is not between"
    );
}

// ---------------------------------------------------------------------------
// Validation
// ---------------------------------------------------------------------------

#[test]
fn a_scenario_is_checked_rather_than_trusted() {
    let cases: &[(Event, Invalid)] = &[
        (Event::start(1.0, 0), Invalid::FaultOutOfRange { fault: 0 }),
        (
            Event::start(1.0, 21),
            Invalid::FaultOutOfRange { fault: 21 },
        ),
        (Event::start(-1.0, 4), Invalid::TimeNotFinite),
        (Event::start(f64::NAN, 4), Invalid::TimeNotFinite),
        (
            Event::new(
                1.0,
                Action::SetMagnitude {
                    fault: 4,
                    magnitude: 1.5,
                },
            ),
            Invalid::MagnitudeOutOfRange,
        ),
        (
            Event::new(
                1.0,
                Action::Setpoint {
                    loop_index: 0,
                    value: 1.0,
                },
            ),
            Invalid::LoopOutOfRange { loop_index: 0 },
        ),
    ];
    for (event, expected) in cases {
        let scenario = Scenario::baseline()
            .with_continuous_disturbances()
            .with_event(*event);
        assert_eq!(scenario.validate(), Err(*expected), "for {event:?}");
    }
    assert_eq!(Scenario::baseline().validate(), Ok(()));
}

// ---------------------------------------------------------------------------
// The content hash
// ---------------------------------------------------------------------------

/// Every field that affects the run affects the digest.
///
/// A digest that missed a field would let two different experiments claim the
/// same identity, which is the one failure this is meant to make impossible.
#[test]
fn the_digest_covers_every_field_that_changes_a_run() {
    let base = Scenario::baseline();
    let variants: Vec<(&str, Scenario)> = vec![
        ("seed", base.with_seed(1.0)),
        ("hours", base.with_hours(1.0)),
        ("sampling", base.sampling_every(60)),
        ("fault", base.with_fault(3)),
        ("open loop", base.open_loop()),
        ("integrator", base.with_integrator(tepsim::Integrator::Rk4)),
        ("event", base.with_event(Event::start(1.0, 5))),
        ("extension", base.with_continuous_disturbances()),
        ("quirk", {
            let mut s = base;
            s.quirks.trip_ends_the_run = false;
            s
        }),
        ("forced idv12", {
            let mut s = base;
            s.driver_forces_idv12 = true;
            s
        }),
    ];

    let reference = base.digest();
    for (what, variant) in &variants {
        assert_ne!(
            variant.digest(),
            reference,
            "changing the {what} did not change the digest"
        );
    }
    // And they are all distinct from each other, not merely from the base.
    for (i, (a_name, a)) in variants.iter().enumerate() {
        for (b_name, b) in variants.iter().skip(i + 1) {
            assert_ne!(a.digest(), b.digest(), "{a_name} collides with {b_name}");
        }
    }

    // Equal scenarios hash equally, including ones built differently.
    assert_eq!(
        Scenario::fault(7).with_hours(3.0).digest(),
        Scenario::baseline().with_hours(3.0).with_fault(7).digest()
    );
}

#[test]
fn the_digest_is_stable_and_printable() {
    let hex = Scenario::baseline().digest_hex();
    let text = core::str::from_utf8(&hex).expect("ascii");
    assert_eq!(text.len(), 16);
    assert!(text.chars().all(|c| c.is_ascii_hexdigit()));
    // Stable across calls and across equal values.
    assert_eq!(Scenario::baseline().digest_hex(), hex);
    println!("baseline scenario digest: {text}");
}

/// Event order matters to the digest, because it matters to the run.
#[test]
fn the_digest_distinguishes_event_order() {
    let a = Scenario::baseline()
        .with_event(Event::stop(5.0, 4))
        .with_event(Event::start(5.0, 4));
    let b = Scenario::baseline()
        .with_event(Event::start(5.0, 4))
        .with_event(Event::stop(5.0, 4));
    assert_ne!(a.digest(), b.digest());
}

// ---------------------------------------------------------------------------
// Nothing above changed the faithful path
// ---------------------------------------------------------------------------

#[test]
fn an_empty_schedule_leaves_the_run_untouched() {
    let run = Simulation::new(Scenario::baseline().with_hours(SHORT)).run();
    assert_eq!(run.outcome, Outcome::Completed);
    assert!(run.scenario.schedule.is_empty());
    // Reactor pressure where Downs and Vogel say it is.
    assert!(
        run.measurement(7)
            .iter()
            .all(|p| (2600.0..2800.0).contains(p))
    );
}
