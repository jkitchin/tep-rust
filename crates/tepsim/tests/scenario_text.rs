//! The canonical scenario text: it round-trips, it is strict, and nothing
//! about a scenario can be added without it noticing.
//!
//! Three properties are being defended here, and they are not the same
//! property.
//!
//! *Round trip.* `from_text(to_text(s)) == s`, with every `f64` equal by bit
//! pattern rather than by value, over the edge cases that a hand-written number
//! format actually gets wrong: negative zero, subnormals, the seed
//! `teprob.f:1187` compiles in, an empty schedule and a full one.
//!
//! *The digest survives.* That is what makes a serialised scenario worth
//! anything. A dataset labelled with its scenario's digest and shipped
//! alongside the scenario text can be checked; without this property it can
//! only be believed.
//!
//! *Nothing is silently defaulted.* Every way of writing a text wrong produces
//! an error that names the thing that was wrong. The reason is the bug this
//! module was written to kill: a field added to `Scenario` that quietly failed
//! to reach a shared link, so the link opened, the run ran, and it was a
//! different run from the one the link claimed.

use tepsim::text::{self, TextError};
use tepsim::{Action, Event, Integrator, Scenario, Schedule};
use tepsim_scenario::{Invalid, MAX_EVENTS};

/// The exact text of the baseline, as a golden.
///
/// Written out rather than computed, so a change to the format, to a key name,
/// to the field order or to how a number is spelled shows up here as a diff
/// and has to be a decision. Every link and every dataset header in circulation
/// is this string with values substituted.
const BASELINE_TEXT: &str = "tepsim.scenario.v1;seed=4651207995;hours=48;\
     step=2.777777777777778e-4;every=180;faults=;controlled=1;idv12=1;trip=0;\
     continuous=0;integrator=euler;events=";

#[test]
fn the_baseline_text_is_every_field_written_out() {
    assert_eq!(Scenario::baseline().to_text(), BASELINE_TEXT);
    assert_eq!(Scenario::from_text(BASELINE_TEXT), Ok(Scenario::baseline()));
}

/// A schedule holding exactly [`MAX_EVENTS`] events, one of every action.
fn full_schedule() -> Schedule {
    let mut schedule = Schedule::new();
    for slot in 0..MAX_EVENTS {
        let at = slot as f64 * 0.25;
        let fault = slot % 20 + 1;
        schedule.add(match slot % 4 {
            0 => Event::start(at, fault),
            1 => Event::stop(at, fault),
            2 => Event::new(
                at,
                Action::SetMagnitude {
                    fault,
                    magnitude: 1.0,
                },
            ),
            _ => Event::new(
                at,
                Action::Setpoint {
                    loop_index: slot % 20 + 1,
                    value: -0.5 * at,
                },
            ),
        });
    }
    schedule
}

/// The scenarios every round-trip claim is made over, each with a name so a
/// failure says which one broke.
fn edge_cases() -> Vec<(&'static str, Scenario)> {
    let mut cases = vec![
        ("baseline", Scenario::baseline()),
        ("empty schedule", Scenario::fault(4).with_hours(8.0)),
        (
            "every fault at once",
            (1..=20).fold(Scenario::baseline(), |s, n| s.with_fault(n)),
        ),
        ("open loop", Scenario::baseline().open_loop()),
        (
            "no driver IDV(12)",
            Scenario {
                driver_forces_idv12: false,
                ..Scenario::baseline()
            },
        ),
        (
            "trip ends the run",
            Scenario {
                quirks: tepsim::tepsim_core::QuirkFixes {
                    trip_ends_the_run: true,
                },
                ..Scenario::baseline()
            },
        ),
        // -0.0 and 0.0 are equal as numbers and differ in one bit. The digest
        // normalises them together on purpose; the text must still put the
        // same bits back, or `to_text` is not a canonical form.
        ("negative zero hours", Scenario::baseline().with_hours(-0.0)),
        (
            "subnormal step",
            Scenario {
                step_hours: f64::from_bits(1),
                ..Scenario::baseline()
            },
        ),
        (
            "smallest normal seed",
            Scenario::baseline().with_seed(f64::MIN_POSITIVE),
        ),
        (
            "largest finite seed",
            Scenario::baseline().with_seed(f64::MAX),
        ),
        (
            "a seed that is not an integer",
            Scenario::baseline().with_seed(4_651_207_995.25),
        ),
        (
            "one sample per step, one step",
            Scenario::baseline()
                .with_hours(1.0 / 3600.0)
                .sampling_every(1),
        ),
        (
            "a full schedule",
            Scenario {
                schedule: full_schedule(),
                extensions: {
                    // `SetMagnitude` is in the schedule above at exactly 1.0,
                    // which the original can express, so the extension is not
                    // strictly needed. It is on anyway, because the flag is a
                    // field of the scenario and this case is here to carry
                    // every field at once.
                    let mut extensions = tepsim::tepsim_core::Extensions::none();
                    extensions.continuous_disturbances = true;
                    extensions
                },
                ..Scenario::baseline()
            },
        ),
        (
            "a fractional magnitude, with the extension",
            Scenario::baseline()
                .with_continuous_disturbances()
                .with_event(Event::new(
                    3.5,
                    Action::SetMagnitude {
                        fault: 13,
                        magnitude: 0.375,
                    },
                )),
        ),
        (
            "two events at the same instant, in order",
            Scenario::baseline()
                .with_event(Event::stop(8.0, 6))
                .with_event(Event::start(8.0, 6)),
        ),
    ];
    for integrator in [
        Integrator::Euler,
        Integrator::Rk4,
        Integrator::DormandPrince,
    ] {
        cases.push((
            match integrator {
                Integrator::Euler => "euler",
                Integrator::Rk4 => "rk4",
                Integrator::DormandPrince => "dopri5",
            },
            Scenario::baseline().with_integrator(integrator),
        ));
    }
    cases
}

#[test]
fn every_scenario_round_trips_bit_for_bit() {
    for (name, scenario) in edge_cases() {
        let text = scenario.to_text();
        let back = Scenario::from_text(&text)
            .unwrap_or_else(|e| panic!("{name}: `{text}` did not parse back: {e}"));

        assert_eq!(back, scenario, "{name}: {text}");
        // `==` on `f64` puts -0.0 and 0.0 together, so the fields that can
        // carry a signed zero or a subnormal are compared by bits as well.
        assert_eq!(back.seed.to_bits(), scenario.seed.to_bits(), "{name}: seed");
        assert_eq!(
            back.hours.to_bits(),
            scenario.hours.to_bits(),
            "{name}: hours"
        );
        assert_eq!(
            back.step_hours.to_bits(),
            scenario.step_hours.to_bits(),
            "{name}: step_hours"
        );
        // Rendering is canonical: whatever came back renders identically, which
        // is the sign bit and the subnormal check for every number inside the
        // schedule as well.
        assert_eq!(back.to_text(), text, "{name}: rendering is not canonical");
    }
}

#[test]
fn the_digest_survives_a_round_trip() {
    for (name, scenario) in edge_cases() {
        let back =
            Scenario::from_text(&scenario.to_text()).unwrap_or_else(|e| panic!("{name}: {e}"));
        assert_eq!(
            back.digest(),
            scenario.digest(),
            "{name}: the digest did not survive"
        );
        assert_eq!(back.digest_hex(), scenario.digest_hex(), "{name}");
    }
}

#[test]
fn distinct_scenarios_have_distinct_texts() {
    let cases = edge_cases();
    for (i, (left_name, left)) in cases.iter().enumerate() {
        for (right_name, right) in cases.iter().skip(i + 1) {
            if left == right {
                continue;
            }
            assert_ne!(
                left.to_text(),
                right.to_text(),
                "`{left_name}` and `{right_name}` differ and render the same"
            );
        }
    }
}

/// Adding a field to `Scenario` must reach the text.
///
/// The destructuring is the teeth. It is exhaustive, so a new field stops this
/// test compiling; the loop below then makes it fail until the field actually
/// changes what `to_text` writes. Without both halves a field can be added and
/// silently fail to reach a shared link, which is the bug this module exists to
/// kill.
#[test]
fn every_field_of_a_scenario_reaches_the_text() {
    let base = Scenario::baseline();
    let Scenario {
        seed: _,
        hours: _,
        step_hours: _,
        sample_every: _,
        disturbances: _,
        controlled: _,
        quirks: _,
        driver_forces_idv12: _,
        schedule: _,
        extensions: _,
        integrator: _,
    } = base;

    let mut with_trip = base;
    with_trip.quirks.trip_ends_the_run = true;

    let changed: [(&str, Scenario); 11] = [
        ("seed", base.with_seed(1.0)),
        ("hours", base.with_hours(1.0)),
        (
            "step_hours",
            Scenario {
                step_hours: 1.0 / 7200.0,
                ..base
            },
        ),
        ("sample_every", base.sampling_every(60)),
        ("disturbances", base.with_fault(1)),
        ("controlled", base.open_loop()),
        ("quirks", with_trip),
        (
            "driver_forces_idv12",
            Scenario {
                driver_forces_idv12: false,
                ..base
            },
        ),
        ("schedule", base.with_event(Event::start(1.0, 1))),
        ("extensions", base.with_continuous_disturbances()),
        ("integrator", base.with_integrator(Integrator::Rk4)),
    ];

    let baseline_text = base.to_text();
    for (field, scenario) in changed {
        assert_ne!(
            scenario.to_text(),
            baseline_text,
            "changing `{field}` did not change the text"
        );
        assert_eq!(
            Scenario::from_text(&scenario.to_text()),
            Ok(scenario),
            "`{field}` did not come back"
        );
    }
}

// ---------------------------------------------------------------------------
// Strictness. Every one of these has to name what was wrong.
// ---------------------------------------------------------------------------

/// Parse and expect a failure, returning the error and its rendered message.
fn rejected(text: &str) -> (TextError, String) {
    match Scenario::from_text(text) {
        Ok(scenario) => panic!("`{text}` was accepted as {scenario:?}"),
        Err(error) => {
            let message = error.to_string();
            (error, message)
        }
    }
}

/// The baseline text with one `key=value` replaced.
fn baseline_with(key: &str, value: &str) -> String {
    BASELINE_TEXT
        .split(';')
        .map(|part| {
            if part.starts_with(&format!("{key}=")) {
                format!("{key}={value}")
            } else {
                part.to_string()
            }
        })
        .collect::<Vec<_>>()
        .join(";")
}

/// The baseline text with one `key=value` removed.
fn baseline_without(key: &str) -> String {
    BASELINE_TEXT
        .split(';')
        .filter(|part| !part.starts_with(&format!("{key}=")))
        .collect::<Vec<_>>()
        .join(";")
}

#[test]
fn an_unknown_field_is_named() {
    // This is the forward-compatibility case seen from the other side: a text
    // written by a build that has a field this one does not.
    let (error, message) = rejected(&format!("{BASELINE_TEXT};wobble=3"));
    assert_eq!(
        error,
        TextError::UnknownField {
            name: "wobble".to_string()
        }
    );
    assert!(message.contains("wobble"), "{message}");
    assert!(message.contains("integrator"), "{message}");
}

#[test]
fn a_later_version_is_named_not_guessed() {
    // A `v2` text names a set of fields this build does not know. Guessing
    // would mean defaulting the ones it cannot see, which is a run that is not
    // the run the text describes.
    let text = BASELINE_TEXT.replacen("v1", "v2", 1);
    let (error, message) = rejected(&text);
    assert_eq!(
        error,
        TextError::Version {
            found: "tepsim.scenario.v2".to_string()
        }
    );
    assert!(message.contains("tepsim.scenario.v1"), "{message}");
    assert!(message.contains("tepsim.scenario.v2"), "{message}");

    // And the same in the other direction: no tag at all is not "assume v1".
    let (bare, _) = rejected("seed=1;hours=1");
    assert!(matches!(bare, TextError::Version { .. }));
}

#[test]
fn every_missing_field_is_named() {
    for key in [
        "seed",
        "hours",
        "step",
        "every",
        "faults",
        "controlled",
        "idv12",
        "trip",
        "continuous",
        "integrator",
        "events",
    ] {
        let (error, message) = rejected(&baseline_without(key));
        assert_eq!(error, TextError::MissingField { name: key });
        assert!(message.contains(key), "{message}");
    }
}

#[test]
fn a_repeated_field_is_refused_rather_than_last_wins() {
    let (error, message) = rejected(&format!("{BASELINE_TEXT};hours=2"));
    assert_eq!(error, TextError::DuplicateField { name: "hours" });
    assert!(message.contains("hours"), "{message}");
}

#[test]
fn a_malformed_number_says_which_field_and_what_it_read() {
    for (key, value) in [("seed", "abc"), ("hours", "two"), ("step", "1/3600")] {
        let (error, message) = rejected(&baseline_with(key, value));
        assert_eq!(
            error,
            TextError::NotANumber {
                field: key,
                text: value.to_string()
            }
        );
        assert!(
            message.contains(key) && message.contains(value),
            "{message}"
        );
    }
    // `every` is an integer field, so a decimal is not a number for it.
    let (error, _) = rejected(&baseline_with("every", "180.0"));
    assert_eq!(
        error,
        TextError::NotANumber {
            field: "every",
            text: "180.0".to_string()
        }
    );
}

#[test]
fn out_of_range_numbers_say_what_was_expected() {
    let cases: [(&str, &str, &str); 8] = [
        ("seed", "0", "greater than zero"),
        ("seed", "-1", "greater than zero"),
        ("seed", "inf", "finite"),
        ("seed", "NaN", "finite"),
        ("hours", "-1", "not be negative"),
        ("step", "0", "greater than zero"),
        ("every", "0", "at least 1"),
        ("faults", "21", "1 to 20"),
    ];
    for (key, value, wanted) in cases {
        let (_, message) = rejected(&baseline_with(key, value));
        assert!(
            message.contains(key),
            "`{key}={value}` did not name the field: {message}"
        );
        let wanted_words: Vec<&str> = wanted.split(' ').collect();
        assert!(
            wanted_words.iter().all(|w| message.contains(w)),
            "`{key}={value}` did not say `{wanted}`: {message}"
        );
    }
}

#[test]
fn a_flag_is_zero_or_one_and_nothing_else() {
    for key in ["controlled", "idv12", "trip", "continuous"] {
        for value in ["true", "yes", "2", ""] {
            let (error, message) = rejected(&baseline_with(key, value));
            assert_eq!(
                error,
                TextError::OutOfRange {
                    field: key,
                    text: value.to_string(),
                    expected: "0 or 1",
                }
            );
            assert!(message.contains("0 or 1"), "{message}");
        }
    }
}

#[test]
fn an_unknown_integrator_is_named() {
    let (error, message) = rejected(&baseline_with("integrator", "leapfrog"));
    assert_eq!(
        error,
        TextError::UnknownIntegrator {
            name: "leapfrog".to_string()
        }
    );
    assert!(
        message.contains("leapfrog") && message.contains("euler"),
        "{message}"
    );

    // `dormand-prince` is the documented alias and must keep working; it
    // normalises to the canonical `dopri5` on the way out.
    let accepted = Scenario::from_text(&baseline_with("integrator", "dormand-prince"))
        .expect("the documented alias must parse");
    assert_eq!(accepted.integrator, Integrator::DormandPrince);
    assert!(accepted.to_text().contains("integrator=dopri5"));
}

#[test]
fn a_malformed_event_says_what_was_wrong_with_it() {
    let cases: [(&str, &str); 7] = [
        ("8", "verb"),
        ("8:start", "index"),
        ("8:jiggle:6", "verb"),
        ("-1:start:6", "negative"),
        ("8:start:21", "1 to 20"),
        ("8:magnitude:6", "value"),
        ("8:magnitude:6:2", "0 to 1"),
    ];
    for (event, wanted) in cases {
        let (error, message) = rejected(&baseline_with("events", event));
        assert!(
            matches!(error, TextError::BadEvent { .. }),
            "`{event}` gave {error:?}"
        );
        assert!(
            message.contains(event) && message.contains(wanted),
            "`{event}` did not say `{wanted}`: {message}"
        );
    }

    // Extra fields are rejected too, so `8:start:6:9` cannot be read as a
    // start with something ignored on the end.
    let (error, _) = rejected(&baseline_with("events", "8:start:6:9"));
    assert!(matches!(error, TextError::BadEvent { .. }));
}

#[test]
fn a_schedule_past_the_limit_is_refused_rather_than_truncated() {
    // `Schedule::add` panics past MAX_EVENTS, so this text cannot be produced
    // by rendering; it can very easily be produced by hand or by a hostile
    // link, and dropping the excess would give a run that is not the run the
    // text describes.
    let events: Vec<String> = (0..=MAX_EVENTS).map(|i| format!("{i}:start:1")).collect();
    let (error, message) = rejected(&baseline_with("events", &events.join(",")));
    assert_eq!(
        error,
        TextError::TooManyEvents {
            found: MAX_EVENTS + 1
        }
    );
    assert!(
        message.contains("33") && message.contains("32"),
        "{message}"
    );
}

#[test]
fn a_fractional_magnitude_without_the_extension_is_refused() {
    // The same rule `Schedule::validate` enforces, applied at the boundary:
    // `teprob.f:341-346` forces every IDV to 0 or 1, so half a fault is not
    // something the original can express, and rounding it to a whole one would
    // give a run that does not match its own description.
    let text = baseline_with("events", "3:magnitude:13:0.5");
    let (error, message) = rejected(&text);
    assert_eq!(
        error,
        TextError::NotRunnable(Invalid::ContinuousDisturbancesNotEnabled)
    );
    assert!(message.contains("continuous=1"), "{message}");

    // With the extension it is accepted, and comes back exactly.
    let allowed = baseline_with("continuous", "1");
    let allowed = allowed.replace("events=", "events=3:magnitude:13:0.5");
    let scenario = Scenario::from_text(&allowed).expect("accepted with the extension");
    assert_eq!(scenario.schedule.len(), 1);
    assert_eq!(scenario.to_text(), allowed);
}

#[test]
fn a_part_that_is_not_a_pair_is_named() {
    let (error, message) = rejected(&format!("{BASELINE_TEXT};nonsense"));
    assert_eq!(
        error,
        TextError::NotAPair {
            part: "nonsense".to_string()
        }
    );
    assert!(message.contains("nonsense"), "{message}");
}

#[test]
fn an_empty_text_is_refused() {
    let (error, message) = rejected("");
    assert_eq!(error, TextError::Empty);
    assert!(message.contains("tepsim.scenario.v1"), "{message}");
    assert_eq!(rejected("   ").0, TextError::Empty);
}

// ---------------------------------------------------------------------------
// Shape of the format.
// ---------------------------------------------------------------------------

#[test]
fn field_order_does_not_matter_but_the_rendered_order_is_canonical() {
    let mut parts: Vec<&str> = BASELINE_TEXT.split(';').collect();
    let tag = parts.remove(0);
    parts.reverse();
    let reordered = format!("{tag};{}", parts.join(";"));
    assert_ne!(reordered, BASELINE_TEXT);

    let scenario = Scenario::from_text(&reordered).expect("any order parses");
    assert_eq!(scenario, Scenario::baseline());
    // Rendering puts it back in the one order, so two texts for one scenario
    // become one text again.
    assert_eq!(scenario.to_text(), BASELINE_TEXT);
}

#[test]
fn the_text_needs_no_percent_encoding_in_a_url_fragment() {
    // `apps/studio/js/share.js` puts this straight into a fragment. Anything
    // outside the fragment-safe set would come back percent-encoded and the
    // link would stop being readable, which is half of what a link is for.
    for (name, scenario) in edge_cases() {
        for byte in scenario.to_text().bytes() {
            let safe = byte.is_ascii_alphanumeric() || b"-._~;=,:+".contains(&byte);
            assert!(
                safe,
                "{name}: byte {byte:#x} (`{}`) is not fragment safe",
                byte as char
            );
        }
    }
}

#[test]
fn a_subnormal_is_short_and_a_common_number_is_readable() {
    // The number rule is "whichever of `{}` and `{:e}` is shorter". Both are
    // exact; the rule is about whether a human can read the result.
    assert!(Scenario::baseline().to_text().contains("hours=48"));
    assert!(Scenario::baseline().to_text().contains("every=180"));
    // One second in hours. `{}` gives 0.0002777777777777778 and `{:e}` gives
    // 2.777777777777778e-4, which is two characters shorter and therefore the
    // one written. Both parse back to the same bits.
    assert!(
        Scenario::baseline()
            .to_text()
            .contains("step=2.777777777777778e-4")
    );

    let subnormal = Scenario {
        step_hours: f64::from_bits(1),
        ..Scenario::baseline()
    };
    assert!(
        subnormal.to_text().contains("step=5e-324"),
        "{}",
        subnormal.to_text()
    );
}

#[test]
fn write_text_and_to_text_agree() {
    // `write_text` is the allocation-free half, for a caller writing into a
    // file or a header rather than into a `String`.
    use core::fmt::Write;
    for (name, scenario) in edge_cases() {
        let mut out = String::new();
        write!(out, "{}", Rendered(&scenario)).expect("writing into a String cannot fail");
        assert_eq!(out, scenario.to_text(), "{name}");
    }

    struct Rendered<'a>(&'a Scenario);
    impl core::fmt::Display for Rendered<'_> {
        fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
            text::write_text(f, self.0)
        }
    }
}
