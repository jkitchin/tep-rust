//! A canonical text form for a [`Scenario`], and a strict parser for it.
//!
//! # Why this exists
//!
//! A run is a pure function of its scenario, so a scenario is the thing worth
//! moving: into a URL fragment, into a dataset header, into a bug report, over
//! a `postMessage` boundary. Before this module every one of those places
//! enumerated the fields by hand, and the browser did it three times over. A
//! field added to [`Scenario`] reached none of them, and the failure was
//! silent: the link still opened, the run still ran, and it was a different run
//! from the one the link claimed.
//!
//! One serialisation, used by every boundary, turns that into a compile error
//! or a test failure. `every_field_of_a_scenario_reaches_the_text` in
//! `tests/scenario_text.rs` destructures [`Scenario`] exhaustively, so adding a
//! field there stops the test compiling until the field is handled here.
//!
//! # Why not serde
//!
//! Serde would work and was rejected on three counts. The surface is eleven
//! scalars and up to [`MAX_EVENTS`] events, so the derive saves perhaps eighty
//! lines and costs a proc-macro dependency in a crate that is `no_std` and
//! ships to a browser with a size budget. The strictness this needs (every
//! field present, unknown fields named in the error, ranges checked at the
//! boundary rather than after it) is `deny_unknown_fields` plus a hand-written
//! `Deserialize` in practice, which is most of the work anyway. And the format
//! has to survive in a URL fragment, where JSON's braces and quotes are
//! percent-encoded into noise; this one is fragment-safe by construction.
//!
//! # The format
//!
//! Semicolon-separated. The first part is the version tag, and the rest are
//! `key=value` pairs. Every key must appear exactly once, in any order, and the
//! renderer always emits them in the order of the table below.
//!
//! ```text
//! tepsim.scenario.v1;seed=4651207995;hours=48;step=2.777777777777778e-4;\
//! every=180;faults=;controlled=1;idv12=1;trip=0;continuous=0;\
//! integrator=euler;events=
//! ```
//!
//! | Key | Field | Value |
//! |---|---|---|
//! | `seed` | [`Scenario::seed`] | a finite number greater than zero |
//! | `hours` | [`Scenario::hours`] | a finite number, not negative |
//! | `step` | `step_hours` | a finite number greater than zero |
//! | `every` | `sample_every` | an integer, at least 1 |
//! | `faults` | `disturbances` | one-based indices, comma separated, ascending |
//! | `controlled` | `controlled` | `0` or `1` |
//! | `idv12` | `driver_forces_idv12` | `0` or `1` |
//! | `trip` | `quirks.trip_ends_the_run` | `0` or `1` |
//! | `continuous` | `extensions.continuous_disturbances` | `0` or `1` |
//! | `integrator` | `integrator` | `euler`, `rk4` or `dopri5` |
//! | `events` | `schedule` | comma-separated events, see below |
//!
//! An event is `time:verb` followed by the verb's own fields:
//!
//! ```text
//! 8:start:6            IDV(6) on at hour 8
//! 20:stop:6            IDV(6) off at hour 20
//! 4:magnitude:13:0.5   IDV(13) at half strength at hour 4
//! 10:setpoint:9:0.25   loop 9 moved to 0.25 at hour 10
//! ```
//!
//! Only the characters `A-Z a-z 0-9 - . _ ~ ; = , :` and `+` appear, all of
//! which a URL fragment carries without percent-encoding, so
//! `apps/studio/js/share.js` can put a whole scenario in a link verbatim.
//!
//! # Numbers
//!
//! Every number is rendered by whichever of Rust's `{}` and `{:e}` is shorter,
//! with `{}` winning a tie. Both produce the shortest decimal that parses back
//! to the same `f64`, so the choice is between two exact renderings and costs
//! nothing: `48` and `180` stay as they are, while a subnormal comes out as
//! `5e-324` rather than as three hundred characters of zeros. Negative zero
//! renders as `-0` and parses back to negative zero, bit for bit.
//!
//! # What round-trips
//!
//! [`Scenario::from_text`] applied to [`Scenario::to_text`] returns an equal
//! scenario, with every `f64` equal *by bit pattern*, for every scenario whose
//! numbers are in the ranges the table above states. That is the whole property
//! and `tests/scenario_text.rs` asserts it over the edge cases: an empty
//! schedule, a full one, every integrator, negative zero and subnormals.
//!
//! It follows that [`Scenario::digest`] survives a round trip, which is the
//! property that makes a serialised scenario worth anything: a dataset labelled
//! with a digest and shipped with its scenario text can be checked, rather than
//! believed.
//!
//! A scenario holding a value outside those ranges renders anyway, and is then
//! *rejected* on the way back in, by name and by reason. That asymmetry is
//! deliberate. A `NaN` duration is not equal to itself, so no round trip could
//! ever hold for it; accepting it here would move the failure to whichever
//! machine ran the file rather than stopping it at the boundary.
//!
//! # Forward compatibility
//!
//! The version tag is [`SCENARIO_VERSION`], the same string
//! [`Scenario::digest`] absorbs, so the two cannot drift: a change to what a
//! scenario contains changes both at once.
//!
//! Adding a field means adding a key and bumping the tag to `v2`. A `v1` text
//! is then rejected by a `v2` build with a message naming both versions, and a
//! `v2` text is rejected by a `v1` build the same way. Nothing is defaulted
//! silently in either direction, which is the point: a link written before a
//! field existed describes a run that no longer has a unique meaning, and
//! saying so is better than picking one. `a_later_version_is_named_not_guessed`
//! and `an_unknown_field_is_named` in `tests/scenario_text.rs` are that
//! behaviour, today, on this build.

use alloc::string::{String, ToString};
use core::fmt::{self, Write};

use tepsim_core::{Extensions, QuirkFixes};
use tepsim_scenario::{Action, Event, Invalid, MAX_EVENTS, Schedule};

use crate::integrator::Integrator;
use crate::scenario::{DISTURBANCES, SCENARIO_VERSION, Scenario};

/// The keys, in the order the renderer emits them.
///
/// Named once so the renderer, the parser and the "did you miss one" check
/// cannot disagree about the set.
const KEYS: [&str; 11] = [
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
];

/// Why a scenario text could not be parsed.
///
/// A closed set rather than a string, so a caller can match on the cause; every
/// variant carries enough to say what was wrong and where. [`fmt::Display`]
/// turns each into the sentence the bindings hand to Python and JavaScript.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum TextError {
    /// The text was empty.
    Empty,
    /// The leading version tag was not [`SCENARIO_VERSION`].
    Version {
        /// What the text actually began with.
        found: String,
    },
    /// A part had no `=` in it.
    NotAPair {
        /// The part as written.
        part: String,
    },
    /// A key that this version knows nothing about.
    UnknownField {
        /// The key as written.
        name: String,
    },
    /// A key that this version requires and the text does not have.
    MissingField {
        /// The key that should have been there.
        name: &'static str,
    },
    /// A key given more than once.
    ///
    /// Rejected rather than last-wins, because two different values for one
    /// field is a text with two meanings and no way to tell which was meant.
    DuplicateField {
        /// The key that repeated.
        name: &'static str,
    },
    /// A value that is not a number at all.
    NotANumber {
        /// Which key.
        field: &'static str,
        /// The value as written.
        text: String,
    },
    /// A number outside what the field admits.
    OutOfRange {
        /// Which key.
        field: &'static str,
        /// The value as written.
        text: String,
        /// What would have been accepted.
        expected: &'static str,
    },
    /// An integrator name this build does not have.
    UnknownIntegrator {
        /// The name as written.
        name: String,
    },
    /// An event that is not `time:verb` with the verb's own fields.
    BadEvent {
        /// The event as written.
        text: String,
        /// What was wrong with it.
        why: &'static str,
    },
    /// More than [`MAX_EVENTS`] events.
    TooManyEvents {
        /// How many the text asked for.
        found: usize,
    },
    /// The scenario parsed, and then failed [`Scenario::validate`].
    NotRunnable(Invalid),
}

impl fmt::Display for TextError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => write!(
                f,
                "a scenario text cannot be empty: it must begin `{SCENARIO_VERSION};`"
            ),
            Self::Version { found } => write!(
                f,
                "this build reads `{SCENARIO_VERSION}` and the text begins `{found}`. \
                 A scenario text is versioned with what a scenario contains, so a \
                 different version describes a run this build cannot reconstruct"
            ),
            Self::NotAPair { part } => {
                write!(f, "`{part}` is not a `key=value` pair")
            }
            Self::UnknownField { name } => write!(
                f,
                "unknown field `{name}`: this version has {}",
                Listed(&KEYS)
            ),
            Self::MissingField { name } => write!(
                f,
                "missing field `{name}`: every field must be written out, so that a \
                 text says what it runs rather than what it leaves to a default"
            ),
            Self::DuplicateField { name } => {
                write!(f, "field `{name}` appears more than once")
            }
            Self::NotANumber { field, text } => {
                write!(f, "field `{field}`: `{text}` is not a number")
            }
            Self::OutOfRange {
                field,
                text,
                expected,
            } => write!(
                f,
                "field `{field}`: `{text}` is out of range, expected {expected}"
            ),
            Self::UnknownIntegrator { name } => write!(
                f,
                "unknown integrator `{name}`: expected euler, rk4 or dopri5"
            ),
            Self::BadEvent { text, why } => write!(f, "event `{text}`: {why}"),
            Self::TooManyEvents { found } => write!(
                f,
                "{found} events, and a schedule holds at most {MAX_EVENTS}"
            ),
            Self::NotRunnable(invalid) => {
                write!(
                    f,
                    "the scenario parsed but cannot run: {}",
                    Reason(*invalid)
                )
            }
        }
    }
}

/// A comma-separated list of names, for an error message.
struct Listed<'a>(&'a [&'a str]);

impl fmt::Display for Listed<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for (index, name) in self.0.iter().enumerate() {
            if index > 0 {
                f.write_str(", ")?;
            }
            f.write_str(name)?;
        }
        Ok(())
    }
}

/// [`Invalid`] as a sentence.
///
/// `Invalid` lives in `tepsim-scenario` and carries no `Display`; writing one
/// there would be a change to a crate this module has no business editing, and
/// the wording wanted here is the wording of a parse error rather than of a
/// validation failure.
struct Reason(Invalid);

impl fmt::Display for Reason {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.0 {
            Invalid::FaultOutOfRange { fault } => {
                write!(
                    f,
                    "IDV({fault}) does not exist; this model has {DISTURBANCES}"
                )
            }
            Invalid::LoopOutOfRange { loop_index } => {
                write!(f, "control loop {loop_index} does not exist; there are 20")
            }
            Invalid::MagnitudeOutOfRange => f.write_str("a magnitude outside 0 to 1"),
            Invalid::ContinuousDisturbancesNotEnabled => f.write_str(
                "a fractional magnitude without `continuous=1`, which the original \
                 cannot express: `teprob.f:341-346` forces every IDV to 0 or 1",
            ),
            Invalid::TimeNotFinite => f.write_str("an event at a negative or non-finite time"),
        }
    }
}

// ---------------------------------------------------------------------------
// Rendering.
// ---------------------------------------------------------------------------

/// Counts bytes instead of keeping them, so the two candidate renderings of a
/// number can be compared without allocating.
struct Counter(usize);

impl Write for Counter {
    fn write_str(&mut self, s: &str) -> fmt::Result {
        self.0 += s.len();
        Ok(())
    }
}

/// Write `value` in whichever of `{}` and `{:e}` is shorter, `{}` winning ties.
///
/// Both are the shortest decimal that parses back to the same `f64`, so this
/// picks between two exact renderings on length alone and can lose no accuracy.
/// Plain notation keeps `48` and `180` as a reader expects them; exponent
/// notation keeps a subnormal from becoming three hundred characters of zeros,
/// which `{}` alone would produce because Rust's `Display` never uses an
/// exponent.
fn write_number<W: Write>(out: &mut W, value: f64) -> fmt::Result {
    let mut plain = Counter(0);
    write!(plain, "{value}")?;
    let mut exponent = Counter(0);
    write!(exponent, "{value:e}")?;
    if exponent.0 < plain.0 {
        write!(out, "{value:e}")
    } else {
        write!(out, "{value}")
    }
}

/// Write the scenario's canonical text into `out`.
///
/// The allocation-free half of [`Scenario::to_text`], for a caller writing
/// straight into a file or a header.
///
/// # Errors
///
/// Whatever `out` returns. Nothing here can fail on its own.
pub fn write_text<W: Write>(out: &mut W, scenario: &Scenario) -> fmt::Result {
    out.write_str(SCENARIO_VERSION)?;

    out.write_str(";seed=")?;
    write_number(out, scenario.seed)?;
    out.write_str(";hours=")?;
    write_number(out, scenario.hours)?;
    out.write_str(";step=")?;
    write_number(out, scenario.step_hours)?;
    write!(out, ";every={}", scenario.sample_every)?;

    out.write_str(";faults=")?;
    for (position, fault) in scenario.active_faults().enumerate() {
        if position > 0 {
            out.write_str(",")?;
        }
        write!(out, "{fault}")?;
    }

    write!(out, ";controlled={}", digit(scenario.controlled))?;
    write!(out, ";idv12={}", digit(scenario.driver_forces_idv12))?;
    write!(out, ";trip={}", digit(scenario.quirks.trip_ends_the_run))?;
    write!(
        out,
        ";continuous={}",
        digit(scenario.extensions.continuous_disturbances)
    )?;
    write!(out, ";integrator={}", scenario.integrator.name())?;

    out.write_str(";events=")?;
    for (position, event) in scenario.schedule.events().enumerate() {
        if position > 0 {
            out.write_str(",")?;
        }
        write_number(out, event.at_hours)?;
        match event.action {
            Action::Start { fault } => write!(out, ":start:{fault}")?,
            Action::Stop { fault } => write!(out, ":stop:{fault}")?,
            Action::SetMagnitude { fault, magnitude } => {
                write!(out, ":magnitude:{fault}:")?;
                write_number(out, magnitude)?;
            }
            Action::Setpoint { loop_index, value } => {
                write!(out, ":setpoint:{loop_index}:")?;
                write_number(out, value)?;
            }
        }
    }

    Ok(())
}

/// `1` or `0`, which is what the format writes a flag as.
const fn digit(flag: bool) -> u8 {
    if flag { 1 } else { 0 }
}

/// The scenario's canonical text.
#[must_use]
pub fn to_text(scenario: &Scenario) -> String {
    let mut out = String::new();
    // Writing into a `String` is infallible, so the result cannot be an error
    // and there is nothing to propagate.
    let _ = write_text(&mut out, scenario);
    out
}

// ---------------------------------------------------------------------------
// Parsing.
// ---------------------------------------------------------------------------

/// Parse a canonical scenario text.
///
/// # Errors
///
/// The first problem found, reading left to right, then the completeness check,
/// then [`Scenario::validate`]. See [`TextError`].
pub fn from_text(text: &str) -> Result<Scenario, TextError> {
    let text = text.trim();
    if text.is_empty() {
        return Err(TextError::Empty);
    }

    let mut parts = text.split(';');
    // `split` on a non-empty string always yields at least one item.
    let tag = parts.next().unwrap_or("");
    if tag != SCENARIO_VERSION {
        return Err(TextError::Version {
            found: tag.to_string(),
        });
    }

    // Slots rather than a map: eleven keys, `no_std`, and this way the
    // completeness check below is a loop over `KEYS` with no allocation.
    let mut values: [Option<&str>; KEYS.len()] = [None; KEYS.len()];
    for part in parts {
        let Some(split) = part.find('=') else {
            return Err(TextError::NotAPair {
                part: part.to_string(),
            });
        };
        let (name, value) = (&part[..split], &part[split + 1..]);
        let Some(index) = KEYS.iter().position(|key| *key == name) else {
            return Err(TextError::UnknownField {
                name: name.to_string(),
            });
        };
        if values[index].is_some() {
            return Err(TextError::DuplicateField { name: KEYS[index] });
        }
        values[index] = Some(value);
    }

    // Every field, always. A text that leaves one out is a text whose meaning
    // depends on the reader's defaults, which is exactly what the digest exists
    // to rule out.
    let mut taken: [&str; KEYS.len()] = [""; KEYS.len()];
    for (index, slot) in values.iter().enumerate() {
        match slot {
            Some(value) => taken[index] = value,
            None => return Err(TextError::MissingField { name: KEYS[index] }),
        }
    }
    // Destructured rather than looked up by name, so the bindings below are
    // pinned to the order of `KEYS`: adding a key there stops this compiling
    // until it is handled here.
    let [
        seed,
        hours,
        step,
        every,
        fault_list,
        controlled,
        idv12,
        trip,
        continuous,
        integrator,
        events,
    ] = taken;

    let mut scenario = Scenario::baseline();
    scenario.seed = positive("seed", seed)?;
    scenario.hours = not_negative("hours", hours)?;
    scenario.step_hours = positive("step", step)?;
    scenario.sample_every = at_least_one("every", every)?;
    scenario.disturbances = faults(fault_list)?;
    scenario.controlled = flag("controlled", controlled)?;
    scenario.driver_forces_idv12 = flag("idv12", idv12)?;
    scenario.quirks = QuirkFixes {
        trip_ends_the_run: flag("trip", trip)?,
    };
    // Built by mutation rather than as a literal: `Extensions` is
    // `#[non_exhaustive]`, so a struct expression is not allowed outside its
    // own crate. `QuirkFixes` is not, and is written as a literal on purpose,
    // so that adding a quirk fails to compile here until it is serialised.
    let mut extensions = Extensions::none();
    extensions.continuous_disturbances = flag("continuous", continuous)?;
    scenario.extensions = extensions;
    scenario.integrator =
        Integrator::parse(integrator).ok_or_else(|| TextError::UnknownIntegrator {
            name: integrator.to_string(),
        })?;
    scenario.schedule = schedule(events)?;

    scenario.validate().map_err(TextError::NotRunnable)?;
    Ok(scenario)
}

/// A finite number.
fn number(name: &'static str, text: &str) -> Result<f64, TextError> {
    let value: f64 = text.parse().map_err(|_| TextError::NotANumber {
        field: name,
        text: text.to_string(),
    })?;
    if value.is_finite() {
        Ok(value)
    } else {
        Err(TextError::OutOfRange {
            field: name,
            text: text.to_string(),
            expected: "a finite number",
        })
    }
}

/// A finite number greater than zero.
fn positive(name: &'static str, text: &str) -> Result<f64, TextError> {
    let value = number(name, text)?;
    if value > 0.0 {
        Ok(value)
    } else {
        Err(TextError::OutOfRange {
            field: name,
            text: text.to_string(),
            expected: "a finite number greater than zero",
        })
    }
}

/// A finite number that is not negative. Negative zero is admitted and kept.
fn not_negative(name: &'static str, text: &str) -> Result<f64, TextError> {
    let value = number(name, text)?;
    if value >= 0.0 {
        Ok(value)
    } else {
        Err(TextError::OutOfRange {
            field: name,
            text: text.to_string(),
            expected: "a finite number that is not negative",
        })
    }
}

/// A positive integer.
fn at_least_one(name: &'static str, text: &str) -> Result<usize, TextError> {
    let value: usize = text.parse().map_err(|_| TextError::NotANumber {
        field: name,
        text: text.to_string(),
    })?;
    if value >= 1 {
        Ok(value)
    } else {
        Err(TextError::OutOfRange {
            field: name,
            text: text.to_string(),
            expected: "an integer of at least 1",
        })
    }
}

/// `0` or `1`, and nothing else.
fn flag(name: &'static str, text: &str) -> Result<bool, TextError> {
    match text {
        "0" => Ok(false),
        "1" => Ok(true),
        _ => Err(TextError::OutOfRange {
            field: name,
            text: text.to_string(),
            expected: "0 or 1",
        }),
    }
}

/// A one-based `IDV` index.
fn fault_index(name: &'static str, text: &str) -> Result<usize, TextError> {
    let value: usize = text.parse().map_err(|_| TextError::NotANumber {
        field: name,
        text: text.to_string(),
    })?;
    if (1..=DISTURBANCES).contains(&value) {
        Ok(value)
    } else {
        Err(TextError::OutOfRange {
            field: name,
            text: text.to_string(),
            expected: "a one-based IDV index in 1 to 20",
        })
    }
}

/// The comma-separated fault list.
fn faults(text: &str) -> Result<[bool; DISTURBANCES], TextError> {
    let mut active = [false; DISTURBANCES];
    if text.is_empty() {
        return Ok(active);
    }
    for part in text.split(',') {
        active[fault_index("faults", part)? - 1] = true;
    }
    Ok(active)
}

/// The comma-separated event list.
fn schedule(text: &str) -> Result<Schedule, TextError> {
    let mut built = Schedule::new();
    if text.is_empty() {
        return Ok(built);
    }
    let count = text.split(',').count();
    if count > MAX_EVENTS {
        return Err(TextError::TooManyEvents { found: count });
    }
    for part in text.split(',') {
        built.add(event(part)?);
    }
    Ok(built)
}

/// One `time:verb[:index[:value]]` event.
fn event(text: &str) -> Result<Event, TextError> {
    let bad = |why: &'static str| TextError::BadEvent {
        text: text.to_string(),
        why,
    };

    let mut fields = text.split(':');
    let at = fields.next().ok_or_else(|| bad("no time"))?;
    let at_hours = match at.parse::<f64>() {
        Ok(value) if value.is_finite() && value >= 0.0 => value,
        _ => {
            return Err(bad(
                "the time must be a finite number of hours, not negative",
            ));
        }
    };
    let verb = fields.next().ok_or_else(|| {
        bad("expected `time:verb`, where verb is start, stop, magnitude or setpoint")
    })?;
    let index = fields
        .next()
        .ok_or_else(|| bad("expected an index after the verb"))?;

    let action = match verb {
        "start" | "stop" => {
            let fault =
                fault_index("events", index).map_err(|_| bad("IDV index must be 1 to 20"))?;
            if verb == "start" {
                Action::Start { fault }
            } else {
                Action::Stop { fault }
            }
        }
        "magnitude" => {
            let fault =
                fault_index("events", index).map_err(|_| bad("IDV index must be 1 to 20"))?;
            let value = fields
                .next()
                .ok_or_else(|| bad("`magnitude` needs a value after the IDV index"))?;
            let magnitude = match value.parse::<f64>() {
                Ok(m) if m.is_finite() && (0.0..=1.0).contains(&m) => m,
                _ => return Err(bad("a magnitude must be a finite number from 0 to 1")),
            };
            Action::SetMagnitude { fault, magnitude }
        }
        "setpoint" => {
            let loop_index: usize = index
                .parse()
                .ok()
                .filter(|n| (1..=20).contains(n))
                .ok_or_else(|| bad("a control loop index must be 1 to 20"))?;
            let value = fields
                .next()
                .ok_or_else(|| bad("`setpoint` needs a value after the loop index"))?;
            let value = match value.parse::<f64>() {
                Ok(v) if v.is_finite() => v,
                _ => return Err(bad("a setpoint must be a finite number")),
            };
            Action::Setpoint { loop_index, value }
        }
        _ => return Err(bad("the verb must be start, stop, magnitude or setpoint")),
    };

    if fields.next().is_some() {
        return Err(bad("too many fields for this verb"));
    }
    Ok(Event::new(at_hours, action))
}
