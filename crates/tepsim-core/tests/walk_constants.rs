//! Re-derives the noise magnitudes and walk spans from the vendored Fortran
//! and compares them against the transcribed tables.
//!
//! The same treatment as `tests/constants.rs` and `tests/nominal_state.rs`,
//! and for the same reason: 101 numbers is well past what anybody checks by
//! eye. A wrong noise magnitude would show up only in Tier 5, as a
//! distribution that is subtly the wrong width, and a wrong span only after a
//! channel had wandered for hours.
//!
//! An integration test because `tepsim-core` is `no_std` and cannot read
//! files.

use std::collections::BTreeMap;
use std::path::PathBuf;

use tepsim_core::CHANNEL_SPANS;
use tepsim_core::constants::MEASUREMENT_NOISE;

fn source() -> String {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../reference/fortran/teprob.f");
    std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("reading {}: {e}", path.display()))
}

/// One `NAME(i)=literal` assignment, as written.
#[derive(Debug)]
struct Parsed {
    literal: String,
    is_double: bool,
    value: f64,
}

/// Parse every assignment to one of the six arrays, keyed by name and index.
fn parse(text: &str) -> BTreeMap<(&'static str, usize), Parsed> {
    const ARRAYS: [&str; 6] = ["XNS", "HSPAN", "HZERO", "SSPAN", "SZERO", "SPSPAN"];
    let mut out = BTreeMap::new();
    for raw in text.lines() {
        let Some(rest) = raw.strip_prefix("      ") else {
            continue;
        };
        let Some((lhs, literal)) = rest.split_once('=') else {
            continue;
        };
        let Some((name, index)) = lhs.split_once('(') else {
            continue;
        };
        let Some(array) = ARRAYS.into_iter().find(|a| *a == name) else {
            continue;
        };
        let Some(index) = index.strip_suffix(')') else {
            continue;
        };
        let Ok(index) = index.parse::<usize>() else {
            continue;
        };
        let literal = literal.trim().to_string();
        let is_double = literal.contains('D');
        // A leading or trailing dot is legal in Fortran and not in Rust's
        // parser, so normalise before parsing rather than trusting either.
        let mut mantissa = literal.replace('D', "E");
        if let Some((m, e)) = mantissa.split_once('E') {
            mantissa = format!("{}E{}", pad(m), e);
        } else {
            mantissa = pad(&mantissa);
        }
        let value = if is_double {
            mantissa.parse::<f64>().expect("a double literal")
        } else {
            f64::from(mantissa.parse::<f32>().expect("a single literal"))
        };
        out.insert(
            (array, index),
            Parsed {
                literal,
                is_double,
                value,
            },
        );
    }
    out
}

/// `.005` and `1.` are Fortran, not Rust.
fn pad(mantissa: &str) -> String {
    let mut s = mantissa.to_string();
    if s.starts_with('.') {
        s.insert(0, '0');
    } else if s.starts_with("-.") {
        s.insert(1, '0');
    }
    if s.ends_with('.') {
        s.push('0');
    }
    s
}

#[test]
fn every_noise_magnitude_matches_the_fortran() {
    let text = source();
    let parsed = parse(&text);
    for index in 1..=41 {
        let entry = parsed
            .get(&("XNS", index))
            .unwrap_or_else(|| panic!("XNS({index}) is not in the Fortran"));
        assert_eq!(
            MEASUREMENT_NOISE[index - 1].to_bits(),
            entry.value.to_bits(),
            "XNS({index}) = {}: transcribed as {}, the Fortran stores {}",
            entry.literal,
            MEASUREMENT_NOISE[index - 1],
            entry.value
        );
    }
}

#[test]
fn every_channel_span_matches_the_fortran() {
    let text = source();
    let parsed = parse(&text);
    for channel in 1..=12 {
        let spans = &CHANNEL_SPANS[channel - 1];
        let cases = [
            ("HSPAN", spans.duration_span),
            ("HZERO", spans.duration_centre),
            ("SSPAN", spans.value_span),
            ("SZERO", spans.value_centre),
            ("SPSPAN", spans.slope_span),
        ];
        for (array, ours) in cases {
            let entry = parsed
                .get(&(array, channel))
                .unwrap_or_else(|| panic!("{array}({channel}) is not in the Fortran"));
            assert_eq!(
                ours.to_bits(),
                entry.value.to_bits(),
                "{array}({channel}) = {}: transcribed as {ours}, the Fortran \
                 stores {}",
                entry.literal,
                entry.value
            );
        }
    }
}

/// All 101 carry a `D` suffix, without exception.
///
/// Checked across every one rather than inferred from a sample. The `273.15`
/// and `1.8` findings both came from the original being inconsistent about
/// exactly this, in blocks that looked uniform.
#[test]
fn all_one_hundred_and_one_are_double_precision() {
    let text = source();
    let parsed = parse(&text);
    let single: Vec<String> = parsed
        .iter()
        .filter(|(_, v)| !v.is_double)
        .map(|((a, i), v)| format!("{a}({i}) = {}", v.literal))
        .collect();
    assert_eq!(parsed.len(), 41 + 5 * 12, "the block changed size");
    assert!(
        single.is_empty(),
        "{} of the 101 are single precision after all: {single:?}",
        single.len()
    );
}

/// `SPSPAN` is zero for every channel, so no segment ever ends with a slope.
///
/// If this ever stops being true, `walk_segment`'s third draw starts affecting
/// the answer rather than only the stream, and the module documentation is
/// wrong.
#[test]
fn no_channel_has_a_nonzero_slope_span() {
    for (channel, spans) in CHANNEL_SPANS.iter().enumerate() {
        assert_eq!(
            spans.slope_span.to_bits(),
            0.0_f64.to_bits(),
            "channel {} has a slope span",
            channel + 1
        );
    }
}

/// The three spike channels have no value or slope span, because the spike
/// rule at `teprob.f:372-396` never reads them.
#[test]
fn the_three_spike_channels_have_only_duration_parameters() {
    for channel in 10..=12 {
        let spans = &CHANNEL_SPANS[channel - 1];
        assert_eq!(
            spans.value_span.to_bits(),
            0.0_f64.to_bits(),
            "SSPAN({channel})"
        );
        assert_eq!(
            spans.value_centre.to_bits(),
            0.0_f64.to_bits(),
            "SZERO({channel})"
        );
        assert!(
            spans.duration_span > 0.0 && spans.duration_centre > 0.0,
            "channel {channel} has no duration either, so it never fires"
        );
    }
    // And the nine walk channels do have value spans, or they would be
    // constants and the disturbances they carry would do nothing.
    for channel in 1..=9 {
        assert!(
            CHANNEL_SPANS[channel - 1].value_span > 0.0,
            "walk channel {channel} cannot move"
        );
    }
}

/// Every duration centre clears its span, so a segment can never have zero or
/// negative length.
///
/// `H` divides the cubic coefficients at `teprob.f:1533-1534` and the original
/// has no guard, so this is what stops the walk producing infinities.
#[test]
fn no_channel_can_draw_a_zero_length_segment() {
    for (channel, spans) in CHANNEL_SPANS.iter().enumerate() {
        assert!(
            spans.duration_centre > spans.duration_span,
            "channel {}: centre {} does not clear span {}, so a draw of -1 \
             gives a segment of length {} and the cubic divides by it",
            channel + 1,
            spans.duration_centre,
            spans.duration_span,
            spans.duration_centre - spans.duration_span
        );
    }
}
