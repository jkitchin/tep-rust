//! Re-derives the nominal initial state from the vendored Fortran and
//! compares it against the transcribed constant.
//!
//! The same approach as `tests/constants.rs`, and for the same reason: fifty
//! numbers is well past what anybody checks reliably by eye, and a wrong one
//! would move the *starting point* rather than the model, so every Tier 1 and
//! Tier 2 number would stay perfect while the trajectory drifted from the
//! published data.
//!
//! An integration test because `tepsim-core` is `no_std` and cannot read
//! files.

use std::path::PathBuf;

use tepsim_core::constants::NOMINAL_STATE;

fn source() -> String {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../reference/fortran/teprob.f");
    std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("reading {}: {e}", path.display()))
}

/// One `YY(i)=literal` line, as written.
struct Parsed {
    index: usize,
    literal: String,
    is_double: bool,
    value: f64,
}

fn parse(text: &str) -> Vec<Parsed> {
    let mut out = Vec::new();
    for raw in text.lines() {
        let Some(rest) = raw.strip_prefix("      YY(") else {
            continue;
        };
        let Some((index, literal)) = rest.split_once(")=") else {
            continue;
        };
        let Ok(index) = index.parse::<usize>() else {
            continue;
        };
        let literal = literal.trim().to_string();
        // Fixed-form Fortran: a `D` exponent means double, anything else is
        // single and is widened on assignment.
        let is_double = literal.contains('D');
        let normalised = literal.replace('D', "E");
        let value = if is_double {
            normalised.parse::<f64>().expect("a double literal")
        } else {
            f64::from(normalised.parse::<f32>().expect("a single literal"))
        };
        out.push(Parsed {
            index,
            literal,
            is_double,
            value,
        });
    }
    out
}

#[test]
fn every_nominal_state_matches_the_fortran_literal_and_its_precision() {
    let text = source();
    let parsed = parse(&text);
    // `TEINIT` sets all fifty, and nothing else in the file writes `YY(i)=`
    // with a literal on the right.
    assert_eq!(parsed.len(), 50, "expected fifty YY assignments");

    for entry in &parsed {
        let ours = NOMINAL_STATE[entry.index - 1];
        assert_eq!(
            ours.to_bits(),
            entry.value.to_bits(),
            "YY({}) = {} ({}): transcribed as {ours:e}, the Fortran stores \
             {:e}. Check the `D` suffix.",
            entry.index,
            entry.literal,
            if entry.is_double { "double" } else { "single" },
            entry.value,
        );
    }
}

/// Three of the fifty are double precision, and they are the three written in
/// exponential form. If that ratio ever changes, the file changed.
#[test]
fn exactly_three_of_the_fifty_are_double_precision() {
    let text = source();
    let parsed = parse(&text);
    let doubles: Vec<usize> = parsed
        .iter()
        .filter(|p| p.is_double)
        .map(|p| p.index)
        .collect();
    assert_eq!(
        doubles,
        vec![20, 22, 24],
        "the double-precision entries moved"
    );
}

/// The single-precision ones really do lose digits. If they did not, the
/// distinction would be untestable and `single` would be ceremony here.
#[test]
fn the_single_precision_entries_are_not_the_nearest_double() {
    let text = source();
    let parsed = parse(&text);
    let lossy = parsed
        .iter()
        .filter(|p| !p.is_double)
        .filter(|p| {
            let exact: f64 = p.literal.parse().expect("a decimal literal");
            exact.to_bits() != p.value.to_bits()
        })
        .count();
    assert!(
        lossy > 40,
        "only {lossy} of the 47 single-precision entries lose anything by \
         being stored as f32, so the precision distinction is nearly free \
         here and the reasoning should be rechecked"
    );
}

/// The valve positions are the last twelve states, and one of them is exactly
/// 50. A transposed block would be invisible in a sum.
#[test]
fn the_last_twelve_entries_are_the_valve_positions() {
    // YY(50) = 50.00000000, the agitator valve, exactly representable.
    assert_eq!(NOMINAL_STATE[49].to_bits(), 50.0_f64.to_bits());
    for position in &NOMINAL_STATE[38..50] {
        assert!(
            (0.0..=100.0).contains(position),
            "{position} is not a valve position"
        );
    }
}
