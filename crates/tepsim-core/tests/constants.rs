//! Re-derives every thermodynamic constant from the vendored Fortran and
//! compares it against the transcribed table.
//!
//! Independent of how the table was produced. The parser here reads
//! `reference/fortran/teprob.f`, decides from the presence or absence of a `D`
//! suffix whether each literal is single or double precision, and computes the
//! `f64` the original would store. 112 constants is well past the number a
//! person checks reliably by eye.
//!
//! The complementary check lives in `tepsim-oracle` and compares against the
//! values gfortran actually stored in `COMMON/CONST/`. That one is decisive,
//! since it does not share this file's reading of Fortran literal semantics.

use std::path::PathBuf;

use tepsim_core::constants;

/// One assignment recovered from the source.
#[derive(Debug)]
struct Parsed {
    array: String,
    index: usize,
    literal: String,
    line: usize,
    /// True when the literal carries a `D` exponent and so is double precision.
    is_double: bool,
    /// What the original stores, derived from the literal and its suffix.
    value: f64,
}

fn source() -> String {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../reference/fortran/teprob.f");
    std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("reading {}: {e}", path.display()))
}

const ARRAYS: [&str; 14] = [
    "AVP", "BVP", "CVP", "AH", "BH", "CH", "AG", "BG", "CG", "AV", "AD", "BD", "CD", "XMW",
];

fn parse_assignments(text: &str) -> Vec<Parsed> {
    let mut out = Vec::new();
    for (n, raw) in text.lines().enumerate() {
        let line = n + 1;
        let Some(rest) = raw.strip_prefix("      ") else {
            continue;
        };
        let Some((lhs, literal)) = rest.split_once('=') else {
            continue;
        };
        let Some((array, idx)) = lhs.split_once('(') else {
            continue;
        };
        let Some(idx) = idx.strip_suffix(')') else {
            continue;
        };
        if !ARRAYS.contains(&array) {
            continue;
        }
        let Ok(index) = idx.parse::<usize>() else {
            continue;
        };
        let literal = literal.trim().to_string();

        // A `D` exponent means double precision. Anything else, including a
        // bare decimal or an `E` exponent, is single.
        let (mantissa, exponent, is_double) = match literal.find(['D', 'd']) {
            Some(pos) => (
                &literal[..pos],
                literal[pos + 1..].parse::<i32>().expect("D exponent"),
                true,
            ),
            None => (literal.as_str(), 0, false),
        };
        let normalized = normalize(mantissa);
        let value = if is_double {
            // Compose in f64 exactly as the compiler would.
            format!("{normalized}e{exponent}")
                .parse::<f64>()
                .expect("double literal")
        } else {
            // Round to f32 first, then widen. This is the whole point.
            let narrow: f32 = normalized.parse().expect("single literal");
            f64::from(narrow)
        };

        out.push(Parsed {
            array: array.to_string(),
            index,
            literal,
            line,
            is_double,
            value,
        });
    }
    out
}

fn normalize(m: &str) -> String {
    let (sign, body) = m.strip_prefix('-').map_or(("", m), |b| ("-", b));
    let mut body = body.to_string();
    if body.ends_with('.') {
        body.push('0');
    }
    if body.starts_with('.') {
        body.insert(0, '0');
    }
    if !body.contains('.') {
        body.push_str(".0");
    }
    format!("{sign}{body}")
}

fn table(array: &str) -> &'static [f64; 8] {
    match array {
        "AVP" => constants::AVP.as_array(),
        "BVP" => constants::BVP.as_array(),
        "CVP" => constants::CVP.as_array(),
        "AH" => constants::AH.as_array(),
        "BH" => constants::BH.as_array(),
        "CH" => constants::CH.as_array(),
        "AG" => constants::AG.as_array(),
        "BG" => constants::BG.as_array(),
        "CG" => constants::CG.as_array(),
        "AV" => constants::AV.as_array(),
        "AD" => constants::AD.as_array(),
        "BD" => constants::BD.as_array(),
        "CD" => constants::CD.as_array(),
        "XMW" => constants::XMW.as_array(),
        other => panic!("unknown array {other}"),
    }
}

#[test]
fn every_constant_matches_the_fortran_literal_and_its_precision() {
    let parsed = parse_assignments(&source());
    assert_eq!(
        parsed.len(),
        112,
        "expected 112 assignments across 14 arrays of 8; the source changed"
    );

    let mut wrong = Vec::new();
    for p in &parsed {
        let ours = table(&p.array)[p.index - 1];
        if ours.to_bits() != p.value.to_bits() {
            wrong.push(format!(
                "  {}({}) at teprob.f:{}: literal {} is {} precision, so it stores {:?}, \
                 but the table has {:?}",
                p.array,
                p.index,
                p.line,
                p.literal,
                if p.is_double { "double" } else { "single" },
                p.value,
                ours
            ));
        }
    }
    assert!(
        wrong.is_empty(),
        "{} of {} constants disagree with the Fortran:\n{}",
        wrong.len(),
        parsed.len(),
        wrong.join("\n")
    );
}

/// The split is a fact about the source, not a guess. If it moves, something
/// changed and the transcription needs revisiting.
#[test]
fn the_precision_split_is_62_single_and_50_double() {
    let parsed = parse_assignments(&source());
    let double = parsed.iter().filter(|p| p.is_double).count();
    let single = parsed.len() - double;
    assert_eq!((single, double), (62, 50));
}

/// The specific canary from B-0004b, kept as its own test so a failure names it.
#[test]
fn xmw_of_b_is_the_single_precision_value() {
    let stored = constants::XMW.as_array()[1];
    assert_eq!(
        stored.to_bits(),
        f64::from(25.4_f32).to_bits(),
        "XMW(2) must be 25.399999618530273"
    );
    assert_ne!(
        stored.to_bits(),
        25.4_f64.to_bits(),
        "XMW(2) must NOT be the nearest double to 25.4; that is wrong by 1.5e-8 \
         relative, five orders of magnitude past the Tier 1 tolerance"
    );
}

/// All 24 Antoine coefficients are single precision, and they feed `exp`, so
/// they are the constants where the distinction matters most.
#[test]
fn every_antoine_coefficient_is_single_precision() {
    let parsed = parse_assignments(&source());
    for p in parsed
        .iter()
        .filter(|p| matches!(p.array.as_str(), "AVP" | "BVP" | "CVP"))
    {
        assert!(
            !p.is_double,
            "{}({}) at teprob.f:{} unexpectedly carries a D suffix",
            p.array, p.index, p.line
        );
    }
    // And at least one of them is genuinely not representable in f32.
    let avp4 = constants::AVP.as_array()[3];
    assert_ne!(avp4.to_bits(), 15.92_f64.to_bits());
    assert_eq!(avp4.to_bits(), f64::from(15.92_f32).to_bits());
}
