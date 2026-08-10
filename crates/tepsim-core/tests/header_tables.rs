//! Cross-checks the measurement and manipulated-variable tables against the
//! header comments in the vendored `teprob.f`.
//!
//! Names and units get copied into plots, dataset columns, dashboards and
//! papers. A mistyped unit is invisible for years and then wrong everywhere, so
//! the tables are compared against the original text by machine rather than by
//! eye. Same discipline as the constants table in B-0006.
//!
//! An integration test because `tepsim-core` is `no_std` and cannot read files.

use std::path::PathBuf;

use tepsim_core::variables::{MANIPULATED, MEASUREMENTS};
use tepsim_core::{MeasIndex, MvIndex};

fn fortran_source() -> String {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../reference/fortran/teprob.f");
    std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("reading {}: {e}", path.display()))
}

/// Collapse runs of whitespace so the header's column alignment does not matter.
fn normalize(text: &str) -> String {
    text.split_whitespace().collect::<Vec<_>>().join(" ")
}

/// The header comment line for `XMEAS(n)` or `XMV(n)`, stripped of its `C`
/// marker and whitespace-normalised.
fn header_line(source: &str, prefix: &str, index: usize) -> String {
    let needle = format!("{prefix}({index})");
    let line = source
        .lines()
        .filter(|line| line.starts_with('C'))
        .map(|line| normalize(line.trim_start_matches('C')))
        .find(|line| line.starts_with(&format!("{needle} ")) || *line == needle)
        .unwrap_or_else(|| panic!("no header comment found for {needle}"));
    line[needle.len()..].trim().to_string()
}

#[test]
fn every_measurement_matches_the_fortran_header() {
    let source = fortran_source();
    let mut checked = 0;

    for info in &MEASUREMENTS {
        let actual = header_line(&source, "XMEAS", info.index);
        let expected = if info.analyzer.is_some() {
            // The sampled block declares "Mole %" once above the group rather
            // than on each line, so those lines carry only the description.
            info.description.to_string()
        } else {
            format!("{} {}", info.description, info.unit.fortran_spelling())
        };
        assert_eq!(
            actual, expected,
            "XMEAS({}) is documented as {actual:?} in teprob.f but the table \
             says {expected:?}",
            info.index
        );
        checked += 1;
    }

    assert_eq!(checked, MeasIndex::COUNT, "all 41 must be checked");
}

#[test]
fn every_manipulated_variable_matches_the_fortran_header() {
    let source = fortran_source();
    let mut checked = 0;

    for info in &MANIPULATED {
        let actual = header_line(&source, "XMV", info.index);
        // The first three carry a "(Corrected Order)" annotation recording the
        // 1991 revision. It is provenance, not part of the name.
        let actual = actual.replace("(Corrected Order)", "");
        assert_eq!(
            normalize(&actual),
            info.description,
            "XMV({}) is documented differently in teprob.f",
            info.index
        );
        checked += 1;
    }

    assert_eq!(checked, MvIndex::COUNT, "all 12 must be checked");
}

/// The annotation really is there on exactly the first three, and dropping it
/// is a deliberate choice rather than an accident of the comparison.
#[test]
fn exactly_the_first_three_manipulated_variables_are_annotated() {
    let source = fortran_source();
    let annotated: Vec<usize> = (1..=MvIndex::COUNT)
        .filter(|i| header_line(&source, "XMV", *i).contains("Corrected Order"))
        .collect();
    assert_eq!(annotated, vec![1, 2, 3]);
}

/// The analyser groupings are declared in prose above each block. Check the
/// stream each analyser samples, since that is what the table encodes.
#[test]
fn the_analyser_headings_name_the_streams_the_table_expects() {
    let source = fortran_source();
    let normalized: Vec<String> = source
        .lines()
        .filter(|line| line.starts_with('C'))
        .map(|line| normalize(line.trim_start_matches('C')))
        .collect();

    for heading in [
        "Reactor Feed Analysis (Stream 6)",
        "Purge Gas Analysis (Stream 9)",
        "Product Analysis (Stream 11)",
    ] {
        assert!(
            normalized.iter().any(|line| line == heading),
            "expected the header to declare {heading:?}"
        );
    }

    // Sampling frequencies, as the table records them.
    let count = |needle: &str| normalized.iter().filter(|l| l.contains(needle)).count();
    assert_eq!(
        count("Sampling Frequency = 0.1 hr"),
        2,
        "two fast analysers"
    );
    assert_eq!(
        count("Sampling Frequency = 0.25 hr"),
        1,
        "one slow analyser"
    );
    assert_eq!(count("Dead Time = 0.1 hr"), 2);
    assert_eq!(count("Dead Time = 0.25 hr"), 1);
}
