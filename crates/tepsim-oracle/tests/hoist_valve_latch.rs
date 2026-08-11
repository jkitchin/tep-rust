//! Proves the valve-command latch may be hoisted out of the end of `TEFUNC`
//! and into the port's pre-phase.
//!
//! B-0012 split the right-hand side into an impure pre-phase, a pure
//! derivative, and an impure post-phase. The latch at `teprob.f:793-804` is the
//! awkward one: the original runs it at the *end* of the routine, immediately
//! before the valve derivatives, but the port needs it at the start, because
//! `derivatives` must be pure and the latch mutates `VCV`.
//!
//! That is only safe if nothing between the `IDV` clamp and the latch writes
//! anything the latch reads. The latch reads `TIME`, `XMV`, `VST`, `IVST` and,
//! through `IVST`, `IDV`. So the claim is:
//!
//! > No statement in `teprob.f:345-792` assigns to `IDV`, `XMV`, `VCV`, `VST`
//! > or `IVST`.
//!
//! That is a claim about four hundred and fifty lines of Fortran, and reading
//! it by eye is exactly how a silent divergence gets into a port. So it is
//! checked mechanically, against the vendored source, and it will fail if the
//! reference is ever replaced with a variant that behaves differently.
//!
//! This test needs no Fortran compiler: it reads the source. It lives here
//! rather than in `tepsim-core` because `reference/` is the oracle crate's
//! business, and the other reference-integrity tests are already here.

use std::path::PathBuf;

/// Everything the hoisted block reads.
const READS: [&str; 5] = ["IDV", "XMV", "VCV", "VST", "IVST"];

/// The window between the `IDV` clamp (`341-344`) and the hoisted block,
/// one-based and inclusive, as `teprob.f` numbers its lines.
///
/// The block begins at 793, not 794. The first draft of this test said 794 and
/// failed, correctly, on `IVST(10)=IDV(14)`: an off-by-one found by the machine
/// rather than by the eye that wrote it, which is the entire argument for
/// checking this mechanically.
const WINDOW: (usize, usize) = (345, 792);

/// The hoisted block: `IVST` from `IDV` at 793-798, then the `VCV` latch at
/// 799-804. Line 805, `YP(I+38)`, shares the `DO 9020` loop with the latch but
/// is the valve derivative and stays in the pure phase, so the port splits the
/// loop in two.
const HOISTED: (usize, usize) = (793, 804);

fn teprob() -> String {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(|p| p.parent())
        .expect("workspace root")
        .join("reference/fortran/teprob.f");
    std::fs::read_to_string(&root).unwrap_or_else(|e| panic!("reading {}: {e}", root.display()))
}

/// The executable part of a fixed-form Fortran line: not a comment, and past
/// the five-column label field and the continuation column.
fn statement(line: &str) -> Option<String> {
    let first = line.chars().next()?;
    if matches!(first, 'C' | 'c' | '*' | '!') {
        return None;
    }
    let body: String = line.chars().skip(6).collect();
    if body.trim().is_empty() {
        None
    } else {
        // Fortran ignores spaces inside statements; removing them means
        // `IVST (10) = ...` cannot slip past a pattern written without them.
        Some(body.replace(' ', "").to_uppercase())
    }
}

/// Does this statement assign to `name`, as a whole array or an element?
///
/// Looks for `NAME=` or `NAME(...)=`, with the name not preceded by another
/// identifier character so that `IVST` does not match inside a longer name.
/// A comparison never produces this shape: Fortran spells equality `.EQ.`, and
/// the relational operators all begin with a dot.
fn assigns_to(statement: &str, name: &str) -> bool {
    let mut from = 0;
    while let Some(at) = statement[from..].find(name) {
        let start = from + at;
        let end = start + name.len();
        from = end;

        let preceded = start
            .checked_sub(1)
            .and_then(|i| statement.as_bytes().get(i))
            .is_some_and(|c| c.is_ascii_alphanumeric() || *c == b'_');
        if preceded {
            continue;
        }

        let rest = &statement[end..];
        let after_subscript = if rest.starts_with('(') {
            rest.find(')').map(|close| &rest[close + 1..])
        } else {
            Some(rest)
        };
        if after_subscript.is_some_and(|r| r.starts_with('=')) {
            return true;
        }
    }
    false
}

#[test]
fn nothing_between_the_clamp_and_the_latch_writes_what_the_latch_reads() {
    let source = teprob();
    let lines: Vec<&str> = source.lines().collect();
    assert!(
        lines.len() >= 1589,
        "teprob.f has {} lines; the ranges in this test assume 1589",
        lines.len()
    );

    let mut writes = Vec::new();
    let mut statements = 0;
    for (number, line) in lines
        .iter()
        .enumerate()
        .map(|(i, l)| (i + 1, l))
        .filter(|(n, _)| (WINDOW.0..=WINDOW.1).contains(n))
    {
        let Some(statement) = statement(line) else {
            continue;
        };
        statements += 1;
        for name in READS {
            if assigns_to(&statement, name) {
                writes.push(format!("  teprob.f:{number}: {}", line.trim_end()));
            }
        }
    }

    println!(
        "scanned {statements} statements over teprob.f:{}-{} for writes to {READS:?}",
        WINDOW.0, WINDOW.1
    );
    assert!(
        writes.is_empty(),
        "the valve-latch hoist is NOT safe: {} statement(s) in the window \
         write something the latch reads, so running it early would change \
         numbers.\n{}",
        writes.len(),
        writes.join("\n")
    );
    assert!(
        statements > 300,
        "only {statements} statements were scanned, so the window or the \
         parser is wrong and this test proves nothing"
    );
}

/// The detector has to be able to see a write, or the test above is vacuous.
///
/// The hoisted block itself is full of them, so pointing the same scan at
/// `teprob.f:793-804` must find several.
#[test]
fn the_scan_finds_the_writes_that_are_actually_there() {
    let source = teprob();
    let lines: Vec<&str> = source.lines().collect();

    let mut found = Vec::new();
    for (number, line) in lines.iter().enumerate().map(|(i, l)| (i + 1, l)) {
        if !(HOISTED.0..=HOISTED.1).contains(&number) {
            continue;
        }
        let Some(statement) = statement(line) else {
            continue;
        };
        for name in READS {
            if assigns_to(&statement, name) {
                found.push((number, name));
            }
        }
    }

    println!("writes inside the hoisted block: {found:?}");
    assert!(
        found.iter().any(|(_, name)| *name == "IVST"),
        "the scan missed the IVST assignments at teprob.f:793-798"
    );
    assert!(
        found.iter().any(|(_, name)| *name == "VCV"),
        "the scan missed the VCV latch at teprob.f:799-804"
    );
    assert!(
        found.len() >= 6,
        "expected at least six writes in the hoisted block, found {}",
        found.len()
    );
}

/// A conditional assignment on one line still counts as a write, and a
/// comparison still does not.
#[test]
fn the_detector_tells_assignment_from_comparison() {
    assert!(assigns_to("IVST(10)=IDV(14)", "IVST"));
    assert!(assigns_to("IF(VCV(I).LT.0.0)VCV(I)=0.0", "VCV"));
    assert!(assigns_to("IDV(I)=1", "IDV"));
    assert!(!assigns_to("IF(IDV(I).GT.0)THEN", "IDV"));
    assert!(!assigns_to("DABS(VCV(I)-XMV(I)).GT.VST(I)*IVST(I)", "VCV"));
    assert!(!assigns_to("DABS(VCV(I)-XMV(I)).GT.VST(I)*IVST(I)", "XMV"));
    assert!(!assigns_to("YP(I+38)=(VCV(I)-VPOS(I))/VTAU(I)", "VCV"));
    // A longer name that merely contains one of ours must not match.
    assert!(!assigns_to("MYVST(3)=1", "VST"));
}
