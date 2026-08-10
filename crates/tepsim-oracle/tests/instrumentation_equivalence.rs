//! Proves the build-time instrumentation changed no numbers.
//!
//! The oracle links a *rewritten* copy of `teprob.f`, because the original
//! keeps the shutdown flag in a local variable that cannot be linked against.
//! Every validation tier is measured against that rewritten copy, so the claim
//! that it computes the same thing as the pristine original is load-bearing and
//! must be measured rather than argued.
//!
//! Both copies define the same symbols, so they cannot be linked into one
//! process. The pristine one is built as a standalone executable by `build.rs`
//! and reports its answers as raw IEEE bit patterns; this test runs the same
//! sequence in-process against the instrumented build and compares bit for bit.

#![cfg(feature = "oracle")]

use std::process::Command;

use tepsim_oracle::{N_STATES, Oracle};

/// The seed compiled into `teprob.f:1187`, matching the probe driver.
const SEED: f64 = 4651207995.0;
const ONE_SECOND_H: f64 = 1.0 / 3600.0;

fn pristine_answers() -> Vec<u64> {
    let probe = env!("TEP_ORACLE_PRISTINE_PROBE");
    let output = Command::new(probe)
        .output()
        .unwrap_or_else(|e| panic!("running the pristine probe at {probe}: {e}"));
    assert!(
        output.status.success(),
        "pristine probe exited with {}",
        output.status
    );
    let text = String::from_utf8(output.stdout).expect("probe output is ASCII");
    text.split_whitespace()
        .map(|tok| {
            tok.parse::<i64>()
                .unwrap_or_else(|e| panic!("parsing {tok:?} from the probe: {e}"))
                as u64
        })
        .collect()
}

#[test]
fn instrumentation_changes_no_numbers() {
    let pristine = pristine_answers();
    assert_eq!(
        pristine.len(),
        N_STATES + 1,
        "the probe reports 50 derivatives and the final generator word"
    );

    let mut oracle = Oracle::lock();
    let (_, yy) = oracle.init();
    oracle.set_rng(SEED);
    let yp = oracle.derivatives(ONE_SECOND_H, &yy);
    let g = oracle.rng();

    let mut differing = Vec::new();
    for (i, (got, want)) in yp.iter().zip(&pristine).enumerate() {
        if got.to_bits() != *want {
            differing.push(format!(
                "  YP({}) instrumented {:?} vs pristine {:?}",
                i + 1,
                got,
                f64::from_bits(*want)
            ));
        }
    }
    assert!(
        differing.is_empty(),
        "the instrumented build disagrees with the pristine original on {} of \
         50 derivatives:\n{}\n\nThe instrumentation is supposed to hoist ISD \
         into a COMMON block and split off an IRAND local, neither of which can \
         change a number. If this fails, the rewrite in instrument.rs is wrong \
         and every validation tier measured against this oracle is suspect.",
        differing.len(),
        differing.join("\n")
    );

    assert_eq!(
        g.to_bits(),
        pristine[N_STATES],
        "both builds must consume exactly the same number of random draws"
    );
}
