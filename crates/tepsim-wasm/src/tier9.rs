//! Tier 9's wasm32 half: the determinism table, exported without glue.
//!
//! [`tepsim::tier9`] holds the table of fixed scenarios and the digest each is
//! committed to produce. A native build checks itself against those constants
//! in `crates/tepsim/tests/tier9.rs`. This module is how a `wasm32` build does
//! the same, and it is the only thing standing between "the port is
//! deterministic across architectures" and an assertion nobody has run.
//!
//! # Why these are `extern "C"` and not `#[wasm_bindgen]`
//!
//! `wasm-bindgen` produces a far nicer interface, and it needs a code
//! generator: `wasm-bindgen-cli`, at a version that matches the crate exactly.
//! That is a reasonable thing to ask of someone building the browser app and an
//! unreasonable thing to put between a continuous-integration run and the one
//! check that says the numbers are the same everywhere.
//!
//! These exports need none of it. Any runtime that can instantiate a
//! module and call a function reads the digests back: `node`, `wasmtime`,
//! `wasmer`, a browser's bare `WebAssembly.instantiate`. They take and return
//! only `u32` and `u64`, so nothing crosses the boundary that needs a memory
//! layout agreed in advance, and the module's linear memory never has to be
//! read from outside.
//!
//! One consequence is worth knowing before pointing a runtime at the module:
//! the rest of this crate is `#[wasm_bindgen]`, so the compiled artifact still
//! *declares* imports from `__wbindgen_placeholder__`. Nothing here calls them,
//! but a runtime must supply something for each before it can instantiate.
//! `cargo xtask tier9` stubs them with functions that throw, which turns "the
//! Tier 9 path touched JavaScript" from a silent success into a loud failure.
//!
//! # What it costs the bundle
//!
//! `PLAN.org` puts a size budget on this crate, so the cost was measured rather
//! than assumed. Adding this module and the Tier 9 table moves the
//! `release-wasm` artifact from 528,736 to 533,025 bytes raw, and from 161,210
//! to 165,310 bytes after `wasm-bindgen`: about 4.1 kB, 2.5 percent. Most of it
//! is the table's `covers` strings, which are worth keeping because a failing
//! digest is useless without knowing what the case exercises. A BLAKE3
//! implementation, the alternative `PLAN.org` names, would cost several times
//! that on its own; see [`tepsim::tier9`] for the rest of that argument.
//!
//! # Names
//!
//! A `no_mangle` symbol is global to the linked artifact, so every one of these
//! carries the `tepsim_wasm_` prefix.

use tepsim::tier9::{self, CASES};

/// How many cases the table holds.
///
/// A caller loops `0..count` rather than hard-coding six, so adding a case to
/// [`tepsim::tier9::CASES`] extends the wasm check with no change here and no
/// change in `xtask`.
//
// SAFETY: `no_mangle` is unsafe because the symbol could collide with another
// of the same name and silently redirect calls. The prefix makes a collision
// implausible. The function takes no arguments, touches no shared state, and
// returns a plain integer, so there is no aliasing or lifetime obligation for a
// caller to uphold. The same argument covers the three below.
#[unsafe(no_mangle)]
pub extern "C" fn tepsim_wasm_tier9_case_count() -> u32 {
    // Saturating rather than `as`: a table longer than `u32::MAX` is
    // impossible, and a silent truncation would under-report the work done.
    u32::try_from(CASES.len()).unwrap_or(u32::MAX)
}

/// Run case `index` here, on this target, and return the digest it produces.
///
/// Zero for an index past the end of the table. Zero is not a digest any run
/// produces, and the caller is expected to have asked
/// [`tepsim_wasm_tier9_case_count`] first.
#[unsafe(no_mangle)]
pub extern "C" fn tepsim_wasm_tier9_case_digest(index: u32) -> u64 {
    match CASES.get(index as usize) {
        Some(case) => case.compute(),
        None => 0,
    }
}

/// The digest case `index` is *committed* to produce, compiled in from the
/// table.
///
/// Exported as well as computed so the module is self-checking: a browser page,
/// or any runtime, can compare the two without being told the answer
/// separately. `xtask` compares three ways, this against its own copy of the
/// table included, which is what would catch a stale `.wasm` left in `target/`
/// from an older commit.
#[unsafe(no_mangle)]
pub extern "C" fn tepsim_wasm_tier9_expected_digest(index: u32) -> u64 {
    match CASES.get(index as usize) {
        Some(case) => case.digest,
        None => 0,
    }
}

/// The whole table as one number, computed here. Compare against
/// [`tepsim::tier9::SUITE_DIGEST`], which
/// [`tepsim_wasm_tier9_expected_suite_digest`] exports.
#[unsafe(no_mangle)]
pub extern "C" fn tepsim_wasm_tier9_suite_digest() -> u64 {
    tier9::suite_digest()
}

/// The committed value of the suite digest, compiled in.
#[unsafe(no_mangle)]
pub extern "C" fn tepsim_wasm_tier9_expected_suite_digest() -> u64 {
    tier9::SUITE_DIGEST
}
