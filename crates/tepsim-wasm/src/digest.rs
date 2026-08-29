//! A small, exact, dependency-free fingerprint over a stream of `f64`s.
//!
//! # Where the hash lives now
//!
//! [`Fnv1a64`] used to be defined here. It is [`tepsim::tier9::Fnv1a64`], and
//! this is a re-export, because the same hash now serves two jobs that must
//! agree exactly: the running checksum these bindings hand a browser chunk by
//! chunk, and the canonical Tier 9 digest every target checks itself against.
//! Two implementations of the same algorithm would eventually stop being the
//! same algorithm, and the failure would look like a cross-platform
//! disagreement.
//!
//! # Why not BLAKE3
//!
//! `PLAN.org` names BLAKE3 for Tier 9. It is FNV-1a, deliberately, and the
//! reasoning is in [`tepsim::tier9`]. The short version: the digest has to be
//! computed inside the wasm module, which has a size budget, and what `xtask`
//! ends up comparing is two 64-bit numbers from two machines, which a stronger
//! hash would not make more comparable.
//!
//! FNV-1a over `f64::to_bits` is integer arithmetic end to end, so the digest
//! cannot itself become a source of the cross-platform disagreement it exists
//! to detect. It is not cryptographic and not collision-resistant against an
//! adversary. It catches drift, which is the failure this project actually has.

pub use tepsim::tier9::Fnv1a64;

/// Lower-case hex, fixed width, so two digests line up when printed one above
/// the other. That is the only way anybody actually compares them.
#[must_use]
pub fn hex64(value: u64) -> String {
    format!("{value:016x}")
}
