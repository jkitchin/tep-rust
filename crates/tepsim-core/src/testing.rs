//! Test helpers shared across the crate.
//!
//! Exists for one reason: this is numerical code, so `clippy::float_cmp` fires
//! on every `assert_eq!` over an `f64`. Suppressing the lint would be wrong,
//! because it catches real mistakes elsewhere, and an approximate comparison
//! would be wrong too, because these values are meant to be *exact*. Comparing
//! bit patterns says precisely what is meant and satisfies the lint honestly.

/// Assert two floats are bit-identical.
///
/// Distinguishes `0.0` from `-0.0` and never treats `NaN` as equal to itself,
/// both of which are the desired behaviour when checking stored constants.
#[track_caller]
pub(crate) fn assert_exact(actual: f64, expected: f64, what: &str) {
    assert_eq!(
        actual.to_bits(),
        expected.to_bits(),
        "{what}: got {actual:?}, expected {expected:?}"
    );
}

/// Assert two floats are not bit-identical.
#[track_caller]
pub(crate) fn assert_not_exact(actual: f64, expected: f64, what: &str) {
    assert_ne!(
        actual.to_bits(),
        expected.to_bits(),
        "{what}: expected {actual:?} to differ from {expected:?}"
    );
}
