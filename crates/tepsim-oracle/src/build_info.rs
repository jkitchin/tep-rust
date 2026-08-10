//! What actually compiled the oracle.
//!
//! Every Tier 1 and Tier 2 number recorded in `LOG.org` is a measurement
//! against a specific compiler and a specific flag set. Recording both
//! alongside the numbers is what lets a later session tell a genuine regression
//! apart from a toolchain upgrade, which otherwise look identical and cost a
//! whole session to distinguish.

/// The flags `build.rs` passed to gfortran, space separated.
///
/// Empty when the `oracle` feature is off, since nothing was compiled.
pub const FORTRAN_FLAGS: &str = match option_env!("TEP_ORACLE_FORTRAN_FLAGS") {
    Some(flags) => flags,
    None => "",
};

/// The gfortran version that compiled the reference, as `major.minor.patch`.
///
/// Empty when the `oracle` feature is off.
pub const GFORTRAN_VERSION: &str = match option_env!("TEP_ORACLE_GFORTRAN_VERSION") {
    Some(version) => version,
    None => "",
};

/// Flags that permit the compiler to reassociate floating-point expressions.
///
/// Any of these would let gfortran rewrite the arithmetic, silently changing
/// the values the entire validation ladder is measured against. Their absence
/// is asserted rather than assumed.
pub const FORBIDDEN_FLAGS: &[&str] = &[
    "-ffast-math",
    "-funsafe-math-optimizations",
    "-fassociative-math",
    "-freciprocal-math",
    "-ffinite-math-only",
    "-Ofast",
];

#[cfg(all(test, feature = "oracle"))]
mod tests {
    use super::*;

    #[test]
    fn the_flag_set_is_exactly_what_we_pinned() {
        let actual: Vec<&str> = FORTRAN_FLAGS.split_whitespace().collect();
        assert_eq!(
            actual,
            vec![
                "-c",
                "-O0",
                "-fno-fast-math",
                "-fno-unsafe-math-optimizations",
                "-fPIC",
                "-std=legacy",
            ],
            "the Fortran flag set changed. Every recorded Tier 1 and Tier 2 \
             number was measured with the old set, so this is a deliberate \
             re-baseline: update the numbers in LOG.org in the same commit."
        );
    }

    #[test]
    fn no_flag_permits_reassociation() {
        for forbidden in FORBIDDEN_FLAGS {
            assert!(
                !FORTRAN_FLAGS.split_whitespace().any(|f| f == *forbidden),
                "{forbidden} lets gfortran reassociate floating-point \
                 expressions, which would change the reference values the port \
                 is validated against"
            );
        }
    }

    #[test]
    fn the_compiler_identified_itself() {
        assert!(
            !GFORTRAN_VERSION.is_empty(),
            "gfortran version unknown; the golden trace in B-0004c needs it to \
             report a toolchain change rather than a regression"
        );
    }
}
