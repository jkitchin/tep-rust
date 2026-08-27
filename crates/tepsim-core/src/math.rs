//! The transcendental functions, and why they do not come from `f64`.
//!
//! Three functions reach this module. `exp`, first needed by the Antoine
//! vapour pressures at `teprob.f:485`; `pow`, first needed by the fractional
//! pressure orders at `teprob.f:508`; and `sqrt`, first needed by the
//! pressure-driven flows at `teprob.f:578`.
//!
//! # `sqrt` is here for a different reason from the other two
//!
//! IEEE-754 specifies `sqrt` exactly: it is correctly rounded, on every
//! conforming platform, and it compiles to a single hardware instruction. It
//! has none of the portability problem that `exp` and `pow` have, and it takes
//! no `libm-system` variant because there is nothing for one to fix.
//!
//! It is routed through this module anyway for a mundane reason: `f64::sqrt`
//! is a `std` method and this crate is `no_std`, so the call has to come from
//! the `libm` crate. That is a build-time detail with no numerical content,
//! and it is spelled out here because the obvious reading, that `sqrt` is
//! another portability hazard, is wrong and would spread caution where none is
//! needed.
//!
//! # Why the vendored crate rather than `f64::exp`
//!
//! `f64::exp` calls whatever libm the platform ships. glibc, musl and Apple's
//! do not agree to the last bit, and a wasm build has no platform libm at all.
//! Tier 9 asserts that the same BLAKE3 digest comes out of x86-64, aarch64 and
//! wasm32, so a host-dependent `exp` is not a rounding detail, it is a failed
//! gate. The vendored [`libm`] crate is the same code everywhere.
//!
//! # The price, measured
//!
//! On this machine gfortran's `DEXP` *is* Apple's libm. Over 1,500,005
//! arguments spanning the whole Antoine range this model reaches
//! (0.587 to 13.08), `f64::exp` matches gfortran on every one, and the
//! vendored crate differs on 9.945% of them, always by exactly one ULP low.
//!
//! So choosing determinism costs bit equality with the Fortran on about one
//! partial pressure in ten, at 1.1e-16 relative. Tier 2's gate is 1e-12, four
//! orders away. That is the trade `PLAN.org` makes deliberately, and B-0018
//! is where it first comes due.
//!
//! `pow` is measured the same way where it is used; see
//! `tests/tier2_kinetics.rs`.
//!
//! # `libm-system`, and what it is for
//!
//! With the non-default `libm-system` feature, [`exp`] becomes `f64::exp`.
//! This is never a shipping configuration: it forfeits the invariant above.
//! It exists so that validation can distinguish two very different findings.
//!
//! A Tier 2 difference under the default build could be transcendental
//! rounding or it could be a mistake in the algebra. Under `libm-system` the
//! transcendental agrees with gfortran exactly, so *any* remaining difference
//! is the algebra. The port therefore keeps a bit-exactness assertion for the
//! whole of Phase 2 rather than losing it to a 1e-12 tolerance at the first
//! call to `exp`, and a regression worth 1e-15 still fails a test instead of
//! hiding under the noise.
//!
//! Enabling it also makes the crate `std`, which is why it is not, and must
//! not become, a default.

/// The exponential.
///
/// See the module documentation: this is the vendored [`libm`] by default and
/// the platform libm under `libm-system`.
#[must_use]
#[cfg(not(feature = "libm-system"))]
pub fn exp(x: f64) -> f64 {
    libm::exp(x)
}

/// The exponential.
///
/// See the module documentation: this is the vendored [`libm`] by default and
/// the platform libm under `libm-system`.
#[must_use]
#[cfg(feature = "libm-system")]
pub fn exp(x: f64) -> f64 {
    x.exp()
}

/// `x` raised to the power `y`.
///
/// See the module documentation: this is the vendored [`libm`] by default and
/// the platform libm under `libm-system`.
///
/// Only for the two *fractional* exponents in the model, `teprob.f:508-509`.
/// Integer powers such as `AGSP**2` and `(FTM(8)/3528.73)**4` must expand to
/// multiplication instead, because that is what gfortran emits for them, and
/// routing them through `pow` would change the last bits.
#[must_use]
#[cfg(not(feature = "libm-system"))]
pub fn pow(x: f64, y: f64) -> f64 {
    libm::pow(x, y)
}

/// `x` raised to the power `y`.
///
/// See the module documentation: this is the vendored [`libm`] by default and
/// the platform libm under `libm-system`.
///
/// Only for the two *fractional* exponents in the model, `teprob.f:508-509`.
/// Integer powers such as `AGSP**2` and `(FTM(8)/3528.73)**4` must expand to
/// multiplication instead, because that is what gfortran emits for them, and
/// routing them through `pow` would change the last bits.
#[must_use]
#[cfg(feature = "libm-system")]
pub fn pow(x: f64, y: f64) -> f64 {
    x.powf(y)
}

/// The square root.
///
/// Unlike [`exp`] and [`pow`] this has a single correct answer on every
/// platform, and it is the same function in both build configurations. See the
/// module documentation for why it goes through here at all.
#[must_use]
pub fn sqrt(x: f64) -> f64 {
    libm::sqrt(x)
}

/// Whether this build uses the platform libm rather than the vendored one.
///
/// Reported by validation harnesses, so that a recorded number always says
/// which `exp` produced it.
pub const USES_SYSTEM_LIBM: bool = cfg!(feature = "libm-system");

#[cfg(test)]
mod tests {
    use super::*;

    /// A handful of values a reader can check independently, so that a
    /// mis-wired feature flag shows up as a wrong number rather than as a
    /// silently different library.
    #[test]
    fn exp_is_the_exponential() {
        assert_eq!(exp(0.0).to_bits(), 1.0_f64.to_bits());
        assert!(exp(f64::NEG_INFINITY) == 0.0);
        assert!(exp(f64::INFINITY).is_infinite());
        assert!(exp(f64::NAN).is_nan());
    }

    /// `sqrt` is exact, so it must agree with the platform in *both*
    /// configurations. If this ever fails, IEEE-754 is not being honoured and
    /// the reasoning in the module documentation no longer holds.
    #[test]
    fn sqrt_is_exact_in_both_configurations() {
        assert_eq!(sqrt(4.0).to_bits(), 2.0_f64.to_bits());
        assert_eq!(sqrt(0.0).to_bits(), 0.0_f64.to_bits());
        assert!(sqrt(-1.0).is_nan());
        // A value whose root is not representable, where a sloppy
        // implementation would drift.
        let x = 2.0_f64;
        assert_eq!(sqrt(x).to_bits(), 0x3FF6A09E667F3BCD);
    }

    #[test]
    fn pow_is_exponentiation() {
        assert_eq!(pow(2.0, 10.0).to_bits(), 1024.0_f64.to_bits());
        assert_eq!(pow(9.0, 0.5).to_bits(), 3.0_f64.to_bits());
        assert_eq!(pow(1.0, f64::NAN).to_bits(), 1.0_f64.to_bits());
        assert_eq!(pow(0.0, 1.1544).to_bits(), 0.0_f64.to_bits());
    }

    /// The two implementations differ, and by how much, pinned at one argument
    /// a reader can check by hand.
    ///
    /// `exp(1)` is `e`, whose correctly rounded double is
    /// [`core::f64::consts::E`]. The platform libm returns it; the vendored
    /// crate returns one ULP above. That is the whole shape of the trade this
    /// module makes, in one assertion, and it is written as an equality on
    /// bits so that a future version of either library changing its rounding
    /// is a test failure rather than a silent shift in every recorded number.
    #[test]
    fn the_vendored_exp_is_one_ulp_high_at_one_and_the_platform_one_is_exact() {
        let correctly_rounded = core::f64::consts::E.to_bits();
        let actual = exp(1.0).to_bits();
        if USES_SYSTEM_LIBM {
            assert_eq!(
                actual, correctly_rounded,
                "the platform libm rounds e correctly"
            );
        } else {
            assert_eq!(
                actual,
                correctly_rounded + 1,
                "the vendored libm is expected to be exactly one ULP high at e"
            );
        }
    }

    /// The default build must be the deterministic one. If this ever fails,
    /// `libm-system` has leaked into `default` and Tier 9 is no longer being
    /// tested by the ordinary gate.
    #[test]
    #[cfg(not(feature = "libm-system"))]
    fn the_default_build_uses_the_vendored_libm() {
        const { assert!(!USES_SYSTEM_LIBM) };
    }
}
