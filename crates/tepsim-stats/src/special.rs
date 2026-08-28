//! The two special functions the t-distribution needs, and the `libm` shims.
//!
//! Both are standard; both are here rather than pulled from a crate so that
//! the equivalence claim rests on code this project tests. Each is checked
//! against exact identities, not against another implementation.

/// `sqrt`, from the vendored `libm`.
///
/// `f64::sqrt` is a `std` method and this crate is `no_std`. Both are the
/// IEEE-754 correctly-rounded square root, so unlike `exp` and `pow` there is
/// no portability question here; the shim exists only to satisfy `no_std`.
#[must_use]
pub fn sqrt(x: f64) -> f64 {
    libm::sqrt(x)
}

/// `exp`, from the vendored `libm`.
#[must_use]
pub fn exp(x: f64) -> f64 {
    libm::exp(x)
}

/// `ln`, from the vendored `libm`.
#[must_use]
pub fn ln(x: f64) -> f64 {
    libm::log(x)
}

/// The natural logarithm of the gamma function, for `x > 0`.
///
/// Lanczos approximation with `g = 7` and nine coefficients, the parameters
/// tabulated by Press et al., *Numerical Recipes in C*, 2nd ed., section 6.1.
/// Accurate to about 15 significant digits over the positive reals, which is
/// all the incomplete beta needs.
///
/// The logarithm rather than the value because `gamma(172)` overflows `f64`
/// while `ln_gamma(172)` is about 719, and the incomplete beta only ever wants
/// the log.
#[must_use]
pub fn ln_gamma(x: f64) -> f64 {
    /// Lanczos coefficients for `g = 7`, `n = 9`, transcribed verbatim from
    /// the published table. Not regrouped and not rounded: a coefficient
    /// retyped to a friendlier number of digits is a different function.
    #[allow(
        clippy::inconsistent_digit_grouping,
        clippy::excessive_precision,
        reason = "transcribed verbatim; regrouping or rounding changes the value"
    )]
    const C: [f64; 9] = [
        0.99999999999980993,
        676.5203681218851,
        -1259.1392167224028,
        771.32342877765313,
        -176.61502916214059,
        12.507343278686905,
        -0.13857109526572012,
        9.9843695780195716e-6,
        1.5056327351493116e-7,
    ];
    const G: f64 = 7.0;

    // The reflection formula for the left half plane. Not reachable from the
    // incomplete beta, whose arguments are always positive, but leaving it out
    // would make the function quietly wrong for a caller who did not know.
    if x < 0.5 {
        return ln(core::f64::consts::PI / libm::sin(core::f64::consts::PI * x))
            - ln_gamma(1.0 - x);
    }

    let x = x - 1.0;
    let mut series = C[0];
    for (index, coefficient) in C.iter().enumerate().skip(1) {
        series += coefficient / (x + index as f64);
    }
    let t = x + G + 0.5;
    0.5 * ln(2.0 * core::f64::consts::PI) + (x + 0.5) * ln(t) - t + ln(series)
}

/// The regularised incomplete beta function `I_x(a, b)`.
///
/// Defined as `B(x; a, b) / B(a, b)`, so it runs from 0 at `x = 0` to 1 at
/// `x = 1` and is the CDF of the Beta(a, b) distribution.
///
/// Evaluated by the continued fraction of Lentz, as given in Press et al.,
/// *Numerical Recipes in C*, 2nd ed., section 6.4. The fraction converges
/// quickly for `x < (a + 1) / (a + b + 2)` and the symmetry
/// `I_x(a, b) = 1 - I_{1-x}(b, a)` covers the rest, which is why the guard
/// below exists.
///
/// # Panics
///
/// Never. Out-of-domain arguments return `NaN`, because a validation routine
/// that panics mid-battery loses every number computed before it.
#[must_use]
pub fn regularized_incomplete_beta(a: f64, b: f64, x: f64) -> f64 {
    // Written so that a NaN argument falls into the guard rather than past it.
    if a.is_nan() || b.is_nan() || a <= 0.0 || b <= 0.0 || !(0.0..=1.0).contains(&x) {
        return f64::NAN;
    }
    // Exact comparisons on purpose: these are the endpoints of the domain, not
    // approximate targets, and the closed forms there are exact.
    if x.to_bits() == 0.0_f64.to_bits() || x.to_bits() == (-0.0_f64).to_bits() {
        return 0.0;
    }
    if x.to_bits() == 1.0_f64.to_bits() {
        return 1.0;
    }

    // Invariant under (a, b, x) -> (b, a, 1-x), which is what lets both
    // branches below share it.
    let front = exp(ln_gamma(a + b) - ln_gamma(a) - ln_gamma(b) + a * ln(x) + b * ln(1.0 - x));

    // The fraction converges quickly only on one side of this point, so the
    // other side is evaluated through the symmetry
    // `I_x(a, b) = 1 - I_{1-x}(b, a)`.
    //
    // Applied to the *fraction* rather than by calling this function again.
    // The recursive form is the obvious one and it does not terminate: at
    // `a == b` the threshold is exactly 1/2, so `I_{1/2}(a, a)` swaps to
    // itself forever. That is not an exotic input, it is Student's t with one
    // degree of freedom at `t = 1`.
    if x < (a + 1.0) / (a + b + 2.0) {
        front * continued_fraction(a, b, x) / a
    } else {
        1.0 - front * continued_fraction(b, a, 1.0 - x) / b
    }
}

/// The modified Lentz evaluation of the beta continued fraction.
///
/// `TINY` guards against a zero denominator, which Lentz's method handles by
/// substituting a number small enough not to matter and large enough not to
/// divide by zero.
fn continued_fraction(a: f64, b: f64, x: f64) -> f64 {
    const TINY: f64 = 1e-30;
    const EPSILON: f64 = 3e-16;
    const MAX_ITERATIONS: usize = 300;

    let qab = a + b;
    let qap = a + 1.0;
    let qam = a - 1.0;

    let mut c = 1.0;
    let mut d = 1.0 - qab * x / qap;
    if libm::fabs(d) < TINY {
        d = TINY;
    }
    d = 1.0 / d;
    let mut h = d;

    for m in 1..=MAX_ITERATIONS {
        let m_f = m as f64;
        let m2 = 2.0 * m_f;

        // The even step.
        let numerator = m_f * (b - m_f) * x / ((qam + m2) * (a + m2));
        d = 1.0 + numerator * d;
        if libm::fabs(d) < TINY {
            d = TINY;
        }
        c = 1.0 + numerator / c;
        if libm::fabs(c) < TINY {
            c = TINY;
        }
        d = 1.0 / d;
        h *= d * c;

        // The odd step.
        let numerator = -(a + m_f) * (qab + m_f) * x / ((a + m2) * (qap + m2));
        d = 1.0 + numerator * d;
        if libm::fabs(d) < TINY {
            d = TINY;
        }
        c = 1.0 + numerator / c;
        if libm::fabs(c) < TINY {
            c = TINY;
        }
        d = 1.0 / d;
        let step = d * c;
        h *= step;

        if libm::fabs(step - 1.0) < EPSILON {
            return h;
        }
    }
    // Did not converge in 300 iterations. `NaN` rather than the last iterate,
    // so a caller cannot mistake a failure for an answer.
    f64::NAN
}
