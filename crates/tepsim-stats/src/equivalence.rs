//! Welch's t-test, and the two one-sided tests built on it.
//!
//! # Why equivalence and not difference
//!
//! `PLAN.org` is explicit: "a failure to reject the null of no difference is
//! not evidence of equivalence". A t-test that fails to reject tells you the
//! data were too noisy or too few to see a difference, which is exactly what a
//! bad port with a small sample produces. TOST inverts the burden: the null is
//! that the difference is *at least* the equivalence margin, and rejecting it
//! in both directions is positive evidence that the difference is smaller than
//! the margin.
//!
//! Schuirmann, D. J. (1987), "A comparison of the two one-sided tests
//! procedure and the power approach for assessing the equivalence of average
//! bioavailability", *Journal of Pharmacokinetics and Biopharmaceutics* 15(6),
//! 657-680.

use crate::distribution::{student_t_quantile, student_t_two_sided_p};
use crate::special::sqrt;
use crate::summary::Summary;

/// Welch's unequal-variances t-test.
///
/// Every field is reported, not only the p-value: `CLAUDE.md`'s "record
/// numbers, not verdicts" applies to a test result as much as to a log entry.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WelchT {
    /// `mean(a) - mean(b)`.
    pub difference: f64,
    /// The standard error of that difference.
    pub standard_error: f64,
    /// The t statistic.
    pub t: f64,
    /// Welch-Satterthwaite degrees of freedom, generally not an integer.
    pub df: f64,
    /// `P(|T| >= |t|)`.
    pub p: f64,
}

/// Welch's t-test for the difference of two means.
///
/// Welch, B. L. (1947), "The generalization of Student's problem when several
/// different population variances are involved", *Biometrika* 34(1-2), 28-35.
///
/// ```text
/// se = sqrt(sa^2/na + sb^2/nb)
/// t  = (mean_a - mean_b) / se
/// df = se^4 / ( (sa^2/na)^2/(na-1) + (sb^2/nb)^2/(nb-1) )
/// ```
///
/// Unlike the pooled test this makes no equal-variance assumption, which
/// matters here: the failure mode Tier 5 exists to catch is the port's variance
/// collapsing relative to the Fortran's, and a test that assumes them equal is
/// the wrong instrument for it.
#[must_use]
pub fn welch_t(a: &Summary, b: &Summary) -> WelchT {
    let va = a.variance() / a.n() as f64;
    let vb = b.variance() / b.n() as f64;
    let standard_error = sqrt(va + vb);
    let difference = a.mean() - b.mean();
    let t = difference / standard_error;
    let df = (va + vb) * (va + vb) / (va * va / (a.n() - 1) as f64 + vb * vb / (b.n() - 1) as f64);
    WelchT {
        difference,
        standard_error,
        t,
        df,
        p: student_t_two_sided_p(t, df),
    }
}

/// The result of two one-sided tests.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Tost {
    /// The underlying Welch test, for its difference, standard error and df.
    pub welch: WelchT,
    /// The equivalence margin, in the units of the measurement.
    pub margin: f64,
    /// The p-value of `H0: difference <= -margin`.
    pub p_lower: f64,
    /// The p-value of `H0: difference >= +margin`.
    pub p_upper: f64,
    /// `max(p_lower, p_upper)`: the TOST p-value.
    pub p: f64,
    /// The `1 - 2*alpha` confidence interval for the difference.
    ///
    /// Not `1 - alpha`: the interval that corresponds to TOST at level `alpha`
    /// is the `1 - 2*alpha` one, and equivalence is declared exactly when it
    /// lies inside the margin. Reporting the wrong interval beside a correct
    /// verdict is a good way to make both untrustworthy.
    pub interval: (f64, f64),
    /// Whether equivalence was declared at `alpha`.
    pub equivalent: bool,
}

/// Two one-sided tests for the equivalence of two means.
///
/// The null is that the difference lies *outside* `[-margin, +margin]`.
/// Rejecting it from both sides at level `alpha` declares equivalence.
///
/// `PLAN.org` sets the margin for Tier 5 at one tenth of the Fortran's standard
/// deviation for that measurement. That choice is the caller's; this function
/// only applies it.
#[must_use]
pub fn tost(a: &Summary, b: &Summary, margin: f64, alpha: f64) -> Tost {
    let welch = welch_t(a, b);
    let se = welch.standard_error;
    let df = welch.df;

    // H0: difference <= -margin, rejected by a large positive statistic.
    let t_lower = (welch.difference + margin) / se;
    let p_lower = 0.5 * student_t_two_sided_p(t_lower, df);
    let p_lower = if t_lower > 0.0 {
        p_lower
    } else {
        1.0 - p_lower
    };

    // H0: difference >= +margin, rejected by a large negative statistic.
    let t_upper = (welch.difference - margin) / se;
    let p_upper = 0.5 * student_t_two_sided_p(t_upper, df);
    let p_upper = if t_upper < 0.0 {
        p_upper
    } else {
        1.0 - p_upper
    };

    let p = if p_lower > p_upper { p_lower } else { p_upper };
    let half_width = student_t_quantile(1.0 - alpha, df) * se;
    let interval = (welch.difference - half_width, welch.difference + half_width);

    Tost {
        welch,
        margin,
        p_lower,
        p_upper,
        p,
        interval,
        // Equivalently `p < alpha`; stated through the interval because that
        // is the form a reader can check against the reported numbers.
        equivalent: interval.0 > -margin && interval.1 < margin,
        // `interval` and `p < alpha` agree by construction; see the test
        // `the_interval_and_the_p_value_agree`.
    }
}

impl core::fmt::Display for Tost {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(
            f,
            "diff={:+.6e} margin=+-{:.6e} CI=[{:+.6e}, {:+.6e}] p={:.4} df={:.2} {}",
            self.welch.difference,
            self.margin,
            self.interval.0,
            self.interval.1,
            self.p,
            self.welch.df,
            if self.equivalent {
                "EQUIVALENT"
            } else {
                "not shown equivalent"
            }
        )
    }
}
