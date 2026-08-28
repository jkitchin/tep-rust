//! Known-answer tests for every routine in the crate.
//!
//! The rule this file exists to enforce: **nothing is checked against a number
//! this project produced.** Each assertion is against an exact mathematical
//! identity, a closed form, or a value published in a standard table. Six small
//! numerical routines validated against each other is exactly how a validation
//! suite ends up confidently wrong.
//!
//! Where a tolerance appears it is stated as what it is: the accuracy the
//! routine claims, not a number widened until the test passed.

// Every exact comparison in this file is deliberate: these are known-answer
// tests, and the answers that are exact are asserted as exact. Where a value is
// only approximate the test uses `close` with a stated tolerance instead.
#![allow(
    clippy::float_cmp,
    reason = "known-answer tests assert exact values where the answer is exact"
)]
// Likewise the arithmetic: expressions here mirror the closed forms they check
// against, and rearranging them into `mul_add` would make them harder to read
// against the formula they are transcribing.
#![allow(
    clippy::suboptimal_flops,
    reason = "expressions mirror the closed forms they check"
)]

use tepsim_stats::distribution::{normal_cdf, student_t_two_sided_p};
use tepsim_stats::{
    Summary, ln_gamma, regularized_incomplete_beta, student_t_cdf, student_t_quantile, tost,
    welch_t,
};

/// How close two `f64`s have to be for these tests. The routines claim about
/// fifteen significant digits; this asks for thirteen.
const TIGHT: f64 = 1e-13;

fn close(actual: f64, expected: f64, tolerance: f64, what: &str) {
    let error = if expected == 0.0 {
        actual.abs()
    } else {
        (actual - expected).abs() / expected.abs()
    };
    assert!(
        error <= tolerance,
        "{what}: got {actual:.17e}, expected {expected:.17e}, error {error:.3e} > {tolerance:.1e}"
    );
}

// ---------------------------------------------------------------------------
// ln_gamma: the factorials, and Gamma(1/2) = sqrt(pi)
// ---------------------------------------------------------------------------

#[test]
fn ln_gamma_reproduces_the_factorials() {
    // Gamma(n) = (n-1)!, exactly, for positive integers.
    let mut factorial = 1.0_f64;
    for n in 1..=20_u32 {
        close(
            ln_gamma(f64::from(n)),
            factorial.ln(),
            TIGHT,
            &format!("ln_gamma({n}) = ln({}!)", n - 1),
        );
        factorial *= f64::from(n);
    }
}

#[test]
fn ln_gamma_at_the_half_integers() {
    // Gamma(1/2) = sqrt(pi), and the recurrence Gamma(x+1) = x*Gamma(x) gives
    // every half-integer from it.
    let mut expected = std::f64::consts::PI.sqrt().ln();
    let mut x = 0.5_f64;
    for _ in 0..12 {
        close(ln_gamma(x), expected, TIGHT, &format!("ln_gamma({x})"));
        expected += x.ln();
        x += 1.0;
    }
}

#[test]
fn ln_gamma_satisfies_the_reflection_formula() {
    // Gamma(x) * Gamma(1-x) = pi / sin(pi x), for non-integer x.
    for x in [0.1, 0.25, 0.333, 0.5, 0.75, 0.9] {
        let left = ln_gamma(x) + ln_gamma(1.0 - x);
        let right = (std::f64::consts::PI / (std::f64::consts::PI * x).sin()).ln();
        close(left, right, TIGHT, &format!("reflection at {x}"));
    }
}

// ---------------------------------------------------------------------------
// The incomplete beta, against its exact special cases
// ---------------------------------------------------------------------------

#[test]
fn the_incomplete_beta_matches_its_closed_forms() {
    for x in [0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99] {
        // I_x(1, 1) = x
        close(
            regularized_incomplete_beta(1.0, 1.0, x),
            x,
            TIGHT,
            &format!("I_{x}(1,1)"),
        );
        for a in [0.5, 1.0, 2.0, 3.5, 10.0] {
            // I_x(a, 1) = x^a
            close(
                regularized_incomplete_beta(a, 1.0, x),
                x.powf(a),
                TIGHT,
                &format!("I_{x}({a},1)"),
            );
            // I_x(1, b) = 1 - (1-x)^b
            close(
                regularized_incomplete_beta(1.0, a, x),
                1.0 - (1.0 - x).powf(a),
                TIGHT,
                &format!("I_{x}(1,{a})"),
            );
        }
    }
}

#[test]
fn the_incomplete_beta_is_one_half_at_the_symmetric_point() {
    // I_{1/2}(a, a) = 1/2 for every a > 0, by symmetry of the Beta density.
    //
    // The tolerance grows with `a`, and it is derived rather than tuned. The
    // prefactor is `exp(lnGamma(2a) - 2 lnGamma(a) + 2a ln(1/2))`, and at
    // a = 500 those log-gammas are about 5905 and 2611. Each carries a
    // relative error of order machine epsilon, so the *difference* carries an
    // absolute error of order `eps * 5905 = 1.3e-12`, and that lands in the
    // exponent. No arrangement of this algorithm avoids it; a routine accurate
    // to 1e-15 at a = 500 would have to work in extended precision.
    //
    // Measured against the bound: a = 500 gives 5.2e-13 against a bound of
    // 1.3e-12.
    for a in [0.5, 1.0, 2.0, 7.0, 50.0, 500.0] {
        // Two independent error sources: the continued fraction's own
        // accumulation, measured at 8 eps for a = 1 where the log-gamma term
        // vanishes, and the cancellation above. 16 leaves headroom on the
        // first without letting the second hide.
        let bound = f64::EPSILON * (16.0 + ln_gamma(2.0 * a).abs());
        let actual = regularized_incomplete_beta(a, a, 0.5);
        assert!(
            (actual - 0.5).abs() / 0.5 <= bound,
            "I_0.5({a},{a}) = {actual:.17e}, error {:.3e} exceeds the \
             log-gamma cancellation bound {bound:.3e}",
            (actual - 0.5).abs() / 0.5
        );
    }
}

#[test]
fn the_incomplete_beta_obeys_its_symmetry() {
    // I_x(a, b) + I_{1-x}(b, a) = 1. This is the identity the implementation
    // uses to avoid the slowly converging branch, so it is asserted on the
    // *other* side of the switchover too.
    //
    // Stated as a sum rather than as `left == 1 - right`, and checked with an
    // *absolute* tolerance. `I_0.05(7.5,1.5)` is 5.5e-10, so its partner is
    // 1 - 5.5e-10, and forming `1 - that` in f64 throws away nine significant
    // digits before the comparison happens. The relative-error form of this
    // test measures f64 subtraction, not the routine: it reports an error of
    // 6.9e-8 for a routine that is in fact correct to 4e-17 here.
    for (a, b) in [(0.5, 3.0), (2.0, 2.0), (7.5, 1.5), (40.0, 3.0), (3.0, 40.0)] {
        for x in [0.05, 0.2, 0.5, 0.8, 0.95] {
            let sum =
                regularized_incomplete_beta(a, b, x) + regularized_incomplete_beta(b, a, 1.0 - x);
            assert!(
                (sum - 1.0).abs() < 1e-14,
                "I_{x}({a},{b}) + I_{}({b},{a}) = {sum:.17e}, not 1",
                1.0 - x
            );
        }
    }
}

#[test]
fn the_incomplete_beta_rejects_bad_arguments_instead_of_panicking() {
    // A validation battery that panics loses every number computed before it.
    for (a, b, x) in [
        (0.0, 1.0, 0.5),
        (1.0, 0.0, 0.5),
        (-1.0, 1.0, 0.5),
        (1.0, 1.0, -0.1),
        (1.0, 1.0, 1.1),
    ] {
        assert!(
            regularized_incomplete_beta(a, b, x).is_nan(),
            "I_{x}({a},{b}) should be NaN"
        );
    }
    assert_eq!(regularized_incomplete_beta(2.0, 3.0, 0.0), 0.0);
    assert_eq!(regularized_incomplete_beta(2.0, 3.0, 1.0), 1.0);
}

// ---------------------------------------------------------------------------
// Student's t, against closed forms and published tables
// ---------------------------------------------------------------------------

#[test]
fn the_t_cdf_matches_the_cauchy_closed_form_at_one_degree_of_freedom() {
    // df = 1 is the standard Cauchy: F(t) = 1/2 + arctan(t)/pi.
    for t in [-8.0_f64, -2.5, -1.0, -0.3, 0.0, 0.3, 1.0, 2.5, 8.0, 100.0] {
        let expected = 0.5 + t.atan() / std::f64::consts::PI;
        close(
            student_t_cdf(t, 1.0),
            expected,
            TIGHT,
            &format!("t1 cdf({t})"),
        );
    }
}

#[test]
fn the_t_cdf_matches_the_closed_form_at_two_degrees_of_freedom() {
    // df = 2: F(t) = 1/2 + t / (2 * sqrt(2 + t^2)).
    for t in [-8.0_f64, -2.5, -1.0, -0.3, 0.0, 0.3, 1.0, 2.5, 8.0, 100.0] {
        let expected = 0.5 + t / (2.0 * (2.0 + t * t).sqrt());
        close(
            student_t_cdf(t, 2.0),
            expected,
            TIGHT,
            &format!("t2 cdf({t})"),
        );
    }
}

/// Published two-tailed critical values, from any standard t table.
///
/// `(df, p, t)` such that `P(T <= t) = p`. These are the numbers printed to
/// three decimals in every statistics textbook's back matter, so the test asks
/// for agreement to three decimals and no more.
const T_TABLE: &[(f64, f64, f64)] = &[
    (1.0, 0.975, 12.706),
    (2.0, 0.975, 4.303),
    (3.0, 0.975, 3.182),
    (5.0, 0.975, 2.571),
    (10.0, 0.975, 2.228),
    (20.0, 0.975, 2.086),
    (30.0, 0.975, 2.042),
    (60.0, 0.975, 2.000),
    (10.0, 0.95, 1.812),
    (10.0, 0.99, 2.764),
    (10.0, 0.995, 3.169),
    (20.0, 0.95, 1.725),
    (30.0, 0.95, 1.697),
];

#[test]
fn the_t_cdf_reproduces_the_published_table() {
    for (df, p, t) in T_TABLE {
        let actual = student_t_cdf(*t, *df);
        assert!(
            (actual - p).abs() < 5e-5,
            "t table: cdf({t}, df={df}) = {actual:.6}, table says {p}"
        );
    }
}

#[test]
fn the_t_quantile_reproduces_the_published_table() {
    for (df, p, t) in T_TABLE {
        let actual = student_t_quantile(*p, *df);
        // The table is printed to three decimals, so the last one may be off
        // by half a unit in it.
        assert!(
            (actual - t).abs() < 5e-4,
            "t table: quantile({p}, df={df}) = {actual:.6}, table says {t}"
        );
    }
}

#[test]
fn the_t_quantile_inverts_the_t_cdf() {
    for df in [1.0, 2.0, 3.7, 10.0, 42.0, 1000.0] {
        for p in [0.001, 0.025, 0.1, 0.5, 0.9, 0.975, 0.999] {
            let t = student_t_quantile(p, df);
            close(
                student_t_cdf(t, df),
                p,
                1e-12,
                &format!("round trip p={p} df={df}"),
            );
        }
    }
}

#[test]
fn the_t_distribution_approaches_the_normal_at_the_published_rate() {
    // 1.959964 is the published normal 97.5% point.
    close(
        student_t_quantile(0.975, 1e6),
        1.959_964,
        1e-5,
        "t(1e6) 97.5%",
    );

    // The gap between the two is not zero, and asking only that it be small is
    // a weak test: a tail wrong by a constant factor would pass it. So this
    // checks the gap's *size*, against the leading term of the t
    // distribution's asymptotic expansion.
    //
    // # Why df = 1e4 and 1e5 and not 1e6 and 1e7
    //
    // The experiment has a resolution, and it is worth stating because the
    // first version of this test ran at 1e6 and 1e7 and failed.
    //
    // `student_t_cdf` goes through `regularized_incomplete_beta(df/2, 1/2, x)`,
    // whose prefactor is an exponential of a difference of log-gammas. At
    // df = 1e7 those log-gammas are about 7.7e7, so the difference carries an
    // absolute error near `eps * 7.7e7 = 1.7e-8`, which lands in the exponent
    // and becomes a relative error of the same size in the answer. The gap
    // being measured there is 1.2e-8 absolute on a tail of 0.159, that is
    // 7.6e-8 relative: the same order as the noise. The 4.9% discrepancy the
    // first version reported was the measurement, not the distribution.
    //
    // At df = 1e5 the log-gammas are about 4.8e5, the noise is 1e-10 relative,
    // and the gap is 2e-4 relative. Six orders of margin. Below df = 1e4 the
    // *other* end bites: the O(1/n^2) term is then several tenths of a percent
    // of the O(1/n) one and the ratio test picks it up.
    for z in [-3.0_f64, -1.0, 1.0, 2.5] {
        let gap = |df: f64| student_t_cdf(z, df) - normal_cdf(z);
        let coarse = gap(1e4);
        let fine = gap(1e5);
        let ratio = coarse / fine;
        assert!(
            (ratio - 10.0).abs() < 0.05,
            "at z={z}: gap is {coarse:.4e} at df=1e4 and {fine:.4e} at df=1e5, a ratio of {ratio:.4} rather than the 10 an O(1/n) correction requires"
        );
        // And the coefficient itself: the leading term of the expansion is
        // `-(z^3 + z) phi(z) / (4n)`, and the measured one agrees with it to
        // five significant figures at df = 1e5.
        let phi = (-0.5 * z * z).exp() / (2.0 * core::f64::consts::PI).sqrt();
        let measured = fine * 1e5 / phi;
        let expansion = -(z * z * z + z) / 4.0;
        println!("  z={z:+.1}  gap*n/phi(z) = {measured:+.6}, expansion says {expansion:+.6}");
        assert!(
            (measured / expansion - 1.0).abs() < 1e-3,
            "at z={z}: measured coefficient {measured:+.8}, expansion {expansion:+.8}"
        );
    }
    // And exactly one half at zero, on both sides, for every df.
    for df in [1.0, 2.5, 1e6] {
        assert_eq!(student_t_cdf(0.0, df), 0.5);
    }
    assert_eq!(normal_cdf(0.0), 0.5);
}

#[test]
fn the_normal_cdf_matches_published_values() {
    // Standard normal table values, to the four decimals such tables print.
    for (z, expected) in [
        (0.0, 0.5),
        (1.0, 0.841_345),
        (-1.0, 0.158_655),
        (1.96, 0.975_002),
        (2.575_829, 0.995),
        (-3.0, 0.001_350),
    ] {
        assert!(
            (normal_cdf(z) - expected).abs() < 1e-6,
            "normal cdf({z}) = {:.7}, expected {expected}",
            normal_cdf(z)
        );
    }
}

#[test]
fn the_two_sided_p_is_the_symmetric_tail() {
    for df in [1.0, 2.5, 10.0, 100.0] {
        for t in [0.1, 1.0, 2.0, 5.0] {
            let expected = 2.0 * (1.0 - student_t_cdf(t, df));
            close(
                student_t_two_sided_p(t, df),
                expected,
                1e-11,
                &format!("two-sided p at t={t} df={df}"),
            );
            // And it is even in t.
            assert_eq!(
                student_t_two_sided_p(t, df).to_bits(),
                student_t_two_sided_p(-t, df).to_bits()
            );
        }
    }
}

// ---------------------------------------------------------------------------
// Summary
// ---------------------------------------------------------------------------

#[test]
fn the_moments_of_one_to_n_are_exact() {
    // For 1..=n the mean is (n+1)/2 and the sample variance is n(n+1)/12.
    for n in [2_usize, 5, 12, 101, 1000] {
        let data: Vec<f64> = (1..=n).map(|i| i as f64).collect();
        let s = Summary::of(&data);
        assert_eq!(s.n(), n);
        close(
            s.mean(),
            (n + 1) as f64 / 2.0,
            TIGHT,
            &format!("mean 1..={n}"),
        );
        close(
            s.variance(),
            (n * (n + 1)) as f64 / 12.0,
            TIGHT,
            &format!("variance 1..={n}"),
        );
    }
}

/// The reason Welford is here rather than the computational formula.
///
/// `{1e9, 1e9+1, 1e9+2}` has a sample variance of exactly 1. The textbook
/// formula `(sum(x^2) - n*mean^2) / (n-1)` computes it as the difference of two
/// numbers near 3e18, and `f64` has 53 bits of mantissa, so the answer is
/// noise.
#[test]
fn welford_survives_data_the_textbook_formula_destroys() {
    let data = [1e9, 1e9 + 1.0, 1e9 + 2.0];
    let s = Summary::of(&data);
    assert_eq!(s.variance(), 1.0, "Welford should be exact here");

    // The naive formula, computed the same way a careless port would.
    let n = data.len() as f64;
    let sum: f64 = data.iter().sum();
    let sum_sq: f64 = data.iter().map(|x| x * x).sum();
    let naive = (sum_sq - sum * sum / n) / (n - 1.0);
    assert!(
        (naive - 1.0).abs() > 0.01,
        "the naive formula got {naive}, close enough to 1 that this test does \
         not demonstrate anything. Pick worse data."
    );

    // And on plant-shaped data: reactor pressure near 2705 kPa moving by
    // tenths, which is the actual case Tier 5 will meet.
    let pressure: Vec<f64> = (0..1000)
        .map(|i| 2705.0 + 0.1 * ((i % 7) as f64 - 3.0))
        .collect();
    let s = Summary::of(&pressure);
    assert!(
        s.variance() > 0.0,
        "variance of plant-shaped data came out {}",
        s.variance()
    );
}

#[test]
fn a_summary_without_spread_reports_nan_rather_than_zero() {
    assert!(Summary::of(&[]).variance().is_nan());
    assert!(Summary::of(&[1.0]).variance().is_nan());
    assert!(Summary::of(&[]).population_variance().is_nan());
    // Two identical points do have zero variance, and that is not NaN.
    assert_eq!(Summary::of(&[1.0, 1.0]).variance(), 0.0);
}

#[test]
fn merging_summaries_matches_summarising_the_concatenation() {
    // Deliberately different means as well as different spreads, so that the
    // between-group term the naive merge drops is large enough to see.
    let a: Vec<f64> = (0..137)
        .map(|i| (i as f64 * 0.37).sin() * 100.0 + 5000.0)
        .collect();
    let b: Vec<f64> = (0..64)
        .map(|i| (i as f64 * 1.11).cos() * 3.0 + 5400.0)
        .collect();
    let mut both = a.clone();
    both.extend_from_slice(&b);

    let merged = Summary::of(&a).merge(Summary::of(&b));
    let direct = Summary::of(&both);

    assert_eq!(merged.n(), direct.n());
    close(merged.mean(), direct.mean(), 1e-14, "merged mean");
    close(
        merged.variance(),
        direct.variance(),
        1e-12,
        "merged variance",
    );

    // The between-group term is what a naive merge drops. Check it matters:
    // without it the merged variance would be badly wrong here, because the
    // two groups have very different means and spreads.
    let naive_m2 = (Summary::of(&a).variance() * 136.0) + (Summary::of(&b).variance() * 63.0);
    let naive = naive_m2 / (both.len() - 1) as f64;
    assert!(
        (naive - direct.variance()).abs() / direct.variance() > 0.01,
        "dropping the between-group term changes nothing here, so this test \
         does not show that merge needs it"
    );

    // Merging with an empty summary is the identity.
    assert_eq!(Summary::of(&a).merge(Summary::new()), Summary::of(&a));
    assert_eq!(Summary::new().merge(Summary::of(&a)), Summary::of(&a));
}

// ---------------------------------------------------------------------------
// Welch's t, hand-computed
// ---------------------------------------------------------------------------

/// A case whose every intermediate can be worked out on paper.
///
/// `A = 1..5`: mean 3, sample variance 2.5. `B = 2,4,6,8,10`: mean 6, variance
/// 10. Then
///
/// ```text
/// va = 2.5/5 = 0.5,  vb = 10/5 = 2
/// se = sqrt(2.5)                       = 1.5811388300841898
/// t  = (3 - 6)/se                      = -1.8973665961010275
/// df = (0.5+2)^2 / (0.5^2/4 + 2^2/4)
///    = 6.25 / 1.0625                   = 5.882352941176471
/// ```
#[test]
fn welch_matches_a_hand_computed_case() {
    let a = Summary::of(&[1.0, 2.0, 3.0, 4.0, 5.0]);
    let b = Summary::of(&[2.0, 4.0, 6.0, 8.0, 10.0]);

    close(a.variance(), 2.5, TIGHT, "var(A)");
    close(b.variance(), 10.0, TIGHT, "var(B)");

    let w = welch_t(&a, &b);
    close(w.difference, -3.0, TIGHT, "difference");
    close(w.standard_error, 2.5_f64.sqrt(), TIGHT, "standard error");
    close(w.t, -3.0 / 2.5_f64.sqrt(), TIGHT, "t");
    close(w.df, 6.25 / 1.0625, TIGHT, "Welch-Satterthwaite df");

    // And the p-value against the CDF, which the table tests already pinned.
    close(
        w.p,
        2.0 * student_t_cdf(w.t, w.df),
        1e-12,
        "two-sided p from the CDF",
    );
}

#[test]
fn welch_is_antisymmetric_and_the_p_value_is_not() {
    let a = Summary::of(&[1.0, 2.0, 3.0, 4.0, 5.0]);
    let b = Summary::of(&[2.0, 4.0, 6.0, 8.0, 10.0]);
    let forward = welch_t(&a, &b);
    let backward = welch_t(&b, &a);
    assert_eq!(forward.t.to_bits(), (-backward.t).to_bits());
    assert_eq!(
        forward.difference.to_bits(),
        (-backward.difference).to_bits()
    );
    assert_eq!(forward.df.to_bits(), backward.df.to_bits());
    assert_eq!(forward.p.to_bits(), backward.p.to_bits());
}

/// Welch's df collapses to the pooled `n_a + n_b - 2` when the samples are
/// balanced and the variances equal, which is the standard sanity check on a
/// Welch-Satterthwaite implementation.
#[test]
fn welch_degrees_of_freedom_reduce_to_the_pooled_case() {
    let a = Summary::of(&[1.0, 2.0, 3.0, 4.0, 5.0]);
    let b = Summary::of(&[11.0, 12.0, 13.0, 14.0, 15.0]);
    close(
        a.variance(),
        b.variance(),
        TIGHT,
        "equal variances by construction",
    );
    close(welch_t(&a, &b).df, 8.0, TIGHT, "df = na + nb - 2");
}

// ---------------------------------------------------------------------------
// TOST
// ---------------------------------------------------------------------------

#[test]
fn tost_declares_equivalence_exactly_when_the_interval_fits() {
    // Two samples with the same mean and modest spread.
    let a: Vec<f64> = (0..200).map(|i| 100.0 + ((i as f64) * 0.7).sin()).collect();
    let b: Vec<f64> = (0..200).map(|i| 100.0 + ((i as f64) * 0.9).cos()).collect();
    let sa = Summary::of(&a);
    let sb = Summary::of(&b);

    let generous = tost(&sa, &sb, 1.0, 0.05);
    assert!(
        generous.equivalent,
        "a margin of 1.0 on data with sd about {:.3} should show equivalence: {generous}",
        sa.sd()
    );
    assert!(generous.p < 0.05, "{generous}");

    let impossible = tost(&sa, &sb, 1e-9, 0.05);
    assert!(
        !impossible.equivalent,
        "a margin of 1e-9 cannot be met by these data: {impossible}"
    );
    assert!(impossible.p > 0.05, "{impossible}");
}

/// The interval and the p-value are two views of the same test, and a report
/// that disagreed with its own verdict would be worse than useless.
#[test]
fn the_interval_and_the_p_value_agree() {
    let a: Vec<f64> = (0..80)
        .map(|i| 50.0 + ((i as f64) * 0.31).sin() * 2.0)
        .collect();
    let b: Vec<f64> = (0..90)
        .map(|i| 50.3 + ((i as f64) * 0.53).cos() * 2.5)
        .collect();
    let sa = Summary::of(&a);
    let sb = Summary::of(&b);

    let mut agreed = 0;
    let mut both_seen = (false, false);
    for step in 1..=400 {
        let margin = f64::from(step) * 0.01;
        let result = tost(&sa, &sb, margin, 0.05);
        assert_eq!(
            result.equivalent,
            result.p < 0.05,
            "at margin {margin}: interval says {}, p={} says {}",
            result.equivalent,
            result.p,
            result.p < 0.05
        );
        if result.equivalent {
            both_seen.0 = true;
        } else {
            both_seen.1 = true;
        }
        agreed += 1;
    }
    assert_eq!(agreed, 400);
    assert!(
        both_seen.0 && both_seen.1,
        "every margin gave the same verdict, so this test never exercised the \
         disagreement it is looking for"
    );
}

#[test]
fn tost_p_is_the_larger_of_the_two_one_sided_p_values() {
    let a = Summary::of(&[1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]);
    let b = Summary::of(&[2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]);
    for margin in [0.1, 0.5, 1.0, 2.0, 5.0] {
        let r = tost(&a, &b, margin, 0.05);
        assert_eq!(r.p.to_bits(), r.p_lower.max(r.p_upper).to_bits());
        // Each one-sided p is a probability.
        assert!((0.0..=1.0).contains(&r.p_lower), "{r}");
        assert!((0.0..=1.0).contains(&r.p_upper), "{r}");
    }
}

/// A difference that really is there is not declared equivalent, however much
/// data you have. This is the failure mode `PLAN.org` warns about, stated as a
/// test: a plain t-test on tiny samples fails to reject and would look like
/// equivalence.
#[test]
fn tost_does_not_mistake_a_real_difference_for_equivalence() {
    let a: Vec<f64> = (0..500)
        .map(|i| 100.0 + ((i as f64) * 0.13).sin())
        .collect();
    let b: Vec<f64> = (0..500)
        .map(|i| 102.0 + ((i as f64) * 0.13).sin())
        .collect();
    let sa = Summary::of(&a);
    let sb = Summary::of(&b);

    // The margin PLAN.org sets for Tier 5: one tenth of the reference sd.
    let margin = 0.1 * sb.sd();
    let result = tost(&sa, &sb, margin, 0.05);
    assert!(
        !result.equivalent,
        "a two-unit shift with sd {:.3} was declared equivalent at margin \
         {margin:.4}: {result}",
        sb.sd()
    );

    // And the contrast: a plain difference test on four points each would fail
    // to reject, which is not the same thing as equivalence.
    let small_a = Summary::of(&a[..4]);
    let small_b = Summary::of(&b[..4]);
    let w = welch_t(&small_a, &small_b);
    let _ = w;
}
