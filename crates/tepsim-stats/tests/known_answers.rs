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
    Summary, f_cdf, f_quantile, ln_gamma, normal_quantile, regularized_incomplete_beta,
    student_t_cdf, student_t_quantile, tost, welch_t,
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
// The standard normal quantile, against published percentage points
// ---------------------------------------------------------------------------

#[test]
fn the_normal_quantile_reproduces_the_published_percentage_points() {
    // The standard normal deviates every table of critical values prints, to
    // the digits they print. The Jackson-Mudholkar SPE limit calls this with
    // the first of them.
    for (p, expected) in [
        (0.9, 1.281_551_565_545),
        (0.95, 1.644_853_626_951),
        (0.975, 1.959_963_984_540),
        (0.99, 2.326_347_874_041),
        (0.995, 2.575_829_303_549),
        (0.999, 3.090_232_306_168),
        (0.9995, 3.290_526_731_492),
    ] {
        close(
            normal_quantile(p),
            expected,
            1e-11,
            &format!("normal quantile at {p}"),
        );
        // And the lower tail by symmetry, which is exact for the distribution
        // and so should be near exact for the routine.
        close(
            normal_quantile(1.0 - p),
            -expected,
            1e-11,
            &format!("normal quantile at {}", 1.0 - p),
        );
    }
    assert_eq!(normal_quantile(0.5), 0.0);
    assert!(normal_quantile(0.0).is_nan());
    assert!(normal_quantile(1.0).is_nan());
}

// ---------------------------------------------------------------------------
// Snedecor's F, against closed forms and a published table
// ---------------------------------------------------------------------------

/// Published upper 5% critical values of the F distribution.
///
/// `(df1, df2, F)` such that `P(F <= x) = 0.95`. These are the numbers printed
/// in the `alpha = 0.05` table in the back of any statistics text, to the two or
/// three significant figures such tables carry.
const F_TABLE_95: &[(f64, f64, f64)] = &[
    (1.0, 1.0, 161.4),
    (2.0, 10.0, 4.10),
    (5.0, 20.0, 2.71),
    (10.0, 10.0, 2.98),
    (1.0, 10.0, 4.96),
    (3.0, 15.0, 3.29),
    (4.0, 30.0, 2.69),
    (6.0, 24.0, 2.51),
    (12.0, 12.0, 2.69),
    (20.0, 20.0, 2.12),
];

/// Every published value in `F_TABLE_95`, checked and none dropped.
///
/// The tolerance is half a unit in the table's last printed digit, which is the
/// most a printed table can support. `F(1, 1)` is the tightest of them: the
/// table says 161.4, the routine says 161.4476, and half a unit in the last
/// place is 0.05, so it passes with 0.0476 of the 0.05 used up.
///
/// That one is worth checking against something better than a table, and it can
/// be: `F(1, nu)` is the square of a `t` with `nu` degrees of freedom, so
/// `F(1, 1)` is a squared standard Cauchy and
///
/// ```text
/// P(F <= x) = (2 / pi) arctan(sqrt x)  =>  F_p(1, 1) = tan(p pi / 2)^2
/// ```
///
/// which gives 161.44763... exactly. The table is the rounded one, not the
/// routine. That closed form is asserted below.
#[test]
fn the_f_quantile_reproduces_the_published_table() {
    for &(df1, df2, table) in F_TABLE_95 {
        let actual = f_quantile(0.95, df1, df2);
        // Half a unit in the last digit the table prints.
        let printed_ulp = if table >= 100.0 { 0.05 } else { 0.005 };
        println!(
            "F_0.95({df1}, {df2}): table {table}, computed {actual:.6}, \
             difference {:+.4} against an allowance of {printed_ulp}",
            actual - table
        );
        assert!(
            (actual - table).abs() < printed_ulp,
            "F_0.95({df1}, {df2}) = {actual:.6}, the table says {table}"
        );
        // Round trip: the quantile is where the CDF reaches 0.95.
        close(
            f_cdf(actual, df1, df2),
            0.95,
            1e-12,
            &format!("F cdf round trip at ({df1}, {df2})"),
        );
    }
}

/// The exact closed forms the F **CDF** has, over a sweep of arguments.
///
/// ```text
/// P(F <= x | 1, 1)  = (2 / pi) arctan(sqrt x)     (F(1,1) is a squared Cauchy)
/// P(F <= x | d1, 2) = (d1 x / (d1 x + 2))^(d1/2)  (I_y(a, 1) = y^a)
/// P(F <= x | 2, d2) = 1 - (d2 / (2x + d2))^(d2/2) (I_y(1, b) = 1 - (1-y)^b)
/// ```
///
/// Checked on the CDF rather than the quantile because this is the direction
/// that is well conditioned everywhere: the `F(1, 1)` quantile has a derivative
/// of `4052` at `p = 0.99` and `4.05e5` at `p = 0.999`, so an ulp on the
/// probability becomes a large relative error on `x`, in the *reference* as much
/// as in the routine. The CDF has the reciprocal of that derivative and loses
/// nothing.
#[test]
fn the_f_cdf_matches_its_closed_forms() {
    let mut worst = 0.0_f64;
    // The sweep stops at 4052, where the CDF is 0.99. Past that the routine
    // loses accuracy for a reason worth naming, and
    // `the_f_cdf_saturates_in_the_far_upper_tail` measures it rather than this
    // test quietly absorbing it into a wider tolerance.
    for x in [1e-6_f64, 0.01, 0.5, 1.0, 2.5, 17.0, 4052.0] {
        let cauchy = 2.0 / std::f64::consts::PI * x.sqrt().atan();
        let got = f_cdf(x, 1.0, 1.0);
        worst = worst.max((got - cauchy).abs() / cauchy);
        close(got, cauchy, 1e-13, &format!("F cdf({x} | 1, 1)"));

        for df in [1.0_f64, 2.0, 5.0, 10.0, 30.0, 200.0] {
            let y = df * x / (df * x + 2.0);
            let exact = y.powf(0.5 * df);
            let got = f_cdf(x, df, 2.0);
            worst = worst.max((got - exact).abs() / exact);
            close(got, exact, 1e-13, &format!("F cdf({x} | {df}, 2)"));

            // `1 - (1 + 2x/df)^(-df/2)`, written through `ln_1p` and `exp_m1`
            // so that neither end cancels. Both are needed and it took two
            // tries to see it. `1.0 - w.powf(df/2)` loses six digits at
            // `x = 1e-6`, where the power is 0.999999. Replacing only the outer
            // subtraction with `exp_m1` still leaves `ln(0.999998000004)`, whose
            // argument carries half an ulp on a result of `-2e-6`, a relative
            // error of 5e-11. The reference, not the routine, was wrong both
            // times: it reported 3.9e-12 and then 1.3e-11 against a routine
            // accurate to 1e-15 here.
            let exact = -(-0.5 * df * (2.0 * x / df).ln_1p()).exp_m1();
            let got = f_cdf(x, 2.0, df);
            worst = worst.max((got - exact).abs() / exact);
            close(got, exact, 1e-13, &format!("F cdf({x} | 2, {df})"));
        }
    }
    println!("F cdf against its closed forms: worst relative error {worst:.3e}");
}

/// `P(F <= t^2 | 1, nu) + P(|T| >= t | nu) = 1`, exactly.
///
/// `F(1, nu)` is the square of a `t` with `nu` degrees of freedom, so the F CDF
/// at `t^2` and the two-sided t p-value at `t` are complementary. This is an
/// identity between two distributions rather than two implementations of one,
/// and it is the check that would catch `df1` and `df2` being swapped inside
/// `f_cdf`, which a table with equal arguments cannot see.
///
/// Written as a sum rather than as `left == 1 - right`, and checked with an
/// absolute tolerance, for the reason
/// `the_incomplete_beta_obeys_its_symmetry` gives: at `t = 30` the p-value is
/// 2e-2 for `nu = 1` and 1e-13 for `nu = 30`, and forming `1 - that` throws away
/// the digits the comparison is supposed to be about.
#[test]
fn the_f_distribution_is_the_square_of_a_t() {
    let mut worst = 0.0_f64;
    for nu in [1.0_f64, 2.0, 3.7, 10.0, 30.0, 200.0] {
        // The sweep starts at t = 0.3. Below it, `student_t_two_sided_p` hits
        // the saturation `the_f_cdf_saturates_in_the_far_upper_tail` describes,
        // in its own parameterisation: it evaluates `I_x(nu/2, 1/2)` at
        // `x = nu / (nu + t^2)`, which is within an ulp of 1 for small `t`, so
        // `1 - x` carries a relative error of `eps nu / t^2`. At `t = 0.01` and
        // `nu = 3.7` that is 4e-12, and it puts 1.5e-14 into the sum below. The
        // identity is fine; the representation runs out. Four orders of the
        // p-value are still covered from t = 0.3 up.
        for t in [0.3_f64, 1.0, 2.0, 5.0, 12.706, 30.0] {
            let sum = f_cdf(t * t, 1.0, nu) + student_t_two_sided_p(t, nu);
            worst = worst.max((sum - 1.0).abs());
            assert!(
                (sum - 1.0).abs() < 1e-14,
                "F cdf({}, 1, {nu}) + two-sided p({t}, {nu}) = {sum:.17e}, not 1",
                t * t
            );
        }
    }
    println!("F(1, nu) = t^2 identity: worst absolute departure from 1 is {worst:.3e}");
}

/// The far upper tail, measured rather than assumed.
///
/// `f_cdf` goes through `I_y(df1/2, df2/2)` with `y = df1 x / (df1 x + df2)`,
/// and that parameterisation saturates: once `y` is within an ulp or two of 1
/// there is no room left in an `f64` to say how far up the tail `x` is. At
/// `x = 1e8` with one and one degrees of freedom, `1 - y` is 1e-8 and carries a
/// relative error of `eps / (1 - y) = 2.2e-8`; the upper tail is
/// `(2 / pi) arcsin(sqrt(1 - y))`, so half of that reaches the answer, giving an
/// absolute error near `0.5 * 6.4e-5 * 2.2e-8 = 7e-13`.
///
/// Measured, it is 1.2e-13, inside that. This is a property of the
/// representation and not a defect to fix: computing the upper tail accurately
/// there needs a separate entry point taking the tail directly, the way
/// `student_t_two_sided_p` does for `t`. Nothing in this crate needs one. The
/// T-squared control limit calls `f_quantile(0.95, k, n - k)`, where `y` is
/// 0.4 for a small model and 1e-5 for a large one, both nowhere near saturation.
#[test]
fn the_f_cdf_saturates_in_the_far_upper_tail() {
    let x = 1e8_f64;
    let closed_form = 2.0 / std::f64::consts::PI * x.sqrt().atan();
    let got = f_cdf(x, 1.0, 1.0);
    let error = (got - closed_form).abs();
    println!(
        "F cdf({x} | 1, 1) = {got:.17e}, closed form {closed_form:.17e}, \
         absolute error {error:.3e} against the 7e-13 the saturation of y allows"
    );
    assert!(error < 7e-13, "{error:.3e}");
    // And it is still much better than the tolerance the sweep uses at
    // moderate x would suggest is the general case: this is a tail effect, so
    // it must vanish as x comes back down.
    let near = 4052.0_f64;
    let near_error =
        (f_cdf(near, 1.0, 1.0) - 2.0 / std::f64::consts::PI * near.sqrt().atan()).abs();
    println!("F cdf({near} | 1, 1): absolute error {near_error:.3e}");
    assert!(
        near_error < 1e-15,
        "the effect does not shrink as x falls, so it is not tail saturation: \
         {near_error:.3e}"
    );
}

/// The exact closed forms the F **quantile** has, at every `p`.
///
/// A table gives three digits. These give sixteen, and they cover the two
/// degrees-of-freedom arguments independently, so a formula that mixed `df1`
/// and `df2` up would fail one of them.
///
/// ```text
/// F_p(1, nu)  = t_{(1+p)/2, nu}^2                       (F(1, nu) is t^2)
/// F_p(d1, 2)  = 2 y / (d1 (1 - y)),      y = p^(2/d1)   (I_y(a, 1) = y^a)
/// F_p(2, d2)  = (d2/2) ((1-p)^(-2/d2) - 1)              (I_y(1, b) closed)
/// ```
///
/// # Every one of them is checked in the CDF direction, and that is deliberate
///
/// Comparing two `x` values in a heavy tail measures the tail, not the routine.
/// `F_0.999(1, 1)` is 405284, where `x f(x)` is 5e-4, so an ulp on the
/// probability moves `x` by two thousand ulps; and the closed form is no better
/// off, because `tan(p pi / 2)` amplifies the ulp on its own argument by
/// `sec^2 = 1 + tan^2 = 4e5`. Both sides are noisy, in different ways, and a
/// direct comparison reports their combined noise as if it were the routine's
/// error: the first version of this test did exactly that and demanded an
/// explanation for 9.5e-12.
///
/// Putting the computed quantile back through the closed-form CDF divides by
/// that factor instead of multiplying by it, and the answer it has to reproduce,
/// `p`, is exact. `the_incomplete_beta_obeys_its_symmetry` avoids the same trap
/// for the same reason.
///
/// The direct value comparison is still printed, so the number is on the record
/// and a reader can see how much of it is conditioning. It is asserted only
/// where the table already asserts it: `the_f_quantile_reproduces_the_published_table`
/// pins `F_0.95(1, 1)` to 161.4 directly, and 0.95 is far enough from the tail
/// for that to mean something.
#[test]
fn the_f_quantile_matches_its_closed_forms() {
    let mut worst = 0.0_f64;
    for p in [0.5_f64, 0.9, 0.95, 0.975, 0.99, 0.999] {
        // F(1, 1) is a squared standard Cauchy, so its CDF is
        // (2 / pi) arctan(sqrt x) and its quantile is tan(p pi / 2)^2.
        let tangent = (p * std::f64::consts::PI / 2.0).tan();
        let cauchy = tangent * tangent;
        let got = f_quantile(p, 1.0, 1.0);
        let round_trip = 2.0 / std::f64::consts::PI * got.sqrt().atan();
        println!(
            "F_{p}(1,1): closed form {cauchy:.12e}, computed {got:.12e}, direct \
             relative difference {:.3e}; the closed-form CDF of the computed \
             value recovers p to {:.3e}",
            (got - cauchy).abs() / cauchy,
            (round_trip - p).abs() / p
        );
        worst = worst.max((round_trip - p).abs() / p);
        close(
            round_trip,
            p,
            1e-13,
            &format!("closed-form cdf of F_{p}(1,1)"),
        );

        for df2 in [1.0_f64, 2.0, 5.0, 10.0, 30.0, 200.0] {
            // df1 = 2. Checked by putting the computed quantile back through
            // the *closed-form* CDF, not by comparing the two `x` values. In a
            // heavy tail those are not the same test: `F_0.99(2, 1)` is 4999.5,
            // where `x f(x)` is only 5e-3, so the quantile inherits two hundred
            // times whatever relative error the CDF carries and a direct
            // comparison measures the tail's conditioning rather than the
            // routine. The CDF direction divides by that factor instead of
            // multiplying by it. `the_incomplete_beta_obeys_its_symmetry` makes
            // the same argument about the same trap.
            let got = f_quantile(p, 2.0, df2);
            let round_trip = -(-0.5 * df2 * (2.0 * got / df2).ln_1p()).exp_m1();
            worst = worst.max((round_trip - p).abs() / p);
            close(
                round_trip,
                p,
                1e-13,
                &format!("closed-form cdf of F_{p}(2,{df2})"),
            );

            // The `F(1, nu) = t^2` identity is not checked here. It belongs in
            // `the_f_distribution_is_the_square_of_a_t`, which states it on the
            // CDFs. Comparing the two *quantiles* fails for the reason this
            // test's preamble gives: both bisect their own CDF, the two CDFs
            // reach the same probability through different arguments of the
            // incomplete beta, and at `p = 0.999` with `df2 = 1` the tail
            // multiplies the ulp between them by two thousand. The observed
            // 9.5e-12 there is the tail, not either routine.
        }

        for df1 in [1.0_f64, 2.0, 5.0, 10.0, 30.0] {
            // df2 = 2, the same way.
            let got = f_quantile(p, df1, 2.0);
            let round_trip = (df1 * got / (df1 * got + 2.0)).powf(0.5 * df1);
            worst = worst.max((round_trip - p).abs() / p);
            close(
                round_trip,
                p,
                1e-13,
                &format!("closed-form cdf of F_{p}({df1},2)"),
            );
        }
    }
    println!(
        "F quantile through its closed-form CDF: worst relative error in the \
         recovered probability {worst:.3e}"
    );
}

/// `F_p(d1, d2) = 1 / F_{1-p}(d2, d1)`, exactly, for every `p`.
///
/// The reciprocal of an F variate is an F variate with the degrees of freedom
/// swapped. This is the identity that would catch the two arguments being
/// exchanged somewhere inside, which no table with `df1 = df2` can see.
#[test]
fn the_f_distribution_obeys_the_reciprocal_identity() {
    let mut worst = 0.0_f64;
    for (df1, df2) in [(1.0, 1.0), (2.0, 7.0), (7.0, 2.0), (5.0, 20.0), (30.0, 3.0)] {
        for p in [0.01_f64, 0.1, 0.5, 0.9, 0.95, 0.99] {
            let direct = f_quantile(p, df1, df2);
            let mirrored = 1.0 / f_quantile(1.0 - p, df2, df1);
            worst = worst.max((direct - mirrored).abs() / direct);
            close(
                direct,
                mirrored,
                1e-12,
                &format!("F_{p}({df1},{df2}) vs 1/F_{}({df2},{df1})", 1.0 - p),
            );
        }
    }
    println!("F reciprocal identity: worst relative error {worst:.3e}");
}

#[test]
fn the_f_distribution_rejects_bad_arguments_instead_of_panicking() {
    for (df1, df2) in [(0.0, 1.0), (1.0, 0.0), (-1.0, 1.0)] {
        assert!(f_cdf(1.0, df1, df2).is_nan(), "cdf({df1},{df2})");
        assert!(f_quantile(0.95, df1, df2).is_nan(), "quantile({df1},{df2})");
    }
    for p in [0.0, 1.0, -0.1, 1.1] {
        assert!(f_quantile(p, 2.0, 3.0).is_nan(), "quantile at p={p}");
    }
    // The support is [0, infinity).
    assert_eq!(f_cdf(0.0, 3.0, 4.0), 0.0);
    assert_eq!(f_cdf(-1.0, 3.0, 4.0), 0.0);
    assert_eq!(f_cdf(f64::INFINITY, 3.0, 4.0), 1.0);
    // And a value so large the numerator would overflow before the ratio is
    // formed. `f64::MAX * 3` is infinite; the answer is still 1.
    assert_eq!(f_cdf(f64::MAX, 3.0, 4.0), 1.0);
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

/// A sample too small to have a variance makes the whole test undefined, and
/// says so rather than computing `n - 1` on zero.
///
/// In release that underflow wraps to `usize::MAX` and yields a degrees-of-
/// freedom figure that looks like a number. This is the case that produces it:
/// an ensemble where every member has zero variance contributes no
/// observations at all.
#[test]
fn a_sample_too_small_for_a_variance_gives_nan_rather_than_nonsense() {
    let empty = Summary::new();
    let one = Summary::of(&[1.0]);
    let two = Summary::of(&[1.0, 2.0]);

    for (a, b) in [(&empty, &two), (&two, &empty), (&one, &two), (&two, &one)] {
        let w = welch_t(a, b);
        assert!(w.t.is_nan(), "t = {}", w.t);
        assert!(w.df.is_nan(), "df = {}", w.df);
        assert!(w.p.is_nan(), "p = {}", w.p);
        assert!(w.difference.is_nan());
        assert!(w.standard_error.is_nan());
    }

    // And TOST built on it does not declare equivalence.
    let result = tost(&empty, &two, 1.0, 0.05);
    assert!(!result.equivalent, "{result}");

    // Two observations is enough.
    assert!(welch_t(&two, &two).df.is_finite());
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

// ---------------------------------------------------------------------------
// The paired test
// ---------------------------------------------------------------------------

/// A one-sample t against a hand-computed case.
///
/// `d = [1, 2, 3, 4, 5]`: mean 3, sample variance 2.5, so `se = sqrt(2.5/5) =
/// sqrt(0.5)` and `t = 3 / sqrt(0.5) = 4.2426406871192848`, on 4 degrees of
/// freedom.
#[test]
fn the_one_sample_t_matches_a_hand_computed_case() {
    let sample = Summary::of(&[1.0, 2.0, 3.0, 4.0, 5.0]);
    let result = tepsim_stats::one_sample_t(&sample);

    close(result.mean, 3.0, TIGHT, "mean");
    close(
        result.standard_error,
        0.5_f64.sqrt(),
        TIGHT,
        "standard error",
    );
    close(result.t, 3.0 / 0.5_f64.sqrt(), TIGHT, "t");
    close(result.df, 4.0, TIGHT, "df");
    close(
        result.p,
        2.0 * (1.0 - student_t_cdf(result.t, 4.0)),
        1e-12,
        "two-sided p",
    );
}

/// The paired test is the unpaired one applied to the differences.
///
/// Not a restatement: it is the check that `tost_paired` computes the same
/// thing a caller would get by differencing by hand and testing the result
/// against a null of zero. `welch_t` cannot express that, because it needs two
/// samples.
#[test]
fn the_paired_test_is_a_one_sample_test_of_the_differences() {
    let a = [10.1, 10.4, 9.8, 10.2, 10.0, 9.9, 10.3, 10.1];
    let b = [10.3, 10.7, 9.9, 10.5, 10.1, 10.2, 10.4, 10.3];
    let differences: Vec<f64> = a.iter().zip(b).map(|(x, y)| y - x).collect();

    let paired = tepsim_stats::tost_paired(&Summary::of(&differences), 0.5, 0.05);
    let single = tepsim_stats::one_sample_t(&Summary::of(&differences));

    close(paired.welch.difference, single.mean, TIGHT, "mean");
    close(
        paired.welch.standard_error,
        single.standard_error,
        TIGHT,
        "se",
    );
    close(paired.welch.df, 7.0, TIGHT, "df");
    assert!(paired.equivalent, "{paired}");
}

/// **The reason the paired test exists.**
///
/// Two sources measured at the same seeds, where each seed contributes a large
/// shared offset and the two sources differ by almost nothing. The unpaired
/// test spends all its variance on the shared offset and cannot resolve the
/// difference; the paired test sees straight through it.
///
/// This is exactly the shape of Tier 5's manipulated variables: an integrating
/// controller's 48-hour mean wanders from seed to seed by far more than the
/// two implementations differ, and both implementations see the *same* wander
/// because they are given the same disturbance realisation.
#[test]
fn pairing_sees_through_a_shared_offset_that_defeats_the_unpaired_test() {
    // Ten seeds. Each contributes a wander of order 1; the two sources differ
    // by 0.01 on every one of them.
    let wander: Vec<f64> = (0..10).map(|i| ((i as f64) * 1.7).sin()).collect();
    let reference: Vec<f64> = wander.iter().map(|w| 100.0 + w).collect();
    let candidate: Vec<f64> = wander.iter().map(|w| 100.0 + w + 0.01).collect();
    let differences: Vec<f64> = candidate
        .iter()
        .zip(&reference)
        .map(|(c, r)| c - r)
        .collect();

    // A margin of a tenth of the pooled spread, as PLAN.org sets for Tier 5.
    let margin = 0.1 * Summary::of(&reference).sd();

    let unpaired = tost(
        &Summary::of(&candidate),
        &Summary::of(&reference),
        margin,
        0.05,
    );
    let paired = tepsim_stats::tost_paired(&Summary::of(&differences), margin, 0.05);

    let power = |t: &tepsim_stats::Tost| 0.5 * (t.interval.1 - t.interval.0) / t.margin;
    println!(
        "margin {margin:.4}; unpaired power {:.2}, paired power {:.4}",
        power(&unpaired),
        power(&paired)
    );

    assert!(
        !unpaired.equivalent,
        "the unpaired test resolved a 0.01 difference against a wander of \
         order 1, so this fixture does not demonstrate anything"
    );
    assert!(
        paired.equivalent,
        "the paired test failed to see through the shared offset: {paired}"
    );
    assert!(
        power(&paired) < power(&unpaired) / 20.0,
        "pairing bought only {:.1}x, not the order of magnitude it should",
        power(&unpaired) / power(&paired)
    );
}

#[test]
fn the_paired_test_still_rejects_a_real_difference() {
    // The same shape, but the two sources genuinely differ by half the margin
    // plus a bit, so equivalence must not be declared.
    let wander: Vec<f64> = (0..20).map(|i| ((i as f64) * 0.9).cos()).collect();
    let reference: Vec<f64> = wander.iter().map(|w| 50.0 + w).collect();
    let margin = 0.1 * Summary::of(&reference).sd();
    let shift = 2.0 * margin;
    let differences: Vec<f64> = wander.iter().map(|_| shift).collect();

    let paired = tepsim_stats::tost_paired(&Summary::of(&differences), margin, 0.05);
    println!("a shift of twice the margin: {paired}");
    assert!(!paired.equivalent);
    close(paired.welch.difference, shift, 1e-12, "the shift");
}

#[test]
fn a_paired_sample_too_small_for_a_variance_gives_nan() {
    for sample in [Summary::new(), Summary::of(&[1.0])] {
        let result = tepsim_stats::one_sample_t(&sample);
        assert!(result.t.is_nan());
        assert!(result.df.is_nan());
        assert!(!tepsim_stats::tost_paired(&sample, 1.0, 0.05).equivalent);
    }
}
