//! Known-answer tests for PCA, dynamic PCA, the two control limits, and the
//! detection metrics.
//!
//! Same rule as `known_answers.rs`: **nothing is checked against a number this
//! project produced.** The references here are
//!
//! - data matrices whose sample correlation matrix is an exact rational, so the
//!   spectrum and the loadings are known on paper;
//! - exact identities that hold for any data at all, such as the sum of
//!   T-squared over the training set being `(n - 1) k`;
//! - published chi-squared and F table values, reached through the limits'
//!   large-sample behaviour;
//! - the Wilson-Hilferty approximation, which the Jackson-Mudholkar SPE limit
//!   provably collapses to when the residual eigenvalues are equal;
//! - Monte Carlo calibration, where the exact answer is 5% by construction and
//!   the tolerance comes from the binomial standard error rather than from
//!   whatever made the test pass.
//!
//! # What this suite is known to catch
//!
//! Passing tests say nothing about what a suite would notice. So the suite was
//! measured: thirty-two deliberate defects were introduced one at a time across
//! `eigen.rs`, `pca.rs`, `dpca.rs`, `detection.rs` and `distribution.rs`, and
//! the whole crate was run against each. Thirty were caught by between one and
//! seven tests. Two were not, and both gaps were real:
//!
//! - Evaluating the annihilated Jacobi pivot instead of assigning zero to it.
//!   Nothing looked at `off(A)` closely enough to see a few ulp there.
//!   `a_rotation_annihilates_its_pivot_exactly` in `eigen_known_answers.rs`
//!   now does.
//! - Computing a correlation as `dot / (sqrt(ss_i) sqrt(ss_j))` rather than
//!   `dot / sqrt(ss_i ss_j)`. The difference is one ulp on exactly collinear
//!   columns and zero elsewhere, so every comparison with a tolerance missed
//!   it. `exactly_collinear_columns_give_a_correlation_of_exactly_one` is the
//!   exact-equality test that does not.
//!
//! The two thinnest survivors after that, each caught by exactly one test, are
//! `F(k, n)` in place of `F(k, n - k)` in the T-squared limit
//! (`the_t_squared_limit_matches_a_hand_computation_at_small_n`) and reversing
//! the lag blocks inside a scored window
//! (`a_scored_window_is_the_augmented_row_it_should_be`). Both are places where
//! deleting one test would leave a defect undetected.

#![allow(
    clippy::float_cmp,
    reason = "known-answer tests assert exact values where the answer is exact"
)]
#![allow(
    clippy::suboptimal_flops,
    reason = "expressions mirror the closed forms they check"
)]
#![allow(
    clippy::needless_range_loop,
    reason = "matrix code indexes row-major storage by (i, j); enumerate() obscures it"
)]

use tepsim_stats::detection::DetectionReport;
use tepsim_stats::pca::{ControlLimits, Statistics};
use tepsim_stats::{
    CorrelationMatrix, Dpca, Pca, Retention, alarms_above, augment_with_lags, detection_delay,
    detection_report, false_alarm_rate, fault_detection_rate, spe_limit, t_squared_limit,
};

/// A deterministic generator, so these tests are reproducible without pulling
/// in a random-number crate. Not a good generator; good enough to make a
/// scatter that is not a pattern.
struct Lcg(u64);

impl Lcg {
    fn new(seed: u64) -> Self {
        Self(seed | 1)
    }
    fn step(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1);
        self.0
    }
    /// Uniform on (0, 1). Strictly open at zero, so `ln` below never sees it.
    fn unit(&mut self) -> f64 {
        ((self.step() >> 11) as f64 + 0.5) / ((1_u64 << 53) as f64)
    }
}

/// Standard normal deviates by the Box-Muller transform.
///
/// Box, G. E. P. and Muller, M. E. (1958), "A note on the generation of random
/// normal deviates", *Annals of Mathematical Statistics* 29(2), 610-611.
///
/// Box-Muller rather than the ziggurat because it is four lines and exactly
/// correct given a uniform source: the ziggurat is faster and has a rejection
/// loop whose correctness is a separate thing to get right, and nothing here is
/// generator bound.
struct Normals {
    lcg: Lcg,
    spare: Option<f64>,
}

impl Normals {
    fn new(seed: u64) -> Self {
        Self {
            lcg: Lcg::new(seed),
            spare: None,
        }
    }
    fn next(&mut self) -> f64 {
        if let Some(z) = self.spare.take() {
            return z;
        }
        let radius = (-2.0 * self.lcg.unit().ln()).sqrt();
        let angle = 2.0 * std::f64::consts::PI * self.lcg.unit();
        self.spare = Some(radius * angle.sin());
        radius * angle.cos()
    }
    fn matrix(&mut self, samples: usize, variables: usize) -> Vec<f64> {
        (0..samples * variables).map(|_| self.next()).collect()
    }
}

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
// Fixtures whose correlation matrix is an exact rational
// ---------------------------------------------------------------------------
//
// Built from three of the four Hadamard columns of order 4, which are mutually
// orthogonal, sum to zero, and have a sum of squares of exactly 4:
//
//     a = ( 1, -1,  1, -1)
//     b = ( 1,  1, -1, -1)
//
// Taking `x3 = 3a + 4b` gives a column with a sum of squares of exactly
// 9*4 + 16*4 = 100, and inner products with `a` and `b` of exactly 12 and 16.
// Since the correlation is `dot / sqrt(ss_i ss_j)` and `sqrt(4 * 100) = 20`
// exactly, the sample correlations are exactly 12/20 and 16/20. Every number
// in the construction is an integer that `f64` holds exactly, so the matrix
// that reaches the eigensolver carries no rounding at all.

/// `x1 = a`, `x2 = 3a + 4b`. Correlation exactly 0.6.
const TWO_VARIABLE: [f64; 8] = [
    1.0, 7.0, //
    -1.0, 1.0, //
    1.0, -1.0, //
    -1.0, -7.0,
];

/// `x1 = a`, `x2 = b`, `x3 = 3a + 4b`. Correlations 0, 0.6 and 0.8, and rank 2:
/// the third standardised column is exactly `0.6 z1 + 0.8 z2`.
const THREE_VARIABLE: [f64; 12] = [
    1.0, 1.0, 7.0, //
    -1.0, 1.0, 1.0, //
    1.0, -1.0, -1.0, //
    -1.0, -1.0, -7.0,
];

/// The 2-by-2 case, worked out on paper.
///
/// `R = [[1, 0.6], [0.6, 1]]` has eigenvalues `1 +- 0.6`, so 1.6 and 0.4, with
/// eigenvectors `(1, 1)/sqrt(2)` and `(1, -1)/sqrt(2)`. Every equicorrelated
/// 2-by-2 has those eigenvectors, whatever the correlation is, which is what
/// makes them checkable without a reference implementation.
#[test]
fn pca_recovers_a_hand_computed_two_variable_spectrum() {
    let pca = Pca::fit(&TWO_VARIABLE, 4, 2, Retention::CumulativeVariance(0.9));

    close(pca.eigenvalues()[0], 1.6, 1e-15, "leading eigenvalue");
    close(pca.eigenvalues()[1], 0.4, 1e-15, "trailing eigenvalue");
    close(pca.total_variance(), 2.0, 1e-15, "trace");
    assert_eq!(pca.rank(), 2);
    assert!(pca.constant_columns().is_empty());

    // 1.6 / 2 = 0.8, short of 0.9, so the rule has to take the second one too.
    assert_eq!(pca.retained(), 2, "cumulative variance at 90%");
    close(pca.explained_variance(), 1.0, 1e-15, "explained variance");

    // Both loading vectors have two components of equal magnitude, so the sign
    // convention has a tie to break and an ulp decides it. Compared up to sign
    // for that reason; `eigen_known_answers.rs` says more about it.
    let root_half = 0.5_f64.sqrt();
    for (j, expected) in [[root_half, root_half], [root_half, -root_half]]
        .iter()
        .enumerate()
    {
        let got = [pca.loadings()[j], pca.loadings()[2 + j]];
        let plus = (got[0] - expected[0])
            .abs()
            .max((got[1] - expected[1]).abs());
        let minus = (got[0] + expected[0])
            .abs()
            .max((got[1] + expected[1]).abs());
        assert!(
            plus.min(minus) < 1e-15,
            "loading {j} is {got:?}, expected +-{expected:?}"
        );
    }

    // With both components retained the model spans the whole space, so the
    // squared prediction error of anything at all vanishes. To rounding, not
    // exactly: the residual is `z - P P' z`, which leaves the `O(eps ||z||)` by
    // which `P` falls short of exactly orthonormal. Squared, that is about
    // 1e-31, and 1e-28 is four orders above it.
    for sample in [[0.0, 0.0], [1.0, 7.0], [-3.0, 11.0]] {
        let spe = pca.spe(&sample);
        assert!(spe < 1e-28, "SPE at {sample:?} is {spe:.3e}");
    }

    // And T-squared is then the full Mahalanobis distance, which for this R
    // has a closed form:
    //
    //   R^-1 = 1/0.64 [[1, -0.6], [-0.6, 1]]
    //   T^2  = (z1^2 - 1.2 z1 z2 + z2^2) / 0.64
    //
    // The standard deviations are sqrt(4/3) and sqrt(100/3), so the
    // standardised values are known too.
    let sd = [(4.0_f64 / 3.0).sqrt(), (100.0_f64 / 3.0).sqrt()];
    for sample in [[1.0, 7.0], [-1.0, 1.0], [0.0, 0.0], [2.5, -4.0]] {
        let z = [sample[0] / sd[0], sample[1] / sd[1]];
        let expected = (z[0] * z[0] - 1.2 * z[0] * z[1] + z[1] * z[1]) / 0.64;
        close(
            pca.t_squared(&sample),
            expected,
            1e-13,
            &format!("T^2 at {sample:?}"),
        );
    }
}

/// The 3-by-3 case, also worked out on paper, including an exact zero
/// eigenvalue.
///
/// `R = [[1, 0, 0.6], [0, 1, 0.8], [0.6, 0.8, 1]]`. Writing `R = I + M`, the
/// eigenvalues of `M` are `+-1` and `0`, so `R` has eigenvalues **2, 1 and 0**
/// with eigenvectors
///
/// ```text
/// lambda = 2:  (0.6, 0.8, 1) / sqrt(2)
/// lambda = 1:  (0.8, -0.6, 0)
/// lambda = 0:  (-0.6, -0.8, 1) / sqrt(2)
/// ```
///
/// The zero is not an accident of the numbers: the third standardised column is
/// exactly `0.6 z1 + 0.8 z2`, so the data really is rank two. That makes this
/// the test for the rank guard, which must refuse to retain the third component
/// however greedy the retention rule is, because T-squared divides by the
/// eigenvalue.
///
/// Each eigenvector here has a clear largest component, so unlike the 2-by-2
/// case the sign convention is unambiguous and the loadings are checked
/// component by component.
#[test]
fn pca_recovers_a_hand_computed_three_variable_spectrum_including_a_zero() {
    let pca = Pca::fit(&THREE_VARIABLE, 4, 3, Retention::CumulativeVariance(0.999));

    println!("3-variable spectrum: {:?}", pca.eigenvalues());
    close(pca.eigenvalues()[0], 2.0, 1e-15, "lambda 1");
    close(pca.eigenvalues()[1], 1.0, 1e-15, "lambda 2");
    assert!(
        pca.eigenvalues()[2].abs() < 1e-15,
        "lambda 3 should be zero, got {:.3e}",
        pca.eigenvalues()[2]
    );
    close(pca.total_variance(), 3.0, 1e-15, "trace");

    assert_eq!(pca.rank(), 2, "the data is exactly rank two");
    assert_eq!(
        pca.retained(),
        2,
        "the rank guard must stop the 99.9% rule at two"
    );
    // Even asking for all three explicitly.
    let greedy = Pca::fit(&THREE_VARIABLE, 4, 3, Retention::Fixed(3));
    assert_eq!(greedy.retained(), 2, "Fixed(3) past the rank");

    let root_half = 0.5_f64.sqrt();
    let expected = [
        [0.6 * root_half, 0.8 * root_half, root_half],
        [0.8, -0.6, 0.0],
        [-0.6 * root_half, -0.8 * root_half, root_half],
    ];
    for j in 0..3 {
        for i in 0..3 {
            close(
                pca.loadings()[i * 3 + j],
                expected[j][i],
                1e-14,
                &format!("loading[{i}][{j}]"),
            );
        }
    }

    // The residual subspace is the exact null direction, so every training row
    // reconstructs perfectly.
    for t in 0..4 {
        let row = &THREE_VARIABLE[t * 3..t * 3 + 3];
        assert!(
            pca.spe(row) < 1e-28,
            "training row {t} has SPE {:.3e} in a model that spans its span",
            pca.spe(row)
        );
    }
}

/// The correlation matrix PCA builds is the one the crate already tests,
/// bit for bit.
///
/// `pca.rs` forms it inline rather than calling [`CorrelationMatrix`], because
/// it needs a different answer for a constant column: zero variance there, not
/// `NaN`. Everywhere else the two must agree.
///
/// The comparison here is to 1e-14 rather than to the bit, because what is
/// available to compare against is `P L P'` rather than the matrix itself, and
/// that reconstruction rounds. The *grouping* inside the formula, which a
/// tolerance cannot see, is pinned by
/// `exactly_collinear_columns_give_a_correlation_of_exactly_one`.
#[test]
fn the_correlation_matrix_agrees_with_the_crate_s_own() {
    let mut rng = Normals::new(0x_C0FF_EE01);
    for (samples, variables) in [(4_usize, 3_usize), (37, 5), (200, 8)] {
        let data = rng.matrix(samples, variables);
        let pca = Pca::fit(&data, samples, variables, Retention::Fixed(1));

        let series: Vec<Vec<f64>> = (0..variables)
            .map(|v| (0..samples).map(|t| data[t * variables + v]).collect())
            .collect();
        let reference = CorrelationMatrix::of(&series);

        // The fitted model does not expose the matrix, so it is reconstructed
        // from the decomposition: `R = P L P'`. That is a stronger check than
        // comparing a stored copy would be, because it also pins the loadings
        // and eigenvalues to the matrix they came from.
        for i in 0..variables {
            for j in 0..variables {
                let mut rebuilt = 0.0;
                for k in 0..variables {
                    rebuilt += pca.loadings()[i * variables + k]
                        * pca.eigenvalues()[k]
                        * pca.loadings()[j * variables + k];
                }
                let want = reference.get(i, j);
                assert!(
                    (rebuilt - want).abs() < 1e-14,
                    "({samples}x{variables}) R[{i}][{j}]: model says {rebuilt:.17e}, \
                     CorrelationMatrix says {want:.17e}"
                );
            }
        }
        // The diagonal is exactly one, which is what makes the trace exactly
        // the variable count and the variance fractions exact.
        close(
            pca.total_variance(),
            variables as f64,
            1e-14,
            &format!("trace of a {variables}-variable correlation matrix"),
        );
    }
}

/// Two exactly proportional columns correlate exactly 1, and the spectrum is
/// exactly `[2, 0]`.
///
/// This is the case `correlation.rs` writes its grouping for. With
/// `x = [1, 2, 3, 4]` and `y = 2x` the centred sums of squares are 5 and 20,
/// whose product is exactly 100 and whose square root is exactly 10, so
/// `dot / sqrt(ss_i ss_j)` returns exactly 1. The other grouping,
/// `dot / (sqrt(ss_i) sqrt(ss_j))`, gives `10 / 10.000000000000002`, which is
/// one ulp below 1, and the eigenvalues come out `1.9999999999999998` and
/// `2.2e-16` instead of 2 and 0.
///
/// So the assertions are `assert_eq!`, deliberately. A tolerance of 1e-15 would
/// pass on both and this test would be worth nothing.
///
/// # Provenance
///
/// A mutation run that swapped the grouping in `pca.rs` produced no failures
/// across the crate's 136 tests. The correlation matrix comparison above uses a
/// tolerance and could not see it, and its name claimed otherwise.
#[test]
fn exactly_collinear_columns_give_a_correlation_of_exactly_one() {
    // x = 1, 2, 3, 4 and y = 2x, as columns of a 4-by-2 matrix.
    let data = [1.0, 2.0, 2.0, 4.0, 3.0, 6.0, 4.0, 8.0];
    let pca = Pca::fit(&data, 4, 2, Retention::Fixed(2));

    assert_eq!(pca.eigenvalues()[0], 2.0, "leading eigenvalue");
    assert_eq!(pca.eigenvalues()[1], 0.0, "trailing eigenvalue");
    assert_eq!(pca.total_variance(), 2.0);
    assert_eq!(pca.rank(), 1, "two proportional columns are rank one");
    assert_eq!(pca.retained(), 1, "Fixed(2) clamped by the rank");

    // The demonstration that the other grouping is a different number, so this
    // test is not asserting something that holds either way.
    let (ss_x, ss_y) = (5.0_f64, 20.0_f64);
    assert_eq!(10.0 / (ss_x * ss_y).sqrt(), 1.0);
    assert!(
        10.0 / (ss_x.sqrt() * ss_y.sqrt()) < 1.0,
        "the two groupings agree here, so this fixture cannot tell them apart"
    );
}

// ---------------------------------------------------------------------------
// Exact identities that hold for any data
// ---------------------------------------------------------------------------

/// Two identities that pin standardisation, the `n - 1` divisor and the
/// projection all at once.
///
/// Let `Z` be the standardised training matrix, so `Z'Z = (n - 1) R`. Then for
/// any `k`,
///
/// ```text
/// sum_t T^2_t   = tr(L_k^-1 P_k' Z'Z P_k) = (n - 1) k
/// sum_t SPE_t   = tr(Z'Z) - tr(P_k' Z'Z P_k) = (n - 1) sum_{j >= k} lambda_j
/// ```
///
/// Both are exact for **any** data, and both are sensitive to exactly the
/// mistakes that are otherwise invisible:
///
/// - Standardising with `n` instead of `n - 1` scales every `z` by
///   `sqrt(n / (n - 1))` and both sums by `n / (n - 1)`. At `n = 100` that is
///   1%, against a tolerance of 1e-12.
/// - Skipping standardisation entirely changes both sums by whatever the
///   variable scales happen to be.
/// - Projecting onto the wrong subspace, or summing the discarded scores
///   instead of measuring the residual, breaks the second one.
///
/// Note that the *correlation matrix itself* is invariant to the `n` versus
/// `n - 1` choice, since `dot / sqrt(ss_i ss_j)` has no divisor in it at all.
/// So the eigenvalues cannot detect that mistake and these sums are the only
/// thing here that can.
#[test]
fn the_training_set_sums_of_both_statistics_are_exact() {
    let mut rng = Normals::new(0x_1DEA_0001);
    for (samples, variables) in [(20_usize, 3_usize), (100, 6), (500, 12)] {
        let data = rng.matrix(samples, variables);
        for k in 1..=variables {
            let pca = Pca::fit(&data, samples, variables, Retention::Fixed(k));
            assert_eq!(pca.retained(), k);

            let mut total_t2 = 0.0;
            let mut total_spe = 0.0;
            for t in 0..samples {
                let row = &data[t * variables..(t + 1) * variables];
                let Statistics { t_squared, spe } = pca.statistics(row);
                total_t2 += t_squared;
                total_spe += spe;
            }

            let expected_t2 = (samples - 1) as f64 * k as f64;
            let discarded: f64 = pca.eigenvalues()[k..].iter().sum();
            let expected_spe = (samples - 1) as f64 * discarded;

            close(
                total_t2,
                expected_t2,
                1e-12,
                &format!("sum of T^2 over {samples} training rows at k={k}"),
            );
            if expected_spe > 0.0 {
                close(
                    total_spe,
                    expected_spe,
                    1e-12,
                    &format!("sum of SPE over {samples} training rows at k={k}"),
                );
            } else {
                assert!(total_spe < 1e-20, "SPE with k = p should vanish");
            }
        }
    }
    println!("training-set sum identities hold at 1e-12 over 3 shapes and every k");
}

/// The scores are an orthogonal change of basis, so they conserve length, and
/// SPE is exactly the discarded part of it.
///
/// ```text
/// sum_j score_j^2       = sum_i z_i^2
/// sum_{j >= k} score_j^2 = SPE
/// ```
///
/// The second is the identity the implementation deliberately does *not* use:
/// it forms the residual by subtracting the reconstruction, because that
/// measures what the model actually failed to explain rather than what an
/// assumed completeness relation says it should have. Checking the two agree is
/// therefore a real check on the loadings' orthonormality, not a tautology.
#[test]
fn the_scores_conserve_length_and_account_for_the_residual() {
    let mut rng = Normals::new(0x5C0E);
    let (samples, variables) = (80_usize, 7_usize);
    let data = rng.matrix(samples, variables);
    let mut worst_parseval = 0.0_f64;
    let mut worst_residual = 0.0_f64;

    for k in [1_usize, 3, 5, 7] {
        let pca = Pca::fit(&data, samples, variables, Retention::Fixed(k));
        for t in 0..samples {
            let row = &data[t * variables..(t + 1) * variables];
            let z = pca.standardise(row);
            let scores = pca.scores(row);
            let length: f64 = z.iter().map(|v| v * v).sum();
            let score_length: f64 = scores.iter().map(|s| s * s).sum();
            worst_parseval = worst_parseval.max((score_length - length).abs() / length);

            let discarded: f64 = scores[k..].iter().map(|s| s * s).sum();
            let spe = pca.spe(row);
            if discarded > 0.0 {
                worst_residual = worst_residual.max((spe - discarded).abs() / discarded);
            } else {
                assert!(spe < 1e-25);
            }
        }
    }
    println!(
        "scores: worst relative Parseval defect {worst_parseval:.3e}, worst \
         relative disagreement between SPE and the discarded scores {worst_residual:.3e}"
    );
    assert!(worst_parseval < 1e-13, "{worst_parseval:.3e}");
    assert!(worst_residual < 1e-12, "{worst_residual:.3e}");
}

/// Standardising is what stops the biggest unit from winning, and here is the
/// case that shows it.
///
/// The same two variables as `TWO_VARIABLE`, with the second multiplied by a
/// thousand. Nothing about the relationship changed, so the correlation matrix,
/// the eigenvalues and the loadings must all be **identical**, and T-squared
/// and SPE with them.
///
/// Covariance PCA would give the opposite answer, and the test computes what it
/// would give rather than asserting that it would be bad. With variances 4/3
/// and 1e6 * 100/3, the leading covariance eigenvector is
/// `(1, 0)` to eleven decimal places in the *scaled* coordinates, that is, the
/// second variable alone. A detector built on that is a scale meter.
#[test]
fn standardising_makes_the_model_immune_to_the_units() {
    let plain = Pca::fit(&TWO_VARIABLE, 4, 2, Retention::Fixed(1));
    let scaled_data: Vec<f64> = TWO_VARIABLE
        .iter()
        .enumerate()
        .map(|(index, x)| if index % 2 == 1 { x * 1000.0 } else { *x })
        .collect();
    let scaled = Pca::fit(&scaled_data, 4, 2, Retention::Fixed(1));

    for j in 0..2 {
        assert_eq!(
            plain.eigenvalues()[j].to_bits(),
            scaled.eigenvalues()[j].to_bits(),
            "eigenvalue {j} moved when the units changed"
        );
    }
    for (a, b) in plain.loadings().iter().zip(scaled.loadings()) {
        assert_eq!(a.to_bits(), b.to_bits(), "a loading moved with the units");
    }
    for t in 0..4 {
        let plain_row = &TWO_VARIABLE[t * 2..t * 2 + 2];
        let scaled_row = &scaled_data[t * 2..t * 2 + 2];
        close(
            scaled.t_squared(scaled_row),
            plain.t_squared(plain_row),
            1e-13,
            &format!("T^2 of row {t} under a change of units"),
        );
    }

    // What a covariance model would have done, from the 2-by-2 closed form.
    // Sample variances 4/3 and 1e6 * 100/3, covariance 12/3 * 1000.
    let (va, vd) = (4.0_f64 / 3.0, 1e6 * 100.0 / 3.0);
    let cov = 1000.0_f64 * 12.0 / 3.0;
    let theta = 0.5 * (2.0 * cov).atan2(va - vd);
    let leading = (theta.cos(), theta.sin());
    println!(
        "covariance PCA on the scaled data would put its first component at \
         ({:.9}, {:.9}); the correlation model puts it at ({:.9}, {:.9})",
        leading.0,
        leading.1,
        scaled.loadings()[0],
        scaled.loadings()[2]
    );
    assert!(
        leading.0.abs() < 0.01,
        "the covariance model is not dominated by the scaled variable here, so \
         this test does not demonstrate anything: leading component {leading:?}"
    );
}

/// A column that never moved is reported and contributes nothing.
#[test]
fn a_constant_column_is_named_and_ignored() {
    // Three variables, the middle one held fixed. TEP does exactly this with
    // XMV(12), the agitator speed.
    let data = [
        1.0, 5.0, 7.0, //
        -1.0, 5.0, 1.0, //
        1.0, 5.0, -1.0, //
        -1.0, 5.0, -7.0,
    ];
    let pca = Pca::fit(&data, 4, 3, Retention::Fixed(3));
    assert_eq!(pca.constant_columns(), &[1]);
    assert_eq!(pca.standard_deviations()[1], 0.0);
    close(pca.means()[1], 5.0, 1e-15, "mean of the constant column");

    // The trace is the number of variables that moved, and the constant one
    // contributes an eigenvalue of exactly zero.
    close(pca.total_variance(), 2.0, 1e-15, "trace");
    assert_eq!(pca.rank(), 2, "one direction carries no variance at all");
    assert_eq!(pca.retained(), 2, "Fixed(3) clamped by the rank");

    // Standardised to exactly zero whatever the observation says, so a value
    // that never appeared in training still adds nothing to either statistic.
    // Not because it is uninteresting, but because there is no scale to measure
    // it against and a raw deviation in the original units cannot be added to a
    // sum of standardised squares.
    let moved = [1.0, 900.0, 7.0];
    assert_eq!(pca.standardise(&moved)[1], 0.0);
    let held = [1.0, 5.0, 7.0];
    assert_eq!(pca.statistics(&moved), pca.statistics(&held));
}

/// Each retention rule does what its documentation says, on a spectrum that
/// distinguishes them.
#[test]
fn the_retention_rules_pick_what_they_claim() {
    // Six variables in two correlated blocks, so the spectrum is spread out.
    let mut rng = Normals::new(0x4E7A);
    let samples = 400;
    let mut data = vec![0.0; samples * 6];
    for t in 0..samples {
        let (f1, f2) = (rng.next(), rng.next());
        for v in 0..3 {
            data[t * 6 + v] = f1 + 0.35 * rng.next();
        }
        for v in 3..6 {
            data[t * 6 + v] = f2 + 0.9 * rng.next();
        }
    }

    let spectrum = Pca::fit(&data, samples, 6, Retention::Fixed(6));
    let eigenvalues = spectrum.eigenvalues().to_vec();
    println!("retention fixture spectrum: {eigenvalues:?}");

    // Kaiser: exactly the count above one.
    let above_one = eigenvalues.iter().filter(|&&v| v > 1.0).count();
    assert_eq!(
        Pca::fit(&data, samples, 6, Retention::Kaiser).retained(),
        above_one
    );
    assert!(
        (1..6).contains(&above_one),
        "the fixture's spectrum makes Kaiser trivial: {above_one} of 6 above one"
    );

    // Cumulative variance: the smallest k reaching the fraction, computed here
    // from the spectrum directly.
    for fraction in [0.5_f64, 0.7, 0.9, 0.95, 0.999] {
        let mut running = 0.0;
        let mut expected = 0;
        while expected < 6 && running / 6.0 < fraction {
            running += eigenvalues[expected];
            expected += 1;
        }
        let pca = Pca::fit(&data, samples, 6, Retention::CumulativeVariance(fraction));
        assert_eq!(
            pca.retained(),
            expected,
            "cumulative variance at {fraction} kept {} of 6",
            pca.retained()
        );
        assert!(
            pca.explained_variance() >= fraction,
            "kept {} components explaining {:.4}, short of {fraction}",
            pca.retained(),
            pca.explained_variance()
        );
        // And the rule is minimal: one fewer would not have reached it.
        if expected > 1 {
            let short: f64 = eigenvalues[..expected - 1].iter().sum::<f64>() / 6.0;
            assert!(short < fraction, "{expected} components was not minimal");
        }
    }

    // Fixed: exactly what was asked for, up to the rank.
    for k in 1..=6 {
        assert_eq!(
            Pca::fit(&data, samples, 6, Retention::Fixed(k)).retained(),
            k
        );
    }
}

// ---------------------------------------------------------------------------
// The T-squared control limit
// ---------------------------------------------------------------------------

/// Published upper 5% points of the chi-squared distribution.
const CHI_SQUARED_95: &[(usize, f64)] = &[
    (1, 3.841_458_820_694),
    (2, 5.991_464_547_108),
    (3, 7.814_727_903_252),
    (5, 11.070_497_693_516),
    (10, 18.307_038_053_275),
    (15, 24.995_790_139_064),
    (20, 31.410_432_844_231),
    (30, 43.772_971_808_505),
];

/// As the training set grows the T-squared limit becomes the chi-squared
/// quantile, and here it does, against published tables.
///
/// ```text
/// k (n + 1)(n - 1) / (n (n - k)) F_alpha(k, n - k)  ->  chi^2_{k, alpha}
/// ```
///
/// because the leading factor tends to 1 and `k F(k, m) -> chi^2_k` as
/// `m -> infinity`. This is the check that pins the whole expression, degrees of
/// freedom included: using `F(k, n)` instead of `F(k, n - k)` survives it (both
/// tend to the same place), but using `F(k - 1, n - k)` or dropping the leading
/// factor does not.
///
/// # Why `n = 1e8` and not more
///
/// Two errors pull in opposite directions. The finite-`n` correction is
/// `O(k / n)`, so it wants `n` large. The incomplete beta's accuracy at
/// `df2 = n - k` degrades like `eps * lnGamma(n / 2)`, which at `n = 1e8` is
/// `eps * 8.4e8 = 1.9e-7`, so it wants `n` small. They cross near 1e8, where the
/// total is a few times 1e-7. `known_answers.rs` documents the same trade for
/// the t distribution approaching the normal, at the same root cause.
#[test]
fn the_t_squared_limit_becomes_the_chi_squared_quantile() {
    for &(k, chi_squared) in CHI_SQUARED_95 {
        let limit = t_squared_limit(k, 100_000_000, 0.95);
        let error = (limit - chi_squared).abs() / chi_squared;
        println!(
            "T^2 limit at k={k}, n=1e8: {limit:.10}, chi^2_{{{k},0.95}} = \
             {chi_squared:.10}, relative difference {error:.3e}"
        );
        assert!(error < 1e-6, "k={k}: {error:.3e}");
    }
}

/// The limit's algebra, against a published F value and a hand multiplication.
///
/// `k = 2`, `n = 12`:
///
/// ```text
/// 2 * 13 * 11 / (12 * 10) = 286 / 120 = 2.3833333...
/// F_0.95(2, 10) = 4.1028210151304  (exact closed form: 5 (0.05^-0.2 - 1))
/// limit = 9.7783900861...
/// ```
///
/// The F value is the 4.10 of the published table, and `known_answers.rs`
/// pins it against both the table and the closed form. Small `n` is the point:
/// this is where every part of the coefficient matters, and where the wrong
/// degrees of freedom or the wrong `n +- 1` shows up as percent, not as 1e-8.
#[test]
fn the_t_squared_limit_matches_a_hand_computation_at_small_n() {
    // The closed form for F_p(2, d2), independent of any table.
    let f = 10.0 / 2.0 * (0.05_f64.powf(-2.0 / 10.0) - 1.0);
    close(f, 4.102_821_015_130_4, 1e-13, "F_0.95(2, 10) closed form");
    let expected = 2.0 * 13.0 * 11.0 / (12.0 * 10.0) * f;
    println!("T^2 limit at k=2, n=12, 95%: hand value {expected:.10}");
    close(
        t_squared_limit(2, 12, 0.95),
        expected,
        1e-14,
        "T^2 limit at k=2, n=12",
    );

    // And a second one at another confidence, to pin the F argument as well as
    // the coefficient: F_0.99(5, 2) has the closed form for df2 = 2.
    let y = 0.99_f64.powf(2.0 / 5.0);
    let f = 2.0 * y / (5.0 * (1.0 - y));
    let expected = 5.0 * 8.0 * 6.0 / (7.0 * 2.0) * f;
    println!("T^2 limit at k=5, n=7, 99%: hand value {expected:.10}");
    close(
        t_squared_limit(5, 7, 0.99),
        expected,
        1e-13,
        "T^2 limit at k=5, n=7",
    );
}

#[test]
fn the_t_squared_limit_reports_its_degenerate_cases() {
    // No components means T-squared is identically zero, so a limit of zero is
    // consistent rather than degenerate.
    assert_eq!(t_squared_limit(0, 100, 0.95), 0.0);
    // As many components as samples is not a fitted model.
    assert!(t_squared_limit(5, 5, 0.95).is_nan());
    assert!(t_squared_limit(5, 4, 0.95).is_nan());
    for confidence in [0.0, 1.0, -0.1, 1.5] {
        assert!(
            t_squared_limit(3, 100, confidence).is_nan(),
            "confidence {confidence}"
        );
    }
    // Monotone in the confidence, which any control limit must be.
    let mut previous = 0.0;
    for confidence in [0.5, 0.9, 0.95, 0.99, 0.999] {
        let limit = t_squared_limit(4, 200, confidence);
        assert!(limit > previous, "limit fell at confidence {confidence}");
        previous = limit;
    }
}

/// The T-squared limit, calibrated by Monte Carlo at a sample size where the
/// exact formula and the large-sample one differ by 17%.
///
/// # The design
///
/// For independent standard normal data, `T^2` of an observation **not** in the
/// training set is exactly `k (n+1)(n-1) / (n (n-k)) F(k, n-k)` distributed.
/// So the exceedance rate of the 95% limit is exactly 5%, not approximately: no
/// asymptotics are involved, and the only error is sampling.
///
/// `n = 60` training samples and `p = k = 5` variables are chosen because they
/// make the exact limit 12.993 while the chi-squared limit is 11.070. A wrong
/// formula is a 17% wrong threshold, which shows up as an 8.9% exceedance rate
/// rather than 5%.
///
/// 500 replicates of 200 test observations each. Replicates rather than one
/// long run because the training set is the thing being tested: within a
/// replicate all 200 observations share one estimated mean and covariance, so
/// they are not independent, and the honest standard error comes from the
/// spread **between** replicates. That spread is measured here rather than
/// assumed, and the test refuses to pass if it is so large that the assertion
/// would have been vacuous.
#[test]
fn the_t_squared_limit_is_calibrated_at_a_sample_size_where_it_matters() {
    const REPLICATES: usize = 500;
    const TRAINING: usize = 60;
    const TESTING: usize = 200;
    const VARIABLES: usize = 5;
    const CONFIDENCE: f64 = 0.95;

    let limit = t_squared_limit(VARIABLES, TRAINING, CONFIDENCE);
    println!(
        "exact T^2 limit at k={VARIABLES}, n={TRAINING}: {limit:.6}; the \
         large-sample chi^2 limit would be 11.070498, a {:.1}% difference",
        (limit / 11.070_497_693_516 - 1.0) * 100.0
    );
    assert!(
        (limit / 11.070_497_693_516 - 1.0).abs() > 0.1,
        "the exact and asymptotic limits agree to within 10% here, so this test \
         cannot tell them apart"
    );

    let mut rng = Normals::new(0x_CA11B);
    let mut rates = Vec::with_capacity(REPLICATES);
    let mut exceedances = 0_usize;
    for _ in 0..REPLICATES {
        let training = rng.matrix(TRAINING, VARIABLES);
        let pca = Pca::fit(&training, TRAINING, VARIABLES, Retention::Fixed(VARIABLES));
        assert_eq!(pca.retained(), VARIABLES);
        let mut hits = 0;
        for _ in 0..TESTING {
            let sample: Vec<f64> = (0..VARIABLES).map(|_| rng.next()).collect();
            if pca.t_squared(&sample) > limit {
                hits += 1;
            }
        }
        exceedances += hits;
        rates.push(hits as f64 / TESTING as f64);
    }

    let total = (REPLICATES * TESTING) as f64;
    let pooled = exceedances as f64 / total;
    let mean = rates.iter().sum::<f64>() / REPLICATES as f64;
    let variance =
        rates.iter().map(|r| (r - mean) * (r - mean)).sum::<f64>() / (REPLICATES - 1) as f64;
    let standard_error = (variance / REPLICATES as f64).sqrt();
    // What it would be if the 100,000 observations were independent, which they
    // are not: printed so the inflation from sharing a training set is visible.
    let binomial = (0.05 * 0.95 / total).sqrt();

    println!(
        "T^2 calibration: {exceedances} exceedances in {total:.0} observations, \
         rate {pooled:.5} against a nominal 0.05000. Between-replicate standard \
         deviation {:.5} over {REPLICATES} replicates, so the standard error of \
         the mean is {standard_error:.5}; the binomial standard error ignoring \
         the shared training set would be {binomial:.5}. Departure: {:.2} standard errors.",
        variance.sqrt(),
        (pooled - 0.05) / standard_error
    );

    // The test has teeth only if the standard error is small compared with the
    // effect a wrong limit would produce. The chi-squared limit would give
    // about 8.9%, which is 0.039 away; require the standard error to be at
    // least ten times smaller than that.
    assert!(
        standard_error < 0.0039,
        "the standard error {standard_error:.5} is too large for this test to \
         distinguish 5% from the 8.9% a wrong limit would give"
    );
    assert!(
        (pooled - 0.05).abs() < 3.0 * standard_error,
        "exceedance rate {pooled:.5} is {:.2} standard errors from the nominal 0.05",
        (pooled - 0.05).abs() / standard_error
    );
}

// ---------------------------------------------------------------------------
// The SPE control limit
// ---------------------------------------------------------------------------

/// With equal residual eigenvalues, Jackson-Mudholkar **is** Wilson-Hilferty.
///
/// Substituting `lambda_j = lambda` for all `m` residual components gives
/// `theta_i = m lambda^i` and `h0 = 1 - 2/3 = 1/3` exactly, and the limit
/// collapses to
///
/// ```text
/// SPE_alpha = m lambda [ 1 - 2/(9m) + z_alpha sqrt(2/(9m)) ]^3
/// ```
///
/// which is exactly Wilson and Hilferty's cube-root approximation to the
/// chi-squared quantile, times `lambda`. That is the right answer: SPE is then
/// `lambda` times a chi-squared with `m` degrees of freedom.
///
/// So this test checks two things at once. First, that the general formula
/// reduces to the special one to machine precision, which pins `h0`, the
/// `theta` power sums and the exponent. Second, that the special one lands near
/// the published chi-squared table, at the accuracy Wilson-Hilferty is known to
/// have, with the discrepancy shrinking as `m` grows, which is the shape a
/// converging approximation has and a wrong constant does not.
///
/// # Two bounds, because Wilson-Hilferty's accuracy depends on `m`
///
/// The transform normalises a chi-squared by taking its cube root, and how much
/// work that has to do depends on how skew the chi-squared is. At `m = 1` the
/// density is unbounded at zero and the 95% point comes out 2.5% low; by
/// `m = 10` the error is 0.08% and by `m = 30` it is 0.014%. So there are two
/// assertions: 3% everywhere, which is the "a few percent at worst" the
/// approximation is documented to give, and 0.1% from `m = 10` up, which is the
/// regime any real residual subspace is in. TEP has 52 measured variables and a
/// 90% model retains a dozen or so, leaving forty residual dimensions.
///
/// The monotone decrease is the third assertion and the structural one: a limit
/// wrong by a constant factor would sail through both bounds at large `m` and
/// fail this.
#[test]
fn the_spe_limit_collapses_to_wilson_hilferty_and_tracks_the_chi_squared_table() {
    let z = 1.644_853_626_951_47_f64;
    let mut previous_gap: Option<(usize, f64)> = None;
    for &(m, chi_squared) in CHI_SQUARED_95 {
        for lambda in [1.0_f64, 0.37, 12.5] {
            let residual = vec![lambda; m];
            let limit = spe_limit(&residual, 0.95);

            let mf = m as f64;
            let wilson_hilferty =
                mf * lambda * (1.0 - 2.0 / (9.0 * mf) + z * (2.0 / (9.0 * mf)).sqrt()).powi(3);
            close(
                limit,
                wilson_hilferty,
                1e-14,
                &format!("Jackson-Mudholkar at m={m}, lambda={lambda}"),
            );

            let gap = (limit - chi_squared * lambda).abs() / (chi_squared * lambda);
            if lambda == 1.0 {
                println!(
                    "SPE limit, m={m} equal residual eigenvalues: {limit:.9}, \
                     lambda chi^2_{{{m},0.95}} = {:.9}, relative gap {gap:.3e}",
                    chi_squared * lambda
                );
                // Wilson-Hilferty's documented accuracy, in two regimes.
                assert!(gap < 0.03, "m={m}: {gap:.3e}");
                if m >= 10 {
                    assert!(gap < 0.001, "m={m}: {gap:.3e}");
                }
                // And it converges: the gap shrinks as the degrees of freedom
                // grow. A limit wrong by a constant factor would not.
                if let Some((previous_m, previous)) = previous_gap {
                    assert!(
                        gap < previous,
                        "the gap grew from {previous:.3e} at m={previous_m} to \
                         {gap:.3e} at m={m}, so this is not a converging approximation"
                    );
                }
                previous_gap = Some((m, gap));
            }
        }
    }
}

/// SPE scales with its eigenvalues, exactly.
///
/// `SPE = sum_j lambda_j g_j^2` for standard normal `g`, so multiplying every
/// residual eigenvalue by `c` multiplies the statistic and its limit by `c`.
/// The Jackson-Mudholkar formula must be homogeneous of degree one, and it is:
/// each `theta_i` picks up `c^i`, `h0` is invariant, and the bracket is
/// invariant, leaving the leading `theta1`.
#[test]
fn the_spe_limit_is_homogeneous_in_the_residual_eigenvalues() {
    let base = [3.1, 2.0, 1.4, 0.9, 0.5, 0.22, 0.05];
    for confidence in [0.9_f64, 0.95, 0.99] {
        let reference = spe_limit(&base, confidence);
        for scale in [1e-6_f64, 0.25, 7.0, 1e6] {
            let scaled: Vec<f64> = base.iter().map(|l| l * scale).collect();
            close(
                spe_limit(&scaled, confidence),
                reference * scale,
                1e-13,
                &format!("SPE limit scaled by {scale} at {confidence}"),
            );
        }
    }
}

#[test]
fn the_spe_limit_reports_its_degenerate_cases() {
    // No residual subspace: SPE is identically zero, so is its limit.
    assert_eq!(spe_limit(&[], 0.95), 0.0);
    assert_eq!(spe_limit(&[0.0, 0.0], 0.95), 0.0);
    // Rounding noise below zero is clamped rather than cubed.
    let clamped = spe_limit(&[-1e-18, -2e-19, 1.0, 1.0], 0.95);
    let clean = spe_limit(&[0.0, 0.0, 1.0, 1.0], 0.95);
    assert_eq!(clamped.to_bits(), clean.to_bits());
    // Out of range confidences.
    for confidence in [-0.1, 1.0, 1.5] {
        assert!(spe_limit(&[1.0, 1.0], confidence).is_nan(), "{confidence}");
    }
    // Monotone in the confidence.
    let mut previous = 0.0;
    for confidence in [0.55, 0.9, 0.95, 0.99, 0.999] {
        let limit = spe_limit(&[2.0, 1.0, 0.5], confidence);
        assert!(limit > previous, "limit fell at {confidence}");
        previous = limit;
    }
}

/// The SPE limit, calibrated against the distribution it is approximating.
///
/// # The design
///
/// SPE is `sum_j lambda_j g_j^2` with independent standard normal `g_j`, when
/// the model is correct. That distribution is sampled here **directly**, from
/// the residual eigenvalues, rather than by fitting a model and reading its
/// residual: fitting would add estimation noise in the eigenvalues and a bias
/// from splitting a noisy spectrum at a fixed `k`, and neither has anything to
/// do with whether the Jackson-Mudholkar formula is right. This isolates the
/// approximation.
///
/// The eigenvalues are deliberately unequal, so the test exercises the general
/// formula rather than the equal-eigenvalue case that
/// `the_spe_limit_collapses_to_wilson_hilferty_and_tracks_the_chi_squared_table`
/// already pins exactly.
///
/// # The tolerance
///
/// Two contributions, both named:
///
/// - Sampling. 400,000 draws at a nominal 5% gives a binomial standard error of
///   `sqrt(0.05 * 0.95 / 400000) = 3.4e-4`.
/// - The approximation itself. Jackson and Mudholkar match three moments, and
///   the equal-eigenvalue case above measures its error on the *quantile* at
///   between 0.1% and 0.9% depending on the degrees of freedom. Near the 95%
///   point the density times the limit is about 0.3, so a 0.5% error on the
///   quantile moves the rate by about 0.0015.
///
/// So the band is `0.05 +- 0.005`, which is three times the second contribution
/// and fifteen times the first. It is wide enough to accommodate an
/// approximation and far too narrow to accommodate a wrong one: dropping the
/// `theta2 h0 (h0 - 1) / theta1^2` term, or using `h0` in place of `1 / h0` as
/// the exponent, moves the limit by tens of percent.
#[test]
fn the_spe_limit_is_calibrated_against_the_distribution_it_approximates() {
    const DRAWS: usize = 400_000;
    let residual = [3.0_f64, 2.0, 1.0, 1.0, 0.5, 0.2];
    let limit = spe_limit(&residual, 0.95);
    println!("SPE limit for {residual:?} at 95%: {limit:.9}");

    let mut rng = Normals::new(0x_5BEE_0007);
    let mut exceedances = 0_usize;
    for _ in 0..DRAWS {
        let statistic: f64 = residual
            .iter()
            .map(|l| {
                let g = rng.next();
                l * g * g
            })
            .sum();
        if statistic > limit {
            exceedances += 1;
        }
    }
    let rate = exceedances as f64 / DRAWS as f64;
    let binomial = (0.05 * 0.95 / DRAWS as f64).sqrt();
    println!(
        "SPE calibration: {exceedances} exceedances in {DRAWS} draws, rate \
         {rate:.5} against a nominal 0.05000; the binomial standard error is \
         {binomial:.5}, so the departure is {:.1} sampling standard errors, \
         which is the Jackson-Mudholkar approximation showing.",
        (rate - 0.05).abs() / binomial
    );
    assert!(
        (rate - 0.05).abs() < 0.005,
        "exceedance rate {rate:.5} is outside the 0.045 to 0.055 band"
    );
}

/// A fitted model's limits are the free functions' limits, with the right
/// arguments.
#[test]
fn a_model_reports_the_limits_of_its_own_shape() {
    let mut rng = Normals::new(0x_11417);
    let data = rng.matrix(300, 6);
    let pca = Pca::fit(&data, 300, 6, Retention::CumulativeVariance(0.8));
    let ControlLimits {
        t_squared,
        spe,
        confidence,
        components,
        samples,
    } = pca.limits(0.99);

    assert_eq!(confidence, 0.99);
    assert_eq!(components, pca.retained());
    assert_eq!(samples, 300);
    assert_eq!(
        t_squared.to_bits(),
        t_squared_limit(pca.retained(), 300, 0.99).to_bits()
    );
    assert_eq!(
        spe.to_bits(),
        spe_limit(&pca.eigenvalues()[pca.retained()..], 0.99).to_bits()
    );
}

// ---------------------------------------------------------------------------
// Dynamic PCA
// ---------------------------------------------------------------------------

/// The augmented layout, entry by entry, against a matrix small enough to write
/// down.
///
/// Off by one in the lag index is the classic dynamic-PCA bug and it is
/// invisible in the fitted spectrum, because a shifted history is still a
/// history. So this checks the layout directly.
#[test]
fn the_lag_augmented_matrix_has_the_layout_it_claims() {
    // Five samples of two variables, numbered so every entry says where it came
    // from: variable v at time t is 10 t + v.
    let data: Vec<f64> = (0..5)
        .flat_map(|t| (0..2).map(move |v| (10 * t + v) as f64))
        .collect();

    let (augmented, rows, columns) = augment_with_lags(&data, 5, 2, 2);
    assert_eq!(rows, 3, "two lags eat the first two samples");
    assert_eq!(columns, 6, "two variables times three time slots");
    // Row r is [ x(r+2) | x(r+1) | x(r) ], present first.
    assert_eq!(
        augmented,
        vec![
            20.0, 21.0, 10.0, 11.0, 0.0, 1.0, //
            30.0, 31.0, 20.0, 21.0, 10.0, 11.0, //
            40.0, 41.0, 30.0, 31.0, 20.0, 21.0,
        ]
    );

    // Zero lags is the identity, entry for entry.
    let (identity, rows, columns) = augment_with_lags(&data, 5, 2, 0);
    assert_eq!((rows, columns), (5, 2));
    assert_eq!(identity, data);
}

/// With no lags, dynamic PCA is PCA, bit for bit.
///
/// The wrapper is thin enough that this should hold exactly rather than
/// approximately, and asserting the bits is what stops it acquiring a
/// difference nobody meant.
#[test]
fn dynamic_pca_with_no_lags_is_static_pca() {
    let mut rng = Normals::new(0x_D9C4);
    let data = rng.matrix(120, 4);
    let rule = Retention::CumulativeVariance(0.9);

    let static_model = Pca::fit(&data, 120, 4, rule);
    let dynamic = Dpca::fit(&data, 120, 4, 0, rule);

    assert_eq!(dynamic.lags(), 0);
    assert_eq!(dynamic.pca().retained(), static_model.retained());
    for (a, b) in dynamic
        .pca()
        .eigenvalues()
        .iter()
        .zip(static_model.eigenvalues())
    {
        assert_eq!(a.to_bits(), b.to_bits());
    }
    for (a, b) in dynamic.pca().loadings().iter().zip(static_model.loadings()) {
        assert_eq!(a.to_bits(), b.to_bits());
    }
    for t in 0..120 {
        let row = &data[t * 4..t * 4 + 4];
        assert_eq!(
            dynamic.t_squared(row).to_bits(),
            static_model.t_squared(row).to_bits()
        );
        assert_eq!(dynamic.spe(row).to_bits(), static_model.spe(row).to_bits());
    }
}

/// A window scored through `Dpca` is the corresponding row of the augmented
/// training matrix.
///
/// The window is given oldest first, the augmented row is present first, and
/// getting that backwards would be silent: the model would still fit, the
/// limits would still compute, and every statistic would be subtly wrong. So
/// the two paths are checked against each other on the training data itself,
/// where the answer is known.
#[test]
fn a_scored_window_is_the_augmented_row_it_should_be() {
    let mut rng = Normals::new(0x_147E);
    let (samples, variables, lags) = (200_usize, 3_usize, 2_usize);
    let data = rng.matrix(samples, variables);
    let dynamic = Dpca::fit(
        &data,
        samples,
        variables,
        lags,
        Retention::CumulativeVariance(0.9),
    );
    let (augmented, rows, columns) = augment_with_lags(&data, samples, variables, lags);

    for r in 0..rows {
        // The window in source order: samples r, r+1, r+2.
        let window = &data[r * variables..(r + lags + 1) * variables];
        assert_eq!(
            dynamic.augment_window(window),
            augmented[r * columns..(r + 1) * columns].to_vec(),
            "row {r}"
        );
        assert_eq!(
            dynamic.t_squared(window).to_bits(),
            dynamic
                .pca()
                .t_squared(&augmented[r * columns..(r + 1) * columns])
                .to_bits()
        );
    }

    // And the shape the limits are read at: the augmented row count, not the
    // record length.
    assert_eq!(dynamic.pca().samples(), samples - lags);
    assert_eq!(dynamic.pca().variables(), variables * (lags + 1));
    assert_eq!(dynamic.limits(0.99).samples, samples - lags);
}

/// Lag augmentation is what lets the model see serial structure, and here is a
/// process where that changes the answer.
///
/// Two variables driven by the same AR(1) factor, but the second lags the first
/// by one sample. A static model sees two variables that correlate weakly at
/// zero lag; a one-lag dynamic model sees an exact linear relation between
/// `x2(t)` and `x1(t-1)` and finds a near-zero eigenvalue for it.
///
/// The claim under test is the qualitative one Ku et al. make, stated as a
/// number: the dynamic model's smallest eigenvalue is far below the static
/// model's, because it has found a relation the static model cannot express.
#[test]
fn lag_augmentation_finds_a_relation_the_static_model_cannot_see() {
    let mut rng = Normals::new(0x_4610);
    let samples = 2000;
    let mut factor = 0.0;
    // Burn in so the series starts from its stationary distribution.
    for _ in 0..500 {
        factor = 0.9 * factor + rng.next();
    }
    let mut data = vec![0.0; samples * 2];
    let mut previous = factor;
    for t in 0..samples {
        factor = 0.9 * factor + rng.next();
        data[t * 2] = factor;
        // Exactly the previous sample of the first variable.
        data[t * 2 + 1] = previous;
        previous = factor;
    }

    let static_model = Pca::fit(&data, samples, 2, Retention::Fixed(2));
    let dynamic = Dpca::fit(&data, samples, 2, 1, Retention::Fixed(4));
    let static_smallest = static_model.eigenvalues()[1];
    let dynamic_smallest = *dynamic.pca().eigenvalues().last().expect("four components");
    println!(
        "static spectrum {:?}\ndynamic spectrum {:?}",
        static_model.eigenvalues(),
        dynamic.pca().eigenvalues()
    );
    assert!(
        dynamic_smallest < 1e-12,
        "the exact lag relation should give the dynamic model a null direction; \
         its smallest eigenvalue is {dynamic_smallest:.3e}"
    );
    assert!(
        static_smallest > 1e-3,
        "the static model should not see the relation at all; its smallest \
         eigenvalue is {static_smallest:.3e}, so this test does not demonstrate \
         anything"
    );
}

// ---------------------------------------------------------------------------
// Detection metrics
// ---------------------------------------------------------------------------

/// The three metrics on a series short enough to count by hand.
///
/// ```text
/// index   0  1  2  3  4  5  6  7  8  9
/// alarm   .  X  .  .  X  .  X  X  X  .
///                  ^ onset = 4
/// ```
///
/// Pre-fault: indices 0 to 3, one alarm, so FAR = 1/4 = 0.25.
/// Post-fault: indices 4 to 9, four alarms, so FDR = 4/6.
/// The first alarm after the onset is at index 4, a delay of 0. The first run of
/// three is at indices 6, 7, 8, a delay of 2. There is no run of four, so the
/// delay at `consecutive = 4` is `None`.
#[test]
fn the_detection_metrics_match_a_hand_counted_series() {
    let alarms = [
        false, true, false, false, true, false, true, true, true, false,
    ];
    let onset = 4;

    close(fault_detection_rate(&alarms, onset), 4.0 / 6.0, 0.0, "FDR");
    close(false_alarm_rate(&alarms, onset), 0.25, 0.0, "FAR");
    assert_eq!(detection_delay(&alarms, onset, 1), Some(0));
    assert_eq!(detection_delay(&alarms, onset, 2), Some(2));
    assert_eq!(detection_delay(&alarms, onset, 3), Some(2));
    assert_eq!(detection_delay(&alarms, onset, 4), None);

    let report = detection_report(&alarms, onset, 3);
    assert_eq!(
        report,
        DetectionReport {
            onset: 4,
            samples: 10,
            pre_fault: 4,
            post_fault: 6,
            false_alarms: 1,
            detections: 4,
            fault_detection_rate: 4.0 / 6.0,
            false_alarm_rate: 0.25,
            detection_delay: Some(2),
            consecutive: 3,
        }
    );
}

/// The boundary at `onset` belongs to the post-fault side, and nothing else
/// does.
///
/// An off-by-one here shifts every reported rate by one sample in a thousand,
/// which no test with a tolerance would catch and which would quietly bias a
/// cross-simulator comparison in whichever direction the fault happens to move
/// the statistic first.
#[test]
fn the_onset_sample_counts_as_post_fault() {
    // Exactly one alarm, at the onset itself.
    let mut alarms = vec![false; 10];
    alarms[4] = true;
    assert_eq!(fault_detection_rate(&alarms, 4), 1.0 / 6.0);
    assert_eq!(false_alarm_rate(&alarms, 4), 0.0);
    assert_eq!(detection_delay(&alarms, 4, 1), Some(0));

    // Move the onset one later and the same alarm becomes a false one.
    assert_eq!(fault_detection_rate(&alarms, 5), 0.0);
    assert_eq!(false_alarm_rate(&alarms, 5), 1.0 / 5.0);
    assert_eq!(detection_delay(&alarms, 5, 1), None);
}

#[test]
fn the_metrics_report_missing_data_as_nan_rather_than_zero() {
    let alarms = vec![true; 5];
    // No fault-free part.
    assert!(false_alarm_rate(&alarms, 0).is_nan());
    // No post-fault part.
    assert!(fault_detection_rate(&alarms, 5).is_nan());
    assert!(fault_detection_rate(&alarms, 99).is_nan());
    assert_eq!(detection_delay(&alarms, 5, 1), None);
    // An onset past the end still counts every sample as fault free, which is
    // the d00_te case: a whole record with no fault in it.
    assert_eq!(false_alarm_rate(&alarms, 99), 1.0);
    let report = detection_report(&alarms, 99, 3);
    assert_eq!(report.pre_fault, 5);
    assert_eq!(report.post_fault, 0);
    assert!(report.fault_detection_rate.is_nan());

    // A run that would not fit is not a detection.
    assert_eq!(detection_delay(&[false, true, true], 1, 3), None);
    assert_eq!(detection_delay(&[false, true, true], 1, 2), Some(0));
}

#[test]
#[should_panic(expected = "at least one alarm")]
fn a_detection_of_zero_consecutive_samples_is_rejected() {
    let _ = detection_delay(&[true, true], 0, 0);
}

/// Thresholding is strict, and `NaN` never alarms.
#[test]
fn alarms_are_raised_strictly_above_the_limit() {
    let statistic = [0.0, 1.0, 1.5, f64::NAN, f64::INFINITY, -2.0];
    assert_eq!(
        alarms_above(&statistic, 1.0),
        vec![false, false, true, false, true, false]
    );
}

// ---------------------------------------------------------------------------
// End to end
// ---------------------------------------------------------------------------

/// The whole chain, on the fault the two statistics were built to tell apart.
///
/// Six variables driven by two latent factors, so the correlation structure is
/// real and a two-component model explains almost all of it. The fault is a bias
/// on variable 5: a step of 1.5, against a residual standard deviation of 0.2.
///
/// That fault lies **in the residual subspace**. It moves variable 5 away from
/// the plane the other five pin it to, and hardly moves the plane at all, so
/// T-squared is nearly blind to it while SPE sees seven residual standard
/// deviations. Demonstrating that split is the point: it is the reason a PCA
/// monitor reports two numbers and not one, and a bug that quietly fed the same
/// projection into both would show up here and nowhere else in this file.
///
/// The false alarm rate on the fault-free part is reported against the nominal
/// 1% too. Its band is wide on purpose, because 200 samples cannot calibrate a
/// 1% rate;
/// `the_t_squared_limit_is_calibrated_at_a_sample_size_where_it_matters` does
/// that properly. Every number is printed, because a Tier 6 comparison is
/// between two runs of this and not against a constant.
#[test]
fn the_detector_separates_a_residual_fault_from_the_model_plane() {
    let mut rng = Normals::new(0x_E2E);
    let variables = 6;
    let build = |rng: &mut Normals, samples: usize, faulty: bool| {
        let mut data = vec![0.0; samples * variables];
        for t in 0..samples {
            let f1 = rng.next();
            let f2 = rng.next();
            for v in 0..3 {
                data[t * variables + v] = f1 + 0.2 * rng.next();
            }
            for v in 3..5 {
                data[t * variables + v] = f2 + 0.2 * rng.next();
            }
            // A sensor bias on the last variable: seven and a half times the
            // residual standard deviation, in a direction the two retained
            // components barely span.
            let bias = if faulty { 1.5 } else { 0.0 };
            data[t * variables + 5] = f2 + 0.2 * rng.next() + bias;
        }
        data
    };

    let training = build(&mut rng, 1000, false);
    let pca = Pca::fit(
        &training,
        1000,
        variables,
        Retention::CumulativeVariance(0.9),
    );
    let limits = pca.limits(0.99);
    println!(
        "model: {} of {variables} components, {:.1}% of the variance, spectrum {:?}",
        pca.retained(),
        pca.explained_variance() * 100.0,
        pca.eigenvalues()
    );
    println!(
        "limits at 99%: T^2 {:.4}, SPE {:.4}",
        limits.t_squared, limits.spe
    );

    // 200 fault-free samples then 800 faulty ones, the shape of a published TEP
    // test file.
    let mut testing = build(&mut rng, 200, false);
    testing.extend(build(&mut rng, 800, true));
    let onset = 200;

    let mut t2_series = Vec::with_capacity(1000);
    let mut spe_series = Vec::with_capacity(1000);
    for t in 0..1000 {
        let Statistics { t_squared, spe } =
            pca.statistics(&testing[t * variables..(t + 1) * variables]);
        t2_series.push(t_squared);
        spe_series.push(spe);
    }

    let t2_report = detection_report(&alarms_above(&t2_series, limits.t_squared), onset, 3);
    let spe_report = detection_report(&alarms_above(&spe_series, limits.spe), onset, 3);
    println!("T^2: {t2_report:?}");
    println!("SPE: {spe_report:?}");

    assert!(
        spe_report.fault_detection_rate > 0.95,
        "SPE detection rate {:.4} on a seven-sigma residual bias",
        spe_report.fault_detection_rate
    );
    assert_eq!(
        spe_report.detection_delay,
        Some(0),
        "a fault this large should be caught on the first three samples"
    );
    // The complementarity, stated as a number. T-squared is measured on the
    // same observations and barely moves, because the fault is orthogonal to
    // the plane it measures distance inside.
    assert!(
        t2_report.fault_detection_rate < 0.2,
        "T-squared detected {:.4} of a fault that lies in the residual \
         subspace, so this fixture does not separate the two statistics and the \
         test does not demonstrate what it claims",
        t2_report.fault_detection_rate
    );
    assert!(
        spe_report.fault_detection_rate > 4.0 * t2_report.fault_detection_rate.max(0.01),
        "SPE {:.4} against T-squared {:.4}",
        spe_report.fault_detection_rate,
        t2_report.fault_detection_rate
    );
    // The false alarm rate on the fault-free part, against the nominal 1%. 200
    // samples at 1% gives a binomial standard error of 0.7%, so this is a wide
    // band by construction: it is here to catch a limit that is out by an order
    // of magnitude, not to calibrate one.
    for report in [&t2_report, &spe_report] {
        assert!(
            report.false_alarm_rate < 0.05,
            "false alarm rate {:.4} against a nominal 0.01",
            report.false_alarm_rate
        );
    }
}
