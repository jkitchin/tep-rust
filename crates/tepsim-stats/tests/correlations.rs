//! Known-answer tests for the correlation matrix.

#![allow(
    clippy::float_cmp,
    reason = "known-answer tests assert exact values where the answer is exact"
)]
#![allow(
    clippy::suboptimal_flops,
    reason = "fixtures mirror the formulas they stand in for"
)]

use tepsim_stats::{CorrelationMatrix, frobenius_distance, worst_correlation_difference};

/// `x = [1,2,3,4]`, `y = [2,4,6,8]`, `z = [4,3,2,1]`, worked out by hand.
///
/// `y = 2x` exactly, so `r(x,y) = 1`. `z = 5 - x` exactly, so `r(x,z) = -1`
/// and `r(y,z) = -1`. These are exact, not approximate: a perfect linear
/// relationship has correlation exactly plus or minus one.
#[test]
fn perfect_linear_relationships_are_exactly_one() {
    let series = vec![
        vec![1.0, 2.0, 3.0, 4.0],
        vec![2.0, 4.0, 6.0, 8.0],
        vec![4.0, 3.0, 2.0, 1.0],
    ];
    let m = CorrelationMatrix::of(&series);
    assert_eq!(m.variables(), 3);
    for i in 0..3 {
        assert_eq!(m.get(i, i), 1.0, "diagonal at {i}");
    }
    assert_eq!(m.get(0, 1), 1.0);
    assert_eq!(m.get(0, 2), -1.0);
    assert_eq!(m.get(1, 2), -1.0);
    // Symmetric.
    for i in 0..3 {
        for j in 0..3 {
            assert_eq!(m.get(i, j), m.get(j, i), "({i}, {j})");
        }
    }
}

/// A hand-computed non-degenerate case.
///
/// `x = [1, 2, 3, 4]`, deviations `[-1.5, -0.5, 0.5, 1.5]`, `sum sq = 5`.
/// `w = [1, 3, 2, 4]`, mean 2.5, deviations `[-1.5, 0.5, -0.5, 1.5]`,
/// `sum sq = 5`.
///
/// ```text
/// dot = 2.25 - 0.25 - 0.25 + 2.25 = 4
/// r   = 4 / 5 = 0.8
/// ```
#[test]
fn a_hand_computed_correlation() {
    let series = vec![vec![1.0, 2.0, 3.0, 4.0], vec![1.0, 3.0, 2.0, 4.0]];
    let m = CorrelationMatrix::of(&series);
    assert!((m.get(0, 1) - 0.8).abs() < 1e-15, "r = {}", m.get(0, 1));
}

/// Correlation is invariant under a positive affine change of units, which is
/// the property that makes it comparable across variables with wildly
/// different scales, as the 53 in Tier 5 are.
#[test]
fn correlation_is_scale_and_offset_invariant() {
    let mut state = 0x0C0F_u64 | 1;
    let mut next = || {
        state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1);
        ((state >> 11) as f64) / ((1_u64 << 53) as f64) - 0.5
    };
    let a: Vec<f64> = (0..500).map(|_| next()).collect();
    let b: Vec<f64> = a.iter().map(|v| 0.7 * v + next() * 0.4).collect();

    let plain = CorrelationMatrix::of(&[a.clone(), b.clone()]);
    // Reactor pressure units against reactor level units, say.
    let rescaled = CorrelationMatrix::of(&[
        a.iter().map(|v| 2705.0 + 61.0 * v).collect(),
        b.iter().map(|v| 0.003 - 1e-4 * v).collect(),
    ]);
    // The second is negated, so the correlation flips sign.
    assert!(
        (plain.get(0, 1) + rescaled.get(0, 1)).abs() < 1e-12,
        "{} against {}",
        plain.get(0, 1),
        rescaled.get(0, 1)
    );
}

/// A constant variable correlates with nothing, and says so rather than
/// claiming zero.
#[test]
fn a_constant_variable_is_nan_off_the_diagonal() {
    let series = vec![vec![1.0, 2.0, 3.0, 4.0], vec![7.0, 7.0, 7.0, 7.0]];
    let m = CorrelationMatrix::of(&series);
    assert_eq!(m.get(0, 0), 1.0);
    assert_eq!(m.get(1, 1), 1.0, "a constant is still itself");
    assert!(m.get(0, 1).is_nan());
    assert!(m.get(1, 0).is_nan());
}

/// Every entry is a correlation, on many random matrices.
#[test]
fn every_entry_is_in_range() {
    let mut state = 0x1234_u64 | 1;
    let mut next = || {
        state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1);
        ((state >> 11) as f64) / ((1_u64 << 53) as f64) - 0.5
    };
    for _ in 0..20 {
        let series: Vec<Vec<f64>> = (0..8).map(|_| (0..200).map(|_| next()).collect()).collect();
        let m = CorrelationMatrix::of(&series);
        for i in 0..8 {
            for j in 0..8 {
                let r = m.get(i, j);
                assert!((-1.0..=1.0).contains(&r), "({i}, {j}) = {r}");
            }
        }
    }
}

/// Near-collinear data is where the ratio can exceed one by an ulp.
#[test]
fn near_collinear_data_does_not_produce_a_correlation_above_one() {
    let a: Vec<f64> = (0..1000).map(|i| 1e6 + i as f64 * 1e-6).collect();
    let b: Vec<f64> = a.iter().map(|v| v * 3.0 + 1.0).collect();
    let m = CorrelationMatrix::of(&[a, b]);
    assert!(m.get(0, 1) <= 1.0, "r = {:.20}", m.get(0, 1));
    assert!(m.get(0, 1) > 0.999_999, "r = {:.20}", m.get(0, 1));
}

// ---------------------------------------------------------------------------
// Frobenius distance
// ---------------------------------------------------------------------------

#[test]
fn a_matrix_is_zero_distance_from_itself() {
    let series = vec![
        vec![1.0, 2.0, 3.0, 4.0, 6.0],
        vec![1.0, 3.0, 2.0, 4.0, 5.0],
        vec![5.0, 1.0, 4.0, 2.0, 3.0],
    ];
    let m = CorrelationMatrix::of(&series);
    let (distance, skipped) = frobenius_distance(&m, &m);
    assert_eq!(distance, 0.0);
    assert_eq!(skipped, 0);
    assert!(worst_correlation_difference(&m, &m).is_some_and(|(_, _, a, b)| a == b));
}

/// A hand-computed distance.
///
/// Two 2-by-2 matrices differing only in the off-diagonal, by `d` in each of
/// the two symmetric positions, are `sqrt(2 d^2) = d sqrt(2)` apart.
#[test]
fn the_frobenius_distance_matches_a_hand_computed_case() {
    // r = 1 for the first pair, r = 0.8 for the second (from the test above).
    let a = CorrelationMatrix::of(&[vec![1.0, 2.0, 3.0, 4.0], vec![2.0, 4.0, 6.0, 8.0]]);
    let b = CorrelationMatrix::of(&[vec![1.0, 2.0, 3.0, 4.0], vec![1.0, 3.0, 2.0, 4.0]]);
    let (distance, skipped) = frobenius_distance(&a, &b);
    let expected = 0.2 * 2.0_f64.sqrt();
    assert!(
        (distance - expected).abs() < 1e-15,
        "distance {distance}, expected {expected}"
    );
    assert_eq!(skipped, 0);

    let (i, j, x, y) = worst_correlation_difference(&a, &b).expect("a difference");
    assert_eq!((i, j), (0, 1));
    assert_eq!(x, 1.0);
    assert!((y - 0.8).abs() < 1e-15);
}

#[test]
fn nan_entries_are_skipped_and_counted() {
    let with_constant = CorrelationMatrix::of(&[
        vec![1.0, 2.0, 3.0, 4.0],
        vec![7.0, 7.0, 7.0, 7.0],
        vec![4.0, 1.0, 3.0, 2.0],
    ]);
    let (distance, skipped) = frobenius_distance(&with_constant, &with_constant);
    assert_eq!(distance, 0.0);
    // The constant variable's row and column, off the diagonal: four entries.
    assert_eq!(skipped, 4);
    // And the worst difference ignores them rather than returning NaN.
    let (_, _, x, y) = worst_correlation_difference(&with_constant, &with_constant)
        .expect("a finite pair remains");
    assert!(x.is_finite() && y.is_finite());
}

#[test]
#[should_panic(expected = "cannot compare")]
fn matrices_of_different_sizes_are_rejected() {
    let a = CorrelationMatrix::of(&[vec![1.0, 2.0], vec![2.0, 1.0]]);
    let b = CorrelationMatrix::of(&[vec![1.0, 2.0]]);
    let _ = frobenius_distance(&a, &b);
}
