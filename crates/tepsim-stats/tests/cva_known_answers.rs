//! Known-answer tests for the generalised eigensolver and canonical variate
//! analysis.
//!
//! The rule `known_answers.rs` and `eigen_known_answers.rs` set: **nothing is
//! checked against a number this project produced**. Every assertion is
//! against a closed form, an exact identity, a construction whose answer is a
//! rational number known before the code runs, or the standard eigensolver
//! that already has its own battery.
//!
//! A generalised eigensolver has failure modes the standard one does not, so
//! the layers are chosen to separate them:
//!
//! - **Cholesky against a textbook factorisation.** Integer in, integer out,
//!   so the whole factorisation is exact and a transposed index or a dropped
//!   term shows up as a mismatch and not as a rounding difference.
//! - **`B = I` reproduces the standard problem bit for bit.** This pins the
//!   plumbing: the reduction, the back-transform and the sign convention all
//!   have to be the identity when `L` is.
//! - **A diagonal `B` gives the known rescaling, exactly.** Powers of four, so
//!   every square root and every division is exact in binary floating point
//!   and the expected answer is a closed form with no error term at all. This
//!   is the layer that catches `L^-1 A L^-T` written as `L^-T A L^-1`.
//! - **A 2-by-2 pencil against the quadratic `det(A - lambda B) = 0`.** A
//!   closed form for a *non-diagonal* `B`, which none of the above is.
//! - **Invariants.** `sum lambda = trace(B^-1 A)`, `prod lambda =
//!   det A / det B`, with both determinants from the 3-by-3 cofactor formula
//!   written out here. And `A V = B V Lambda` with `V' B V = I`, which needs
//!   no reference at all.
//! - **Canonical correlations of a construction built from Pythagorean
//!   triples.** The answers are 3/5 and 5/13, known before the test runs, and
//!   every intermediate quantity is a small integer or an exact square root,
//!   so the assertion is bit-exact rather than approximate.

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

use tepsim_stats::cva::{Cva, cholesky, generalized_symmetric_eigen, past_future};
use tepsim_stats::symmetric_eigen;

/// The same deterministic generator `eigen_known_answers.rs` uses, so these
/// tests need no random-number crate and are reproducible.
struct Lcg(u64);

impl Lcg {
    fn new(seed: u64) -> Self {
        Self(seed | 1)
    }
    fn next(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1);
        ((self.0 >> 11) as f64) / ((1_u64 << 53) as f64) - 0.5
    }
}

fn random_symmetric(rng: &mut Lcg, n: usize, scale: f64) -> Vec<f64> {
    let mut a = vec![0.0; n * n];
    for i in 0..n {
        for j in i..n {
            let value = scale * rng.next();
            a[i * n + j] = value;
            a[j * n + i] = value;
        }
    }
    a
}

/// A symmetric positive definite matrix: `M M' + n I`, which is definite for
/// any `M` because the shift is strictly positive.
fn random_definite(rng: &mut Lcg, n: usize, scale: f64) -> Vec<f64> {
    let mut m = vec![0.0; n * n];
    for value in &mut m {
        *value = scale * rng.next();
    }
    let mut b = vec![0.0; n * n];
    for i in 0..n {
        for j in i..n {
            let mut total = 0.0;
            for k in 0..n {
                total += m[i * n + k] * m[j * n + k];
            }
            if i == j {
                total += n as f64;
            }
            // Written into both triangles from one sum, so the result is
            // exactly symmetric and the solver's strict check has nothing to
            // trip on.
            b[i * n + j] = total;
            b[j * n + i] = total;
        }
    }
    b
}

/// The 16-by-16 Sylvester Hadamard matrix, as rows of `+-1`.
///
/// Rows 1 through 15 sum to exactly zero and are mutually orthogonal with
/// `h_i . h_i = 16`. Every entry, every inner product and every square root of
/// an inner product is exact in `f64`, which is what makes the canonical
/// correlation test below exact rather than approximate.
fn hadamard16() -> Vec<Vec<f64>> {
    let mut h = vec![vec![1.0; 16]; 16];
    for i in 0..16 {
        for j in 0..16 {
            // The Sylvester construction: the sign is the parity of the number
            // of bit positions where `i` and `j` are both one.
            h[i][j] = if (i & j).count_ones() % 2 == 0 {
                1.0
            } else {
                -1.0
            };
        }
    }
    h
}

// ---------------------------------------------------------------------------
// Cholesky
// ---------------------------------------------------------------------------

/// The textbook 3-by-3 factorisation, exact.
///
/// ```text
/// [   4   12  -16 ]     [  2  0  0 ]
/// [  12   37  -43 ]  =  [  6  1  0 ] L'
/// [ -16  -43   98 ]     [ -8  5  3 ]
/// ```
///
/// This is the worked example every numerical linear algebra text and the
/// Cholesky article carry, and it is chosen here because every entry of both
/// the input and the answer is a small integer: the factorisation is exact in
/// `f64`, so the assertion is equality and not a tolerance.
#[test]
fn the_textbook_cholesky_factorisation_is_exact() {
    let b = [
        4.0, 12.0, -16.0, //
        12.0, 37.0, -43.0, //
        -16.0, -43.0, 98.0,
    ];
    let expected = [
        2.0, 0.0, 0.0, //
        6.0, 1.0, 0.0, //
        -8.0, 5.0, 3.0,
    ];
    let l = cholesky(&b, 3).expect("the textbook matrix is positive definite");
    assert_eq!(l, expected.to_vec(), "got {l:?}");

    // And `L L'` is the input again, exactly, since all three are integers.
    for i in 0..3 {
        for j in 0..3 {
            let mut total = 0.0;
            for k in 0..3 {
                total += l[i * 3 + k] * l[j * 3 + k];
            }
            assert_eq!(total, b[i * 3 + j], "L L' differs at ({i}, {j})");
        }
    }
}

/// A matrix that is not positive definite gets `None`, not a repaired answer.
#[test]
fn an_indefinite_matrix_has_no_factorisation() {
    // Eigenvalues 1 and -1: symmetric, not definite.
    assert!(cholesky(&[0.0, 1.0, 1.0, 0.0], 2).is_none());
    // Singular: the second pivot is exactly zero, which is not positive.
    assert!(cholesky(&[1.0, 1.0, 1.0, 1.0], 2).is_none());
    // Negative on the diagonal.
    assert!(cholesky(&[-1.0], 1).is_none());
    // A `NaN` pivot fails rather than propagating.
    assert!(cholesky(&[f64::NAN], 1).is_none());
    // And the identity factors into itself, exactly.
    let i3 = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
    assert_eq!(cholesky(&i3, 3).expect("the identity is definite"), i3);
}

/// `L L' = B` over matrices with no closed form, to machine precision.
#[test]
fn the_factorisation_reconstructs_its_input() {
    let mut rng = Lcg::new(0xC401E5);
    let mut worst = 0.0_f64;
    for n in [1_usize, 2, 5, 11, 24] {
        let b = random_definite(&mut rng, n, 2.0);
        let l = cholesky(&b, n).expect("M M' + n I is positive definite");
        let scale = b.iter().map(|x| x.abs()).fold(0.0_f64, f64::max);
        for i in 0..n {
            for j in 0..n {
                let mut total = 0.0;
                for k in 0..n {
                    total += l[i * n + k] * l[j * n + k];
                }
                worst = worst.max((total - b[i * n + j]).abs() / scale);
            }
            // Strictly upper entries are exactly zero, not nearly.
            for j in (i + 1)..n {
                assert_eq!(l[i * n + j], 0.0, "L is not lower triangular at ({i}, {j})");
            }
        }
    }
    println!("Cholesky: worst scaled ||L L' - B||_max = {worst:.3e}");
    assert!(worst < 1e-14, "{worst:.3e}");
}

// ---------------------------------------------------------------------------
// The generalised problem, against closed forms
// ---------------------------------------------------------------------------

/// `B = I` gives the standard problem's answer, **bit for bit**.
///
/// With `B = I` the Cholesky factor is the identity, so the reduction and the
/// back-transform are both no-ops and there is nothing for the generalised
/// path to do that the standard path did not. Anything less than bit equality
/// would mean a stray rounding step somewhere in the reduction, most likely a
/// normalisation. This crate deliberately does not renormalise the
/// back-transformed vectors, precisely so that this identity is exact and can
/// be asserted as one.
#[test]
fn a_generalised_problem_with_b_equal_to_the_identity_is_the_standard_one() {
    let mut rng = Lcg::new(0x1DE7717);
    for n in [1_usize, 2, 3, 7, 16] {
        let a = random_symmetric(&mut rng, n, 6.0);
        let mut identity = vec![0.0; n * n];
        for i in 0..n {
            identity[i * n + i] = 1.0;
        }

        let standard = symmetric_eigen(&a, n);
        let generalised =
            generalized_symmetric_eigen(&a, &identity, n).expect("the identity is definite");

        assert_eq!(generalised.sweeps(), standard.sweeps(), "n={n}");
        assert_eq!(
            generalised.off_diagonal_norm().to_bits(),
            standard.off_diagonal_norm().to_bits(),
            "n={n}"
        );
        for j in 0..n {
            assert_eq!(
                generalised.values()[j].to_bits(),
                standard.values()[j].to_bits(),
                "n={n}: eigenvalue {j} is {} against {}",
                generalised.values()[j],
                standard.values()[j]
            );
            for i in 0..n {
                assert_eq!(
                    generalised.component(i, j).to_bits(),
                    standard.component(i, j).to_bits(),
                    "n={n}: V[{i}][{j}] is {} against {}",
                    generalised.component(i, j),
                    standard.component(i, j)
                );
            }
        }
        // The Cholesky factor is the identity, so every pivot is one.
        assert_eq!(generalised.smallest_pivot(), 1.0);
        assert_eq!(generalised.largest_pivot(), 1.0);
        println!(
            "B = I at n={n}: {} sweeps, identical bit for bit, \
             ||V'BV - I||_max = {:.3e}",
            generalised.sweeps(),
            generalised.b_orthonormality()
        );
    }
}

/// A diagonal `B` gives the known rescaling, exactly.
///
/// For `A = diag(a)` and `B = diag(b)` the pencil is `a_i v = lambda b_i v`
/// componentwise, so
///
/// ```text
/// lambda_i = a_i / b_i,      v_i = e_i / sqrt(b_i)
/// ```
///
/// with `V' B V = I`. Every `b_i` here is an even power of two, so `sqrt(b_i)`
/// and its reciprocal are exact, and every `a_i / b_i` is a ratio of powers of
/// two and therefore exact too. The expected answer carries no error term at
/// all, so the assertion is equality.
///
/// This is the layer that catches the reduction written the wrong way round.
/// `L^-T A L^-1` is also symmetric and also has real eigenvalues; with a
/// diagonal `B` it produces `a_i / b_i` as well, so this test alone would not
/// separate them. What it does separate is `L^-1 A L^-T` from any form that
/// divides by `b_i` once rather than twice, which is the more common slip and
/// which shows up here as eigenvalues out by a factor of `b_i`.
#[test]
fn a_diagonal_pencil_returns_the_exact_rescaling() {
    // Powers of two throughout, `b` at even exponents so its square roots are
    // exact: 4 = 2^2, 16 = 2^4, 1 = 2^0, 64 = 2^6.
    let a = [8.0, 2.0, 6.0, 4.0];
    let b = [4.0, 16.0, 1.0, 64.0];
    let n = 4;
    let mut am = vec![0.0; n * n];
    let mut bm = vec![0.0; n * n];
    for i in 0..n {
        am[i * n + i] = a[i];
        bm[i * n + i] = b[i];
    }

    let eigen = generalized_symmetric_eigen(&am, &bm, n).expect("a positive diagonal is definite");
    // a/b is 2, 0.125, 6, 0.0625; descending that is 6, 2, 0.125, 0.0625, from
    // original positions 2, 0, 1, 3.
    let expected_values = [6.0, 2.0, 0.125, 0.0625];
    let expected_source = [2_usize, 0, 1, 3];
    assert_eq!(eigen.values(), &expected_values, "got {:?}", eigen.values());
    for (j, &source) in expected_source.iter().enumerate() {
        for i in 0..n {
            let expected = if i == source { 1.0 / b[i].sqrt() } else { 0.0 };
            assert_eq!(
                eigen.component(i, j),
                expected,
                "V[{i}][{j}] for eigenvalue {}",
                eigen.values()[j]
            );
        }
    }
    // Diagonal in, diagonal out: the reduced matrix is already diagonal, so no
    // rotation should fire at all.
    assert_eq!(eigen.sweeps(), 0);
    assert_eq!(eigen.b_orthonormality(), 0.0);
    println!(
        "diagonal pencil: eigenvalues {:?} exactly, pivots {} to {}",
        eigen.values(),
        eigen.smallest_pivot(),
        eigen.largest_pivot()
    );
}

/// `B = 4 I` divides every eigenvalue by four and halves every eigenvector,
/// exactly.
///
/// Scaling `B` by a power of two scales `C = L^-1 A L^-T` by the same power,
/// and every step of a Jacobi rotation is invariant under an exact power-of-two
/// scaling of its input: the rotation angle comes from a ratio, and the
/// updates are all products with the scaled entries. So the answer is the
/// standard problem's, divided by four, with no error term. A non-diagonal `A`
/// here, so this is not the previous test again.
#[test]
fn scaling_b_by_a_power_of_two_scales_the_spectrum_exactly() {
    let mut rng = Lcg::new(0x5CA1E);
    for n in [2_usize, 5, 9] {
        let a = random_symmetric(&mut rng, n, 3.0);
        let mut four_i = vec![0.0; n * n];
        for i in 0..n {
            four_i[i * n + i] = 4.0;
        }
        let standard = symmetric_eigen(&a, n);
        let scaled = generalized_symmetric_eigen(&a, &four_i, n).expect("4I is definite");
        for j in 0..n {
            assert_eq!(
                scaled.values()[j],
                standard.values()[j] / 4.0,
                "n={n}, eigenvalue {j}"
            );
            for i in 0..n {
                assert_eq!(
                    scaled.component(i, j),
                    standard.component(i, j) / 2.0,
                    "n={n}, V[{i}][{j}]"
                );
            }
        }
        assert_eq!(scaled.smallest_pivot(), 2.0);
        println!("B = 4I at n={n}: spectrum divided by four, exactly");
    }
}

/// The 2-by-2 pencil against `det(A - lambda B) = 0`.
///
/// Expanding the determinant gives a quadratic in `lambda`,
///
/// ```text
/// (b11 b22 - b12^2) lambda^2
///   - (a11 b22 + a22 b11 - 2 a12 b12) lambda
///   + (a11 a22 - a12^2) = 0
/// ```
///
/// which is a closed form for a **non-diagonal** `B`, the case every other
/// exact test above avoids.
///
/// # The discriminant has to be written the long way, and finding that out
/// cost a failing run
///
/// The stable quadratic formula, `q = -(c + sign(c) sqrt(c^2 - 4 p r)) / 2`
/// with roots `q / p` and `r / q`, protects the *small root* from
/// cancellation. It does nothing for the discriminant itself, and here the
/// discriminant is where the cancellation is: at `a11 = a22 = 7`,
/// `a12 = -1e-6`, `b12 = 0`, `c^2` is 196 and `4 p r` is `196 - 4e-12`, so
/// `c^2 - 4 p r` keeps four significant digits of a quantity that should have
/// sixteen. The reference then reported eigenvalues 9.3e-10 away from the
/// truth while the solver had them exactly right, and the test failed on the
/// reference.
///
/// With `b11 = b22 = 1` the discriminant expands and regroups into
///
/// ```text
/// (a11 - a22)^2 + 4 (a12 - b12 a11)(a12 - b12 a22)
/// ```
///
/// which is the form used below. The `a12^2 b12^2` terms cancel algebraically
/// rather than numerically, and on the case above it evaluates to `4e-12`
/// exactly. This is the same lesson `eigen_known_answers.rs` records about the
/// textbook eigenvector formula: a closed form is only a reference if it is
/// the accurate closed form.
#[test]
fn the_two_by_two_pencil_matches_the_determinant_quadratic() {
    let mut worst = 0.0_f64;
    let mut worst_case = None;
    let mut checked = 0;
    for a11 in [-3.0_f64, -0.25, 0.0, 1.0, 7.0] {
        for a22 in [-3.0_f64, -0.25, 0.0, 1.0, 7.0] {
            for a12 in [-2.0_f64, -1e-6, 0.5, 3.0] {
                for b12 in [-0.9_f64, -0.1, 0.0, 0.3, 0.95] {
                    // `B = [[1, b12], [b12, 1]]` is positive definite exactly
                    // when `|b12| < 1`, so every case here has a factorisation.
                    let a = [a11, a12, a12, a22];
                    let b = [1.0, b12, b12, 1.0];
                    let eigen =
                        generalized_symmetric_eigen(&a, &b, 2).expect("|b12| < 1 is definite");

                    let p = 1.0 - b12 * b12;
                    let c = -(a11 + a22 - 2.0 * a12 * b12);
                    let r = a11 * a22 - a12 * a12;
                    // The regrouped discriminant. See the doc comment: the
                    // literal `c * c - 4 p r` cancels to four digits here.
                    let discriminant =
                        (a11 - a22) * (a11 - a22) + 4.0 * (a12 - b12 * a11) * (a12 - b12 * a22);
                    assert!(
                        discriminant >= 0.0,
                        "a definite pencil has real eigenvalues; got {discriminant}"
                    );
                    let q = -0.5 * (c + c.signum() * discriminant.sqrt());
                    // Both roots, from `q/p` and `r/q`, so neither is formed by
                    // a subtraction of nearly equal numbers. `q` is zero only
                    // when `c` and the discriminant both are, which needs
                    // `r = 0` too and leaves a double root at zero.
                    let (mut low, mut high) = if q == 0.0 { (0.0, 0.0) } else { (q / p, r / q) };
                    if low > high {
                        core::mem::swap(&mut low, &mut high);
                    }

                    let scale = high.abs().max(low.abs()).max(1.0);
                    let here = ((eigen.values()[0] - high).abs())
                        .max((eigen.values()[1] - low).abs())
                        / scale;
                    if here > worst {
                        worst = here;
                        worst_case = Some((a11, a22, a12, b12));
                    }
                    checked += 2;
                }
            }
        }
    }
    println!(
        "2x2 pencil: {checked} eigenvalues against det(A - lambda B) = 0, \
         worst scaled error {worst:.3e} at {worst_case:?}"
    );
    assert!(worst < 1e-13, "{worst:.3e} at {worst_case:?}");
}

/// `sum lambda = trace(B^-1 A)` and `prod lambda = det A / det B`.
///
/// Both determinants come from the 3-by-3 cofactor formula written out here,
/// and the trace of `B^-1 A` from Cramer's rule on the same cofactors, so
/// neither side of either identity touches the solver.
#[test]
fn the_generalised_spectrum_reproduces_the_trace_and_the_determinant() {
    let mut rng = Lcg::new(0x07AC_EDE7);
    let mut worst_trace = 0.0_f64;
    let mut worst_determinant = 0.0_f64;
    for _ in 0..200 {
        let a = random_symmetric(&mut rng, 3, 4.0);
        let b = random_definite(&mut rng, 3, 1.0);
        let eigen = generalized_symmetric_eigen(&a, &b, 3).expect("M M' + 3I is definite");

        let det = |m: &[f64]| -> f64 {
            m[0] * (m[4] * m[8] - m[5] * m[7]) - m[1] * (m[3] * m[8] - m[5] * m[6])
                + m[2] * (m[3] * m[7] - m[4] * m[6])
        };
        // `trace(B^-1 A) = sum_ij cofactor(B)_ij A_ij / det B`, which is
        // Cramer's rule and needs no matrix inverse.
        let cofactor = |m: &[f64], i: usize, j: usize| -> f64 {
            let rows: Vec<usize> = (0..3).filter(|r| *r != i).collect();
            let cols: Vec<usize> = (0..3).filter(|c| *c != j).collect();
            let minor = m[rows[0] * 3 + cols[0]] * m[rows[1] * 3 + cols[1]]
                - m[rows[0] * 3 + cols[1]] * m[rows[1] * 3 + cols[0]];
            if (i + j) % 2 == 0 { minor } else { -minor }
        };
        let det_b = det(&b);
        let mut trace = 0.0;
        for i in 0..3 {
            for j in 0..3 {
                // `(B^-1)_ij = cofactor(B)_ji / det B`, and
                // `trace(B^-1 A) = sum_ij (B^-1)_ij A_ji`.
                trace += cofactor(&b, j, i) * a[j * 3 + i] / det_b;
            }
        }

        let sum: f64 = eigen.values().iter().sum();
        let product: f64 = eigen.values().iter().product();
        let scale = eigen
            .values()
            .iter()
            .map(|v| v.abs())
            .fold(1.0_f64, f64::max);
        worst_trace = worst_trace.max((sum - trace).abs() / (scale * 3.0));
        worst_determinant =
            worst_determinant.max((product - det(&a) / det_b).abs() / scale.powi(3));
    }
    println!(
        "generalised trace and determinant over 200 random 3x3 pencils: \
         worst scaled trace error {worst_trace:.3e}, worst scaled determinant \
         error {worst_determinant:.3e}"
    );
    assert!(worst_trace < 1e-13, "{worst_trace:.3e}");
    assert!(worst_determinant < 1e-13, "{worst_determinant:.3e}");
}

/// `A V = B V Lambda` and `V' B V = I`, over pencils with no closed form.
///
/// The defining equation of the problem, checkable to machine precision with
/// no reference at all. The `B`-orthonormality is the quantity the module
/// declines to enforce by rescaling, so this is where its size is established
/// rather than assumed.
#[test]
fn the_pencil_residual_and_b_orthonormality_stay_small() {
    let mut rng = Lcg::new(0xBE9C1);
    for n in [1_usize, 2, 4, 12, 25] {
        let a = random_symmetric(&mut rng, n, 5.0);
        let b = random_definite(&mut rng, n, 1.5);
        let eigen = generalized_symmetric_eigen(&a, &b, n).expect("M M' + nI is definite");
        assert!(eigen.converged(), "n={n}");

        let scale = a.iter().map(|x| x.abs()).fold(0.0_f64, f64::max)
            * eigen
                .values()
                .iter()
                .map(|v| v.abs())
                .fold(1.0_f64, f64::max);
        let mut residual = 0.0_f64;
        for i in 0..n {
            for j in 0..n {
                let mut av = 0.0;
                let mut bv = 0.0;
                for k in 0..n {
                    av += a[i * n + k] * eigen.component(k, j);
                    bv += b[i * n + k] * eigen.component(k, j);
                }
                residual = residual.max((av - eigen.values()[j] * bv).abs() / scale);
            }
        }
        // `V' B V - I`, recomputed here rather than read off the struct, and
        // then compared against what the struct reports. Asserting only that
        // the reported number is small would pass for a diagnostic stuck at
        // zero, which is the more likely way for it to be wrong: a residual
        // nobody computes reads as a perfect residual. A mutation run confirmed
        // it, so this checks the value as well as the bound.
        let mut orthonormality = 0.0_f64;
        for i in 0..n {
            for j in 0..n {
                let mut total = 0.0;
                for k in 0..n {
                    let mut bv = 0.0;
                    for m in 0..n {
                        bv += b[k * n + m] * eigen.component(m, j);
                    }
                    total += eigen.component(k, i) * bv;
                }
                let expected = if i == j { 1.0 } else { 0.0 };
                orthonormality = orthonormality.max((total - expected).abs());
            }
        }
        println!(
            "pencil n={n}: {} sweeps, ||AV - BV L||_max / scale = {residual:.3e}, \
             ||V'BV - I||_max = {:.3e} (recomputed {orthonormality:.3e}), \
             pivots {:.3e} to {:.3e}",
            eigen.sweeps(),
            eigen.b_orthonormality(),
            eigen.smallest_pivot(),
            eigen.largest_pivot()
        );
        assert!(residual < 1e-13, "n={n}: {residual:.3e}");
        assert!(
            eigen.b_orthonormality() < 1e-13,
            "n={n}: {:.3e}",
            eigen.b_orthonormality()
        );
        assert!(
            (eigen.b_orthonormality() - orthonormality).abs() < 1e-16,
            "n={n}: the reported B-orthonormality {:.3e} is not the recomputed \
             {orthonormality:.3e}",
            eigen.b_orthonormality()
        );
        if n > 1 {
            // On random data the deviation is genuinely nonzero, so a stuck
            // zero is a failure and not a lucky exact answer.
            assert!(
                orthonormality > 0.0,
                "n={n}: V'BV came out exactly the identity, which on random \
                 data means it was not computed"
            );
        }
        for pair in eigen.values().windows(2) {
            assert!(pair[0] >= pair[1], "n={n}: not descending");
        }
    }
}

/// An indefinite `B` yields `None` from the pencil solver, not a spectrum.
#[test]
fn a_pencil_with_an_indefinite_b_declines_to_answer() {
    let a = [2.0, 1.0, 1.0, 3.0];
    assert!(generalized_symmetric_eigen(&a, &[0.0, 1.0, 1.0, 0.0], 2).is_none());
    // Singular `B`: the second pivot is exactly zero.
    assert!(generalized_symmetric_eigen(&a, &[1.0, 1.0, 1.0, 1.0], 2).is_none());
}

#[test]
#[should_panic(expected = "not symmetric")]
fn an_asymmetric_pencil_is_rejected() {
    let _ = generalized_symmetric_eigen(&[1.0, 2.0, 3.0, 4.0], &[1.0, 0.0, 0.0, 1.0], 2);
}

/// The sign convention survives the back-transform.
///
/// [`crate::symmetric_eigen`] imposes "largest-magnitude component positive" on
/// the reduced problem's vectors, and `L^-T` can move which component is
/// largest, so the convention has to be re-imposed afterwards. Nothing else in
/// this file notices if it is not: a mutation that dropped the re-imposition
/// left every other test passing, because every other test compares vectors up
/// to sign or uses a `B` diagonal enough that the leader cannot move.
#[test]
fn the_largest_component_of_every_generalised_eigenvector_is_positive() {
    let mut rng = Lcg::new(0x0005_169E);
    for n in [1_usize, 2, 4, 9, 17] {
        let a = random_symmetric(&mut rng, n, 3.0);
        let b = random_definite(&mut rng, n, 4.0);
        let eigen = generalized_symmetric_eigen(&a, &b, n).expect("definite");
        for j in 0..n {
            let vector = eigen.eigenvector(j);
            let mut leader = 0;
            for (i, value) in vector.iter().enumerate() {
                if value.abs() > vector[leader].abs() {
                    leader = i;
                }
            }
            assert!(
                vector[leader] > 0.0,
                "n={n} vector {j} peaks at index {leader} with {}",
                vector[leader]
            );
        }
    }
}

/// The same input gives the same bits.
#[test]
fn the_generalised_decomposition_is_reproducible_bit_for_bit() {
    let mut rng = Lcg::new(0xDE7E12);
    let a = random_symmetric(&mut rng, 9, 2.0);
    let b = random_definite(&mut rng, 9, 1.0);
    let first = generalized_symmetric_eigen(&a, &b, 9).expect("definite");
    let second = generalized_symmetric_eigen(&a, &b, 9).expect("definite");
    assert_eq!(first, second);
}

// ---------------------------------------------------------------------------
// Canonical correlations
// ---------------------------------------------------------------------------

/// A construction whose canonical correlations are 3/5 and 5/13, exactly.
///
/// Take four mutually orthogonal `+-1` vectors of length sixteen from the
/// Sylvester Hadamard matrix, `h1` through `h4`, each with `h . h = 16` and
/// each summing to zero so that centring is a no-op. Set
///
/// ```text
/// X = [ h1, h2 ]          Y = [ 3 h1 + 4 h3,  5 h2 + 12 h4 ]
/// ```
///
/// Then every cross-product matrix is diagonal:
/// `Sxx = diag(16, 16)`, `Syy = diag(400, 2704)`, `Sxy = diag(48, 80)`, and
/// the canonical correlations are `48 / sqrt(16 * 400) = 3/5` and
/// `80 / sqrt(16 * 2704) = 5/13`.
///
/// The numbers are chosen so that nothing rounds. The Pythagorean triples
/// 3-4-5 and 5-12-13 make `sqrt(400) = 20` and `sqrt(2704) = 52` exact, so the
/// Cholesky factor is `diag(4, 4, 20, 52)` exactly, the reduced matrix has
/// `48 / 80 = 3/5` in it exactly, and the single Jacobi rotation on a
/// `[[0, c], [c, 0]]` block returns `+-c` exactly. So this asserts equality
/// rather than a tolerance, and 3/5 and 5/13 were known before the code ran.
///
/// It also pins the `+rho, -rho` pairing: the four eigenvalues have to be
/// `3/5, 5/13, -5/13, -3/5`, and nothing in the implementation enforces that.
#[test]
fn the_canonical_correlations_of_a_pythagorean_construction_are_exact() {
    let h = hadamard16();
    let samples = 16;
    let mut x = vec![0.0; samples * 2];
    let mut y = vec![0.0; samples * 2];
    for t in 0..samples {
        x[t * 2] = h[1][t];
        x[t * 2 + 1] = h[2][t];
        y[t * 2] = 3.0 * h[1][t] + 4.0 * h[3][t];
        y[t * 2 + 1] = 5.0 * h[2][t] + 12.0 * h[4][t];
    }

    let cva = Cva::fit(&x, &y, samples, 2, 2).expect("the blocks are full rank");
    assert_eq!(cva.pairs(), 2);
    assert_eq!(
        cva.correlations(),
        &[3.0 / 5.0, 5.0 / 13.0],
        "got {:?}",
        cva.correlations()
    );
    // Unclamped too, so the clamp is not what is being tested here.
    assert_eq!(cva.raw_correlations(), &[3.0 / 5.0, 5.0 / 13.0]);
    // The spectrum is `+-rho`, exactly.
    assert_eq!(
        cva.eigen().values(),
        &[3.0 / 5.0, 5.0 / 13.0, -(5.0 / 13.0), -(3.0 / 5.0)]
    );
    assert_eq!(cva.spectrum_symmetry(), 0.0);
    println!(
        "Pythagorean construction: canonical correlations {:?}, exactly 3/5 and \
         5/13; spectrum symmetry {:.3e}, ||V'BV - I||_max = {:.3e}",
        cva.correlations(),
        cva.spectrum_symmetry(),
        cva.eigen().b_orthonormality()
    );
}

/// With one variable on each side, the canonical correlation is the absolute
/// Pearson correlation.
///
/// The reference is the definition, written out here from the centred sums,
/// not anything this crate computes. `p = q = 1` leaves the canonical problem
/// with a single pair and no rotation to choose, so the two must agree.
#[test]
fn a_single_pair_reduces_to_the_pearson_correlation() {
    let mut rng = Lcg::new(0x9EA250);
    let mut worst = 0.0_f64;
    for case in 0..50 {
        let samples = 40;
        let mut x = vec![0.0; samples];
        let mut y = vec![0.0; samples];
        for t in 0..samples {
            let common = rng.next();
            x[t] = 3.0 + 10.0 * common + rng.next();
            // The sign alternates by case, so both signs of the underlying
            // correlation are exercised and the absolute value in the claim is
            // doing work.
            y[t] = -700.0 + if case % 2 == 0 { 1.0 } else { -1.0 } * common + 0.4 * rng.next();
        }

        let mean = |v: &[f64]| v.iter().sum::<f64>() / v.len() as f64;
        let (mx, my) = (mean(&x), mean(&y));
        let mut sxy = 0.0;
        let mut sxx = 0.0;
        let mut syy = 0.0;
        for t in 0..samples {
            sxy += (x[t] - mx) * (y[t] - my);
            sxx += (x[t] - mx) * (x[t] - mx);
            syy += (y[t] - my) * (y[t] - my);
        }
        let pearson = sxy / (sxx * syy).sqrt();

        let cva = Cva::fit(&x, &y, samples, 1, 1).expect("neither column is constant");
        assert_eq!(cva.pairs(), 1);
        worst = worst.max((cva.correlations()[0] - pearson.abs()).abs());
    }
    println!("p = q = 1 against the Pearson definition: worst error {worst:.3e}");
    assert!(worst < 1e-14, "{worst:.3e}");
}

/// When one block is an exact linear function of the other, every canonical
/// correlation is one.
///
/// `Y = X M` with `M` non-singular means the two blocks span the same space,
/// so every direction in `X` is perfectly predicted by some direction in `Y`
/// and all `min(p, q)` correlations are one.
///
/// **How exactly.** A canonical correlation is a cosine, so it cannot exceed
/// one, and [`Cva::correlations`] clamps at one for that reason. The
/// unclamped value is what this test reports, and it is one to within a few
/// ulp rather than bit-exactly: the route to it runs through a Cholesky
/// factorisation of `Sxx` and of `Syy = M' Sxx M`, and those two round
/// differently even though the ratio they form is exactly one in exact
/// arithmetic. The size of the gap is measured and printed here rather than
/// hidden behind the clamp, which is why [`Cva::raw_correlations`] exists.
#[test]
fn an_exactly_linear_pair_has_canonical_correlations_of_one() {
    let mut rng = Lcg::new(0x11EA12);
    let samples = 60;
    let mut worst_gap = 0.0_f64;
    let mut worst_ulps = 0_i64;
    for (p, q) in [(1_usize, 1_usize), (2, 2), (3, 3), (2, 4), (5, 3)] {
        let mut x = vec![0.0; samples * p];
        for value in &mut x {
            *value = 100.0 + 20.0 * rng.next();
        }
        // `Y = X M` with `M` a `p`-by-`q` matrix. The rank of the pair is
        // `min(p, q)`, so that many correlations must be one whichever way the
        // rectangle points; the extra columns of the wider block are exact
        // linear combinations of the narrower one, which is still an exact
        // linear function.
        let mut m = vec![0.0; p * q];
        for value in &mut m {
            *value = rng.next();
        }
        // Push the diagonal off zero so the map is not accidentally singular.
        for i in 0..p.min(q) {
            m[i * q + i] += 1.5;
        }
        let mut y = vec![0.0; samples * q];
        for t in 0..samples {
            for j in 0..q {
                let mut total = 0.0;
                for i in 0..p {
                    total += x[t * p + i] * m[i * q + j];
                }
                y[t * q + j] = total;
            }
        }

        let Some(cva) = Cva::fit(&x, &y, samples, p, q) else {
            // A rectangular `M` makes the wider block's cross-product matrix
            // singular, which `cholesky` reports rather than repairing. That
            // is the documented behaviour, so it is a pass for this shape and
            // the square shapes carry the claim.
            println!("p={p} q={q}: the wider block is rank deficient, as expected");
            continue;
        };
        for (j, &rho) in cva.raw_correlations().iter().enumerate() {
            let gap = (rho - 1.0).abs();
            worst_gap = worst_gap.max(gap);
            let ulps = (rho.to_bits() as i64 - 1.0_f64.to_bits() as i64).abs();
            worst_ulps = worst_ulps.max(ulps);
            assert!(
                gap < 1e-13,
                "p={p} q={q}: canonical correlation {j} is {rho:.17e}, not one"
            );
        }
        // Clamped, every one of them is at most exactly one.
        for &rho in cva.correlations() {
            assert!(
                rho <= 1.0,
                "a correlation cannot exceed one; got {rho:.17e}"
            );
        }
    }
    println!(
        "exactly linear pairs: worst |rho - 1| = {worst_gap:.3e} ({worst_ulps} ulp \
         from one), before the clamp to [0, 1]"
    );
}

/// Canonical correlations are invariant under a non-singular linear
/// transformation of either block.
///
/// The property that distinguishes CVA from PCA, and the reason
/// [`Cva::fit`] does not standardise its columns. Asserted rather than
/// asserted-in-prose because a fit that quietly standardised would still pass
/// every other test in this file.
#[test]
fn canonical_correlations_ignore_a_change_of_basis() {
    let mut rng = Lcg::new(0xBA515);
    let samples = 50;
    let (p, q) = (3_usize, 4_usize);
    let mut x = vec![0.0; samples * p];
    let mut y = vec![0.0; samples * q];
    for t in 0..samples {
        let a = rng.next();
        let b = rng.next();
        for i in 0..p {
            x[t * p + i] = a * (i as f64 + 1.0) + rng.next();
        }
        for j in 0..q {
            y[t * q + j] = 0.5 * a - 0.3 * b + rng.next();
        }
    }
    let plain = Cva::fit(&x, &y, samples, p, q).expect("full rank");

    // Rescale each `X` column by wildly different factors, and mix the `Y`
    // columns with a non-singular triangular map. Both are changes of basis.
    let mut x2 = x.clone();
    for t in 0..samples {
        for i in 0..p {
            x2[t * p + i] *= 10f64.powi(3 * i as i32 - 3);
        }
    }
    let mut y2 = vec![0.0; samples * q];
    for t in 0..samples {
        for j in 0..q {
            let mut total = 2.0 * y[t * q + j];
            for k in 0..j {
                total += 0.7 * y[t * q + k];
            }
            y2[t * q + j] = total;
        }
    }
    let transformed = Cva::fit(&x2, &y2, samples, p, q).expect("full rank");

    let mut worst = 0.0_f64;
    for (a, b) in plain.correlations().iter().zip(transformed.correlations()) {
        worst = worst.max((a - b).abs());
    }
    println!(
        "change of basis: correlations {:?} against {:?}, worst difference \
         {worst:.3e}",
        plain.correlations(),
        transformed.correlations()
    );
    assert!(worst < 1e-13, "{worst:.3e}");
}

/// The canonical variates really are correlated at the reported rate, and
/// uncorrelated across pairs.
///
/// The defining property, checked on the training data itself: `corr(u_j, v_j)`
/// is `rho_j`, and `corr(u_j, u_k)` and `corr(u_j, v_k)` are zero for `j != k`.
/// Computed here from the variates by the Pearson definition written out, so
/// the check does not reuse the machinery that produced them.
#[test]
fn the_canonical_variates_have_the_correlations_they_claim() {
    let mut rng = Lcg::new(0xA21A7E5);
    let samples = 200;
    let (p, q) = (3_usize, 3_usize);
    let mut x = vec![0.0; samples * p];
    let mut y = vec![0.0; samples * q];
    for t in 0..samples {
        let a = rng.next();
        let b = rng.next();
        x[t * p] = a + 0.1 * rng.next();
        x[t * p + 1] = b + 0.1 * rng.next();
        x[t * p + 2] = rng.next();
        y[t * q] = 0.9 * a + 0.2 * rng.next();
        y[t * q + 1] = 0.4 * b + 0.6 * rng.next();
        y[t * q + 2] = rng.next();
    }
    let cva = Cva::fit(&x, &y, samples, p, q).expect("full rank");

    let mut u = vec![vec![0.0; samples]; cva.pairs()];
    let mut v = vec![vec![0.0; samples]; cva.pairs()];
    for t in 0..samples {
        for (j, (a, b)) in cva
            .variates(&x[t * p..(t + 1) * p], &y[t * q..(t + 1) * q])
            .into_iter()
            .enumerate()
        {
            u[j][t] = a;
            v[j][t] = b;
        }
    }
    let correlate = |a: &[f64], b: &[f64]| -> f64 {
        let mean = |s: &[f64]| s.iter().sum::<f64>() / s.len() as f64;
        let (ma, mb) = (mean(a), mean(b));
        let mut sab = 0.0;
        let mut saa = 0.0;
        let mut sbb = 0.0;
        for t in 0..a.len() {
            sab += (a[t] - ma) * (b[t] - mb);
            saa += (a[t] - ma) * (a[t] - ma);
            sbb += (b[t] - mb) * (b[t] - mb);
        }
        sab / (saa * sbb).sqrt()
    };

    let mut worst_pair = 0.0_f64;
    let mut worst_cross = 0.0_f64;
    for j in 0..cva.pairs() {
        worst_pair = worst_pair.max((correlate(&u[j], &v[j]).abs() - cva.correlations()[j]).abs());
        for k in 0..cva.pairs() {
            if j != k {
                worst_cross = worst_cross.max(correlate(&u[j], &u[k]).abs());
                worst_cross = worst_cross.max(correlate(&u[j], &v[k]).abs());
                worst_cross = worst_cross.max(correlate(&v[j], &v[k]).abs());
            }
        }
    }
    // Every variate has mean zero over the training data, which is what the
    // centring in `variates` is for. A correlation is invariant under a
    // constant offset, so the two checks above hold whether or not the
    // observation is centred; a mutation run showed exactly that, and this is
    // the assertion that notices. The `X` block here is built around zero and
    // the `Y` block around zero too, so the check is made deliberately sharp by
    // shifting one of them far away first.
    let mut worst_mean = 0.0_f64;
    for j in 0..cva.pairs() {
        let scale = u[j].iter().map(|z| z.abs()).fold(1.0_f64, f64::max);
        worst_mean = worst_mean.max((u[j].iter().sum::<f64>() / samples as f64).abs() / scale);
        let scale = v[j].iter().map(|z| z.abs()).fold(1.0_f64, f64::max);
        worst_mean = worst_mean.max((v[j].iter().sum::<f64>() / samples as f64).abs() / scale);
    }
    let mut shifted = x.clone();
    for t in 0..samples {
        for i in 0..p {
            shifted[t * p + i] += 1000.0;
        }
    }
    let offset = Cva::fit(&shifted, &y, samples, p, q).expect("a shift changes no correlation");
    let mut worst_shift = 0.0_f64;
    for t in 0..samples {
        for (a, _) in offset.variates(&shifted[t * p..(t + 1) * p], &y[t * q..(t + 1) * q]) {
            worst_shift = worst_shift.max(a.abs());
        }
    }
    println!(
        "canonical variates: correlations {:?}, worst |corr(u_j, v_j)| - rho_j = \
         {worst_pair:.3e}, worst off-pair correlation {worst_cross:.3e}, \
         worst scaled variate mean {worst_mean:.3e}, largest variate after a \
         1000-unit shift of X {worst_shift:.3e}",
        cva.correlations()
    );
    assert!(worst_pair < 1e-13, "{worst_pair:.3e}");
    assert!(worst_cross < 1e-13, "{worst_cross:.3e}");
    assert!(worst_mean < 1e-13, "{worst_mean:.3e}");
    // A variate that is not centred on the training mean carries the 1000-unit
    // offset straight through, so it is orders larger than the deviations it is
    // supposed to describe.
    assert!(
        worst_shift < 100.0,
        "shifting X by 1000 units produced variates up to {worst_shift:.3e}, so \
         the observation is not being centred on the training mean"
    );
}

/// Exactly uncorrelated blocks give canonical correlations of exactly zero,
/// never a small negative number.
///
/// Two blocks built from disjoint Hadamard rows have `Sxy = 0` exactly, so
/// every canonical correlation is zero. That is the case the lower half of the
/// clamp in [`Cva::correlations`] exists for: the `+rho, -rho` pairing puts the
/// smallest reported correlation next to the largest unreported *negative*
/// eigenvalue, and at exactly zero an ulp of rounding decides which side of the
/// pair the sort puts first. A correlation is a magnitude and cannot be
/// negative, so the clamp is the right answer and this pins it.
#[test]
fn exactly_uncorrelated_blocks_give_exactly_zero() {
    let h = hadamard16();
    let samples = 16;
    let mut x = vec![0.0; samples * 2];
    let mut y = vec![0.0; samples * 2];
    for t in 0..samples {
        x[t * 2] = h[1][t];
        x[t * 2 + 1] = h[2][t];
        // Different rows entirely, so every `X`-`Y` inner product is exactly
        // zero and there is no linear relationship of any size to find.
        y[t * 2] = h[4][t];
        y[t * 2 + 1] = h[8][t];
    }
    let cva = Cva::fit(&x, &y, samples, 2, 2).expect("both blocks are full rank");
    println!(
        "orthogonal blocks: correlations {:?}, unclamped {:?}, spectrum {:?}",
        cva.correlations(),
        cva.raw_correlations(),
        cva.eigen().values()
    );
    for &rho in cva.correlations() {
        assert!(
            (0.0..1e-15).contains(&rho),
            "orthogonal blocks should give zero, got {rho:.17e}"
        );
        // Not merely non-negative: a correlation reported as -0.0 would satisfy
        // `>= 0.0` and print as a negative number.
        assert!(
            rho.is_sign_positive(),
            "a canonical correlation of {rho:.17e} carries a negative sign"
        );
    }
}

/// A constant column makes the cross-product matrix singular, and that is
/// reported rather than repaired.
#[test]
fn a_rank_deficient_block_returns_none() {
    let samples = 10;
    let mut x = vec![0.0; samples * 2];
    let mut y = vec![0.0; samples];
    for t in 0..samples {
        x[t * 2] = t as f64;
        // A constant second column: zero variance, so `Sxx` is singular.
        x[t * 2 + 1] = 7.0;
        y[t] = 2.0 * t as f64 + 1.0;
    }
    assert!(Cva::fit(&x, &y, samples, 2, 1).is_none());
    // A duplicated column is the same failure by a different route.
    let mut duplicated = vec![0.0; samples * 2];
    for t in 0..samples {
        duplicated[t * 2] = t as f64;
        duplicated[t * 2 + 1] = t as f64;
    }
    assert!(Cva::fit(&duplicated, &y, samples, 2, 1).is_none());
}

// ---------------------------------------------------------------------------
// The past/future construction
// ---------------------------------------------------------------------------

/// The past/future layout, checked index by index.
///
/// Off-by-one here is the classic dynamic-monitoring bug: it produces a model
/// that predicts the present from the present, which reports enormous
/// canonical correlations and detects nothing. It is invisible in any fitted
/// statistic, so the layout is checked directly against a matrix whose entries
/// encode their own coordinates.
#[test]
fn past_and_future_windows_line_up_where_they_should() {
    let samples = 9;
    let variables = 2;
    let mut data = vec![0.0; samples * variables];
    for t in 0..samples {
        for v in 0..variables {
            // `t.v`, so an entry names the observation and channel it came
            // from and a transposition is legible in the failure message.
            data[t * variables + v] = t as f64 + v as f64 / 10.0;
        }
    }

    let (past, future, rows) = past_future(&data, samples, variables, 3, 2);
    assert_eq!(rows, 9 + 1 - 3 - 2);
    let past_width = variables * 3;
    let future_width = variables * 2;
    for r in 0..rows {
        let present = 3 + r;
        for l in 0..3 {
            for v in 0..variables {
                assert_eq!(
                    past[r * past_width + l * variables + v],
                    (present - 1 - l) as f64 + v as f64 / 10.0,
                    "past row {r}, block {l}, variable {v}"
                );
            }
        }
        for m in 0..2 {
            for v in 0..variables {
                assert_eq!(
                    future[r * future_width + m * variables + v],
                    (present + m) as f64 + v as f64 / 10.0,
                    "future row {r}, block {m}, variable {v}"
                );
            }
        }
    }

    // A record exactly long enough for one window gives one row, and one
    // shorter gives none rather than a panic.
    let (_, _, one) = past_future(&data[..5 * variables], 5, variables, 3, 2);
    assert_eq!(one, 1);
    let (past, future, none) = past_future(&data[..4 * variables], 4, variables, 3, 2);
    assert_eq!(none, 0);
    assert!(past.is_empty() && future.is_empty());
}

/// On a first-order autoregressive series the leading canonical correlation
/// between past and future is the process's own lag-one correlation.
///
/// For `z_t = phi z_{t-1} + e_t` with independent `e`, the whole state is one
/// number, so a single canonical pair carries everything and its correlation
/// is `corr(z_{t-1}, z_t) = phi`. Constructed here from an exactly
/// reproducible generator with `phi = 0.8` and a long record, and compared
/// against `phi` itself, which is a property of the process and not of this
/// code.
///
/// The agreement is statistical, not numerical: a finite record estimates
/// `phi` with a standard error of about `sqrt((1 - phi^2) / n)`, which is
/// 0.008 here, and the autoregressive estimator is biased low by about
/// `2 phi / n`. So the bound is 0.05, three times the sampling error, and it
/// is a check that the construction points the right way rather than a
/// precision claim.
#[test]
fn a_first_order_process_shows_its_own_lag_one_correlation() {
    let mut rng = Lcg::new(0xA210);
    let phi = 0.8;
    let samples = 8000;
    let mut data = vec![0.0; samples];
    let mut z = 0.0;
    for t in 0..samples {
        // The generator is uniform, not normal; canonical correlation is a
        // second-moment quantity and does not care.
        z = phi * z + rng.next();
        data[t] = z;
    }

    let (past, future, rows) = past_future(&data, samples, 1, 1, 1);
    let cva = Cva::fit(&past, &future, rows, 1, 1).expect("neither block is constant");
    let estimate = cva.correlations()[0];
    println!(
        "AR(1) with phi = {phi}: leading canonical correlation over {rows} rows \
         is {estimate:.4}, against a sampling error of {:.4}",
        ((1.0 - phi * phi) / rows as f64).sqrt()
    );
    assert!(
        (estimate - phi).abs() < 0.05,
        "expected about {phi}, got {estimate}"
    );

    // With two lags of past against two leads of future the process is still
    // first order, so the *second* canonical correlation has to be far smaller
    // than the first: there is only one state to carry.
    let (past, future, rows) = past_future(&data, samples, 1, 2, 2);
    let cva = Cva::fit(&past, &future, rows, 2, 2).expect("full rank");
    let (first, second) = (cva.correlations()[0], cva.correlations()[1]);
    println!("AR(1) with two lags: canonical correlations {first:.4} and {second:.4}");
    assert!(
        second < 0.25 * first,
        "a first-order process should show one dominant canonical pair; got \
         {first:.4} and {second:.4}"
    );
}
