//! A symmetric eigensolver: the cyclic Jacobi method.
//!
//! Everything Tier 6 needs sits on top of this. PCA is an eigendecomposition of
//! a correlation matrix, Hotelling's T-squared divides scores by eigenvalues,
//! and SPE sums the ones PCA threw away. Get this wrong and every detector
//! number downstream is wrong in a way no amount of testing further down will
//! localise.
//!
//! # Why Jacobi and not the good algorithm
//!
//! For a general symmetric eigenproblem the right answer is tridiagonalise by
//! Householder reflections, then QL or QR with implicit shifts: `O(n^3)` with a
//! small constant, versus Jacobi's `O(n^3)` per sweep. Jacobi is here anyway,
//! for three reasons that outrank speed in this project.
//!
//! It is **deterministic by construction**. The sweep order is fixed (row-major
//! over the strict upper triangle), there is no pivot search whose tie-breaking
//! could differ, no deflation criterion, no shift strategy, and no iteration
//! count that depends on how the arithmetic rounded. `CLAUDE.md` makes bit
//! identity across x86-64, aarch64 and wasm32 a hard invariant, and a shifted QR
//! implementation would need a great deal of care to keep it.
//!
//! It is **small enough to read**. The whole method is one rotation, applied in
//! a loop. A reviewer can check it against the closed form on this page. The QL
//! implicit-shift routine is not something a reviewer checks by reading.
//!
//! It is **accurate in the way that matters here**. Jacobi computes the small
//! eigenvalues of a positive definite matrix to high *relative* accuracy, which
//! QR does not (Demmel and Veselic 1992). Those small eigenvalues are exactly
//! what SPE is built from: the residual subspace is the discarded tail of the
//! spectrum, and the Jackson-Mudholkar limit is a function of their first three
//! power sums. An absolute-accuracy method would deliver them as noise.
//!
//! The matrices Tier 6 diagonalises are at most a few hundred square (52
//! measured variables, or 52 times the lag count for dynamic PCA), and it is fit
//! once per detector. The cost never becomes the issue.
//!
//! # The rotation
//!
//! One step annihilates a single off-diagonal entry `a_pq` by a rotation in the
//! `(p, q)` plane. With `J` the identity except for
//! `J[p][p] = J[q][q] = c`, `J[p][q] = s`, `J[q][p] = -s`, and `B = J' A J`,
//!
//! ```text
//! B[p][q] = c s (a_pp - a_qq) + (c^2 - s^2) a_pq
//! ```
//!
//! Setting that to zero and writing `t = s / c` gives
//!
//! ```text
//! t^2 + 2 tau t - 1 = 0,     tau = (a_qq - a_pp) / (2 a_pq)
//! ```
//!
//! whose two roots are the rotations through `theta` and `theta + 90` degrees.
//! Taking the **smaller** root keeps `|theta| <= 45` degrees, which is the
//! choice that makes the method stable: it keeps each rotation close to the
//! identity, so the diagonal entries move by `t * a_pq` rather than being
//! rebuilt from scratch, and the eigenvalues stay near the diagonal entries they
//! started on. The larger root also annihilates `a_pq`, so a sign slip here is
//! not caught by "did it converge"; it is caught by the eigenvalues being wrong.
//!
//! Golub, G. H. and Van Loan, C. F. (2013), *Matrix Computations*, 4th ed.,
//! section 8.5. The formulation, including the `t = sign(tau) / (|tau| +
//! sqrt(1 + tau^2))` root and the `d_pp = a_pp - t a_pq` diagonal update, is
//! Algorithm 8.5.1 and 8.5.2.
//!
//! Demmel, J. and Veselic, K. (1992), "Jacobi's method is more accurate than
//! QR", *SIAM Journal on Matrix Analysis and Applications* 13(4), 1204-1245.

use alloc::vec;
use alloc::vec::Vec;
use core::cmp::Ordering;

use crate::special::{not_positive, sqrt};

/// The sweep budget.
///
/// Cyclic Jacobi converges quadratically once the off-diagonal norm is small,
/// and in practice needs six to ten sweeps for any matrix this project will
/// meet, at any size. Sixty is not a tuned number: it is far enough above the
/// observed requirement that hitting it means the method is not converging at
/// all, which [`SymmetricEigen::converged`] then reports rather than hiding.
const MAX_SWEEPS: usize = 60;

/// The eigenvalues and eigenvectors of a real symmetric matrix.
///
/// Eigenvalues are in **descending** order, which is the order PCA wants, and
/// the eigenvectors are permuted to match.
#[derive(Clone, Debug, PartialEq)]
pub struct SymmetricEigen {
    n: usize,
    values: Vec<f64>,
    /// `n` by `n`, row-major. Column `j` is the unit eigenvector for
    /// `values[j]`, so `vectors[i * n + j]` is component `i` of eigenvector `j`.
    vectors: Vec<f64>,
    sweeps: usize,
    off_diagonal_norm: f64,
    frobenius_norm: f64,
}

impl SymmetricEigen {
    /// The order of the matrix.
    #[must_use]
    pub const fn n(&self) -> usize {
        self.n
    }

    /// The eigenvalues, descending.
    #[must_use]
    pub fn values(&self) -> &[f64] {
        &self.values
    }

    /// The eigenvectors as columns of an `n` by `n` row-major matrix.
    ///
    /// `vectors()[i * n + j]` is component `i` of the eigenvector belonging to
    /// `values()[j]`.
    #[must_use]
    pub fn vectors(&self) -> &[f64] {
        &self.vectors
    }

    /// Component `i` of eigenvector `j`.
    ///
    /// # Panics
    ///
    /// If either index is out of range.
    #[must_use]
    pub fn component(&self, i: usize, j: usize) -> f64 {
        assert!(
            i < self.n && j < self.n,
            "({i}, {j}) is outside a {}-by-{} decomposition",
            self.n,
            self.n
        );
        self.vectors[i * self.n + j]
    }

    /// Eigenvector `j`, copied out.
    ///
    /// # Panics
    ///
    /// If `j` is out of range.
    #[must_use]
    pub fn eigenvector(&self, j: usize) -> Vec<f64> {
        assert!(j < self.n, "no eigenvector {j} in a {}-vector set", self.n);
        (0..self.n).map(|i| self.vectors[i * self.n + j]).collect()
    }

    /// How many sweeps it took.
    #[must_use]
    pub const fn sweeps(&self) -> usize {
        self.sweeps
    }

    /// The Frobenius norm of the off-diagonal part that survived, `off(A)`.
    ///
    /// The residual of the whole method, reported rather than swallowed: this
    /// is the number that says how nearly diagonal the transformed matrix
    /// actually became. `CLAUDE.md` asks for numbers rather than verdicts, and
    /// [`converged`](Self::converged) is the verdict.
    #[must_use]
    pub const fn off_diagonal_norm(&self) -> f64 {
        self.off_diagonal_norm
    }

    /// The Frobenius norm of the input, the scale `off_diagonal_norm` is
    /// measured against.
    #[must_use]
    pub const fn frobenius_norm(&self) -> f64 {
        self.frobenius_norm
    }

    /// Whether the off-diagonal norm fell below `eps` times the matrix norm
    /// inside the sweep budget.
    #[must_use]
    pub fn converged(&self) -> bool {
        self.off_diagonal_norm <= f64::EPSILON * self.frobenius_norm
    }

    /// Rebuild `V D V'` from the decomposition.
    ///
    /// Exists for the caller who wants to see the residual, and for the test
    /// that checks it: `A = V D V'` is the whole claim an eigensolver makes,
    /// and it is checkable to machine precision without a second eigensolver to
    /// compare against.
    #[must_use]
    pub fn reconstruct(&self) -> Vec<f64> {
        let n = self.n;
        let mut out = vec![0.0; n * n];
        for i in 0..n {
            for j in i..n {
                let mut total = 0.0;
                for k in 0..n {
                    total += self.vectors[i * n + k] * self.values[k] * self.vectors[j * n + k];
                }
                out[i * n + j] = total;
                // Written into both triangles from one sum, so the result is
                // exactly symmetric. Computing the lower half separately would
                // let the two halves disagree in the last bit and would make a
                // symmetry check on the output meaningless.
                out[j * n + i] = total;
            }
        }
        out
    }
}

/// Diagonalise a real symmetric matrix by cyclic Jacobi.
///
/// `matrix` is `n` by `n`, row-major.
///
/// # Sign convention
///
/// An eigenvector is only defined up to sign. This fixes it: the component of
/// largest magnitude is made positive, and the first such component wins a tie.
/// Without a convention the sign falls out of the rotation order, which makes
/// PCA loadings impossible to compare between two fits of nearly the same data
/// and makes a hand-checked test depend on an implementation detail.
///
/// The convention is reproducible, which is what it is for, but it is not
/// *stable*: when two components have nearly equal magnitude, an ulp decides
/// which of them leads and therefore which sign the whole vector takes. That is
/// not special to this rule. Every sign convention has a degenerate set, and
/// "largest component positive" has the smallest one of the usual choices. Do
/// not compare loadings between two fits component by component without
/// allowing for it; compare `|<u, v>|`.
///
/// # Ordering
///
/// Descending by eigenvalue, by a **stable** sort, so exactly equal eigenvalues
/// keep the order the sweeps left them in. A repeated eigenvalue has no
/// distinguished eigenvector, but it should at least have a reproducible one.
///
/// # Panics
///
/// If `matrix.len()` is not `n * n`, or if the matrix is not symmetric bit for
/// bit. The strict symmetry check is deliberate. Reading only one triangle, as
/// LAPACK does, would silently answer a question the caller did not ask when
/// the two halves disagree, and every matrix this crate builds is assembled
/// symmetric on purpose.
#[must_use]
pub fn symmetric_eigen(matrix: &[f64], n: usize) -> SymmetricEigen {
    assert_eq!(
        matrix.len(),
        n * n,
        "a {n}-by-{n} matrix needs {} entries, not {}",
        n * n,
        matrix.len()
    );
    for i in 0..n {
        for j in (i + 1)..n {
            assert!(
                matrix[i * n + j].to_bits() == matrix[j * n + i].to_bits(),
                "not symmetric: a[{i}][{j}] = {:.17e} but a[{j}][{i}] = {:.17e}",
                matrix[i * n + j],
                matrix[j * n + i]
            );
        }
    }

    let mut a = matrix.to_vec();
    let mut v = vec![0.0; n * n];
    for i in 0..n {
        v[i * n + i] = 1.0;
    }

    let frobenius_norm = frobenius(&a);
    let mut off_diagonal_norm = off_diagonal(&a, n);
    let mut sweeps = 0;
    while sweeps < MAX_SWEEPS && off_diagonal_norm > f64::EPSILON * frobenius_norm {
        // The fixed cyclic order: row-major over the strict upper triangle. Not
        // "largest off-diagonal first", which is the classical Jacobi method and
        // converges in fewer rotations but needs a search whose tie-breaking is
        // one more thing to pin down, and whose search cost dominates the
        // rotation cost for anything but a tiny matrix.
        for p in 0..n {
            for q in (p + 1)..n {
                rotate(&mut a, &mut v, n, p, q);
            }
        }
        sweeps += 1;
        off_diagonal_norm = off_diagonal(&a, n);
    }

    // Descending, stably.
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&i, &j| {
        a[j * n + j]
            .partial_cmp(&a[i * n + i])
            // A NaN eigenvalue means a NaN went in. Comparing equal keeps the
            // sort a valid total order instead of producing a garbage
            // permutation on top of the garbage values.
            .unwrap_or(Ordering::Equal)
    });

    let values: Vec<f64> = order.iter().map(|&j| a[j * n + j]).collect();
    let mut vectors = vec![0.0; n * n];
    for (new_column, &old_column) in order.iter().enumerate() {
        let mut leader = 0;
        let mut largest = 0.0;
        for i in 0..n {
            let magnitude = libm::fabs(v[i * n + old_column]);
            if magnitude > largest {
                largest = magnitude;
                leader = i;
            }
        }
        // Negation flips one bit and rounds nothing, so the flipped vector is
        // exactly as orthonormal as the one the sweeps produced.
        let flip = n > 0 && v[leader * n + old_column] < 0.0;
        for i in 0..n {
            let value = v[i * n + old_column];
            vectors[i * n + new_column] = if flip { -value } else { value };
        }
    }

    SymmetricEigen {
        n,
        values,
        vectors,
        sweeps,
        off_diagonal_norm,
        frobenius_norm,
    }
}

/// One Jacobi rotation, annihilating `a[p][q]`.
///
/// Updates `a` in place on both triangles and accumulates the rotation into
/// `v`, so that `v` ends up holding the product of every rotation applied and
/// `a` ends up diagonal.
#[allow(
    clippy::suboptimal_flops,
    reason = "these expressions transcribe the closed forms in the module docs; \
              folding them into fused multiply-adds would change the rounding \
              for a speed gain this crate has no use for"
)]
fn rotate(a: &mut [f64], v: &mut [f64], n: usize, p: usize, q: usize) {
    let apq = a[p * n + q];
    // Already zero, so there is nothing to annihilate and the rotation would be
    // the identity. Phrased as a negated `>` rather than `== 0.0` so that a NaN
    // pivot also falls out here instead of producing a NaN rotation that then
    // spreads across the whole matrix.
    if not_positive(libm::fabs(apq)) {
        return;
    }
    let app = a[p * n + p];
    let aqq = a[q * n + q];

    let tau = 0.5 * (aqq - app) / apq;
    // The smaller root of `t^2 + 2 tau t - 1 = 0`, written so that neither
    // branch subtracts two nearly equal numbers. The naive `tau - sqrt(1 + tau
    // * tau)` loses every significant digit for large positive `tau`, which is
    // the common case: `tau` is large exactly when the pivot is already small,
    // which is exactly when the sweeps are close to done.
    let t = if libm::fabs(tau) > 1e150 {
        // `1 + tau * tau` overflows above about 1.3e154. There `t` is `1 / (2
        // tau)` to full precision, which is below 3e-151: the rotation is the
        // identity to every bit that exists, and `a[p][q]`, being smaller than
        // the diagonal by that factor, is below any convergence threshold.
        0.5 / tau
    } else if tau >= 0.0 {
        1.0 / (tau + sqrt(1.0 + tau * tau))
    } else {
        -1.0 / (-tau + sqrt(1.0 + tau * tau))
    };
    let c = 1.0 / sqrt(1.0 + t * t);
    let s = t * c;

    // `d_pp = a_pp - t a_pq` and `d_qq = a_qq + t a_pq` are exact consequences
    // of the quadratic `t` solves, and they are what the expanded forms
    // `c^2 a_pp - 2 s c a_pq + s^2 a_qq` reduce to. Using the reduced form is
    // not a shortcut: the expanded one cancels `c^2 a_pp` against `s^2 a_qq`
    // when the two diagonal entries are close, which is the case a nearly
    // degenerate pair of eigenvalues produces.
    let shift = t * apq;
    a[p * n + p] = app - shift;
    a[q * n + q] = aqq + shift;
    // Set to exactly zero rather than computing the rounded value. The rotation
    // was chosen to annihilate this entry, so zero is the answer; letting
    // rounding leave a few ulp there would stop the off-diagonal norm from
    // decreasing monotonically and could keep the sweeps alive forever.
    a[p * n + q] = 0.0;
    a[q * n + p] = 0.0;

    for k in 0..n {
        if k == p || k == q {
            continue;
        }
        let akp = a[k * n + p];
        let akq = a[k * n + q];
        let rotated_p = c * akp - s * akq;
        let rotated_q = s * akp + c * akq;
        // Written into both triangles from the same value, so `a` stays exactly
        // symmetric through every rotation.
        a[k * n + p] = rotated_p;
        a[p * n + k] = rotated_p;
        a[k * n + q] = rotated_q;
        a[q * n + k] = rotated_q;
    }

    for k in 0..n {
        let vkp = v[k * n + p];
        let vkq = v[k * n + q];
        v[k * n + p] = c * vkp - s * vkq;
        v[k * n + q] = s * vkp + c * vkq;
    }
}

/// The Frobenius norm of a matrix, `sqrt(sum a_ij^2)`.
fn frobenius(a: &[f64]) -> f64 {
    sqrt(a.iter().map(|x| x * x).sum::<f64>())
}

/// `off(A)`: the Frobenius norm of the off-diagonal part.
///
/// The quantity cyclic Jacobi drives to zero, and the one Golub and Van Loan's
/// convergence proof is stated in terms of.
fn off_diagonal(a: &[f64], n: usize) -> f64 {
    let mut total = 0.0;
    for i in 0..n {
        for j in 0..n {
            if i != j {
                let value = a[i * n + j];
                total += value * value;
            }
        }
    }
    sqrt(total)
}
