//! Canonical variate analysis, and the generalised symmetric eigenproblem it
//! is built on.
//!
//! `PLAN.org` names CVA beside PCA and DPCA as a Tier 6 detector. It does not
//! fit on the machinery those two use. PCA diagonalises one correlation
//! matrix, which [`crate::eigen`] does directly; CVA asks a two-sided question,
//! "which direction in one block of variables is most correlated with which
//! direction in another", and that is a **generalised** symmetric eigenproblem
//!
//! ```text
//! A v = lambda B v,     A symmetric, B symmetric positive definite
//! ```
//!
//! which cyclic Jacobi does not solve. This module reduces it to one that
//! Jacobi does solve, and then builds canonical correlations on top.
//!
//! Hotelling, H. (1936), "Relations between two sets of variates",
//! *Biometrika* 28(3-4), 321-377. The generalised eigenproblem below is his
//! formulation of the problem.
//!
//! Golub, G. H. and Van Loan, C. F. (2013), *Matrix Computations*, 4th ed.,
//! section 8.7.1, for the Cholesky reduction and the conditioning caveat that
//! comes with it.
//!
//! Larimore, W. E. (1990), "Canonical variate analysis in identification,
//! filtering, and adaptive control", *Proceedings of the 29th IEEE Conference
//! on Decision and Control*, 596-604, for the past/future construction
//! [`past_future`] implements.
//!
//! Russell, E. L., Chiang, L. H. and Braatz, R. D. (2000), "Fault detection in
//! industrial processes using canonical variate analysis and dynamic principal
//! component analysis", *Chemometrics and Intelligent Laboratory Systems*
//! 51(1), 81-93, for CVA as a TEP fault detector.
//!
//! # The reduction, and why this one
//!
//! `B` is positive definite, so it has a Cholesky factorisation `B = L L'`
//! with `L` lower triangular and positive on the diagonal. Substituting and
//! writing `y = L' v`,
//!
//! ```text
//! A v = lambda B v
//! A L^-T y = lambda L L' L^-T y = lambda L y
//! (L^-1 A L^-T) y = lambda y
//! ```
//!
//! so the eigenvalues of the pencil are exactly the eigenvalues of
//! `C = L^-1 A L^-T`, which is symmetric because `A` is, and the eigenvectors
//! come back as `v = L^-T y`.
//!
//! The obvious alternative, forming `B^-1 A` and diagonalising that, is wrong
//! in a way that does not announce itself. `B^-1 A` is **not symmetric**, so a
//! symmetric eigensolver applied to it is being lied to, and a general one
//! gives up the guarantee that the eigenvalues are real and the eigenvectors
//! orthogonal. The eigenvalues of a real pencil with `B` definite *are* real;
//! an unsymmetric solver would return them with small imaginary parts and a
//! reader would have to decide which of those are rounding. The Cholesky route
//! keeps symmetry the whole way, so the answer is real by construction.
//!
//! # The numerical hazards, and what is done about each
//!
//! **`B` near singular.** The reduction multiplies by `L^-1` twice, so an
//! error in `A` of size `eps ||A||` becomes an error in `C` of size
//! `eps ||A|| / lambda_min(B)`. The method is accurate to about
//! `eps * cond(B)`, and there is no arrangement of it that is not: the problem
//! itself is that ill conditioned. This is not hidden. [`cholesky`] returns
//! `None` rather than a factorisation when `B` is not positive definite, and
//! [`GeneralizedSymmetricEigen::smallest_pivot`] reports the smallest diagonal
//! entry of `L`, whose square is the scale the reduction divides by. A caller
//! whose `B` is a covariance matrix of nearly collinear columns will see it
//! there.
//!
//! **Getting the reduction backwards.** `L^-1 A L^-T` and `L^-T A L^-1` are
//! both symmetric and both look right; only the first pairs with
//! `v = L^-T y`. Mixing them gives the spectrum of a different pencil, with no
//! symptom other than wrong numbers. The `B = I` known-answer test cannot see
//! the difference, because `L` is then the identity, which is exactly why the
//! test battery also carries a non-diagonal `B` with an independent closed
//! form.
//!
//! **The computed `C` is only symmetric to rounding.** Mathematically
//! `L^-1 A L^-T` is symmetric, but it is computed as two triangular solves and
//! the two triangles round differently, by an ulp or so.
//! [`crate::symmetric_eigen`] asserts bit-exact symmetry on purpose, so the
//! two halves are averaged before it is called. `(x + x) / 2` is exact in
//! binary floating point, so a genuinely symmetric `C`, which is what the
//! `B = I` case produces, survives the averaging unchanged.
//!
//! **Normalisation that is not.** The back-transformed vectors satisfy
//! `V' B V = I` in exact arithmetic, because `y' y = 1` and `v' B v = y' y`.
//! They are **not** rescaled afterwards to force it, and that is deliberate:
//! the deviation is the reduction's own accuracy, reported by
//! [`GeneralizedSymmetricEigen::b_orthonormality`], and dividing it away would
//! hide the one number that says how much `cond(B)` cost. It also keeps
//! `B = I` returning the standard problem's answer bit for bit, since
//! rescaling would divide by a number that is one only to within an ulp.
//!
//! # What CVA does with it
//!
//! Given two blocks of variables observed together, `X` with `p` columns and
//! `Y` with `q`, Hotelling's canonical correlations are the stationary values
//! of `corr(X a, Y b)`. Writing `Sxx`, `Syy` and `Sxy` for the centred
//! cross-product matrices, they solve
//!
//! ```text
//! [  0   Sxy ] [a]          [ Sxx   0  ] [a]
//! [ Syx   0  ] [b]  = rho * [  0   Syy ] [b]
//! ```
//!
//! which is a generalised symmetric eigenproblem of order `p + q` with a
//! block-diagonal, positive definite `B`. Eliminating `b` recovers the
//! textbook `Sxy Syy^-1 Syx a = rho^2 Sxx a`, so the eigenvalues of this
//! pencil are the canonical correlations themselves rather than their squares,
//! and they come in `+rho, -rho` pairs. That pairing is a strong free check:
//! [`Cva::spectrum_symmetry`] measures how far it is from exact, and nothing
//! in the computation was told to enforce it.
//!
//! This form is used rather than building `Sxy Syy^-1 Syx` because that
//! product needs an explicit inverse, squares the condition number, and
//! returns `rho^2`, from which a small `rho` is recovered with half the digits
//! it went in with.

use alloc::vec;
use alloc::vec::Vec;

use crate::eigen::symmetric_eigen;
use crate::special::{not_positive, sqrt};

/// The Cholesky factorisation `B = L L'` of a symmetric positive definite
/// matrix.
///
/// `matrix` is `n` by `n`, row-major. Returns the **lower** triangular `L`,
/// also `n` by `n` row-major, with exact zeros above the diagonal, or `None`
/// if the matrix is not positive definite.
///
/// The right-looking, column-oriented form: for `j <= i`,
///
/// ```text
/// L[j][j] = sqrt( B[j][j] - sum_{k<j} L[j][k]^2 )
/// L[i][j] = ( B[i][j] - sum_{k<j} L[i][k] L[j][k] ) / L[j][j]
/// ```
///
/// # Why `None` rather than a panic or a repair
///
/// A non-positive pivot means the matrix is not positive definite, which for a
/// covariance matrix means the columns are linearly dependent: a constant
/// column, a duplicated variable, or more variables than observations. Adding
/// a small multiple of the identity until the factorisation succeeds is the
/// usual patch and it is a tuned number wearing a numerical disguise, which
/// this project does not allow. Returning `None` makes the caller decide, and
/// makes the failure visible in a report rather than absorbed into a spectrum.
///
/// The test is `pivot > 0`, phrased through this crate's `not_positive` so
/// that a `NaN` pivot, which is what a `NaN` in the input produces, also fails
/// rather than propagating. Writing it as `pivot <= 0.0` would be *false* for
/// `NaN` and would turn the guard into a passthrough.
///
/// # Panics
///
/// If `matrix.len()` is not `n * n`, or if the matrix is not symmetric bit for
/// bit. Same rule and same reason as [`crate::symmetric_eigen`]: reading one
/// triangle would answer a question the caller did not ask.
#[must_use]
pub fn cholesky(matrix: &[f64], n: usize) -> Option<Vec<f64>> {
    assert_symmetric(matrix, n, "cholesky");
    let mut l = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..=i {
            let mut sum = matrix[i * n + j];
            for k in 0..j {
                sum -= l[i * n + k] * l[j * n + k];
            }
            if i == j {
                if not_positive(sum) {
                    return None;
                }
                l[i * n + i] = sqrt(sum);
            } else {
                l[i * n + j] = sum / l[j * n + j];
            }
        }
    }
    Some(l)
}

/// Solve `L X = M` for `X`, by forward substitution.
///
/// `l` is `n` by `n` lower triangular, `m` is `n` by `cols`, both row-major.
///
/// The identity case matters and is why this is written out rather than
/// specialised: with `L = I` every subtracted term is `0.0 * x`, so `sum` is
/// the right-hand side unchanged and the division is by exactly one. The
/// result is then `M` bit for bit, signed zeros included, which is what makes
/// [`generalized_symmetric_eigen`] return the standard problem's answer
/// exactly when `B` is the identity.
fn forward_substitute(l: &[f64], n: usize, m: &[f64], cols: usize) -> Vec<f64> {
    let mut x = vec![0.0; n * cols];
    for i in 0..n {
        for c in 0..cols {
            let mut sum = m[i * cols + c];
            for k in 0..i {
                sum -= l[i * n + k] * x[k * cols + c];
            }
            x[i * cols + c] = sum / l[i * n + i];
        }
    }
    x
}

/// Solve `L' v = y` for `v`, by back substitution.
///
/// `l` is the same lower triangular factor; the system solved is against its
/// transpose, which is upper triangular. This is the back-transform
/// `v = L^-T y` that turns an eigenvector of the reduced problem into one of
/// the pencil.
fn back_substitute_transpose(l: &[f64], n: usize, y: &[f64]) -> Vec<f64> {
    let mut v = y.to_vec();
    for i in (0..n).rev() {
        let mut sum = v[i];
        for k in (i + 1)..n {
            // `L'[i][k]` is `L[k][i]`.
            sum -= l[k * n + i] * v[k];
        }
        v[i] = sum / l[i * n + i];
    }
    v
}

/// Transpose an `n` by `n` row-major matrix.
fn transpose(m: &[f64], n: usize) -> Vec<f64> {
    let mut out = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..n {
            out[j * n + i] = m[i * n + j];
        }
    }
    out
}

/// Panic unless `matrix` is `n` by `n` and symmetric bit for bit.
fn assert_symmetric(matrix: &[f64], n: usize, who: &str) {
    assert_eq!(
        matrix.len(),
        n * n,
        "{who}: a {n}-by-{n} matrix needs {} entries, not {}",
        n * n,
        matrix.len()
    );
    for i in 0..n {
        for j in (i + 1)..n {
            assert!(
                matrix[i * n + j].to_bits() == matrix[j * n + i].to_bits(),
                "{who}: not symmetric: a[{i}][{j}] = {:.17e} but a[{j}][{i}] = {:.17e}",
                matrix[i * n + j],
                matrix[j * n + i]
            );
        }
    }
}

/// The solution of `A v = lambda B v`.
///
/// Eigenvalues descending, eigenvectors permuted to match, and the diagnostics
/// that say how much the reduction cost.
#[derive(Clone, Debug, PartialEq)]
pub struct GeneralizedSymmetricEigen {
    n: usize,
    values: Vec<f64>,
    /// `n` by `n`, row-major. Column `j` is the eigenvector for `values[j]`.
    vectors: Vec<f64>,
    sweeps: usize,
    off_diagonal_norm: f64,
    frobenius_norm: f64,
    smallest_pivot: f64,
    largest_pivot: f64,
    b_orthonormality: f64,
}

impl GeneralizedSymmetricEigen {
    /// The order of the pencil.
    #[must_use]
    pub const fn n(&self) -> usize {
        self.n
    }

    /// The generalised eigenvalues, descending.
    #[must_use]
    pub fn values(&self) -> &[f64] {
        &self.values
    }

    /// The eigenvectors as columns of an `n` by `n` row-major matrix.
    ///
    /// `vectors()[i * n + j]` is component `i` of the eigenvector belonging to
    /// `values()[j]`. They are `B`-orthonormal, `V' B V = I`, to within
    /// [`b_orthonormality`](Self::b_orthonormality).
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

    /// How many Jacobi sweeps the reduced problem needed.
    #[must_use]
    pub const fn sweeps(&self) -> usize {
        self.sweeps
    }

    /// `off(C)` for the reduced matrix: the eigensolver's own residual.
    #[must_use]
    pub const fn off_diagonal_norm(&self) -> f64 {
        self.off_diagonal_norm
    }

    /// `||C||_F`, the scale [`off_diagonal_norm`](Self::off_diagonal_norm) is
    /// measured against.
    #[must_use]
    pub const fn frobenius_norm(&self) -> f64 {
        self.frobenius_norm
    }

    /// Whether the reduced problem's off-diagonal norm fell below `eps` times
    /// its Frobenius norm inside the sweep budget.
    #[must_use]
    pub fn converged(&self) -> bool {
        self.off_diagonal_norm <= f64::EPSILON * self.frobenius_norm
    }

    /// The smallest diagonal entry of the Cholesky factor `L`.
    ///
    /// Its square is roughly the smallest eigenvalue of `B`, and the reduction
    /// divides by it twice. Read together with
    /// [`largest_pivot`](Self::largest_pivot): the ratio of their squares is a
    /// lower bound on `cond(B)`, and `eps` times that is about the accuracy
    /// available. This is the number to look at before believing a spectrum.
    #[must_use]
    pub const fn smallest_pivot(&self) -> f64 {
        self.smallest_pivot
    }

    /// The largest diagonal entry of the Cholesky factor.
    #[must_use]
    pub const fn largest_pivot(&self) -> f64 {
        self.largest_pivot
    }

    /// The largest deviation of `V' B V` from the identity.
    ///
    /// Measured, not corrected. See the module docs.
    #[must_use]
    pub const fn b_orthonormality(&self) -> f64 {
        self.b_orthonormality
    }
}

/// Solve `A v = lambda B v` for symmetric `A` and symmetric positive definite
/// `B`.
///
/// Both matrices are `n` by `n`, row-major. Returns `None` when `B` is not
/// positive definite, which is [`cholesky`]'s verdict and not a threshold
/// applied here.
///
/// # Sign convention and ordering
///
/// Descending by eigenvalue, and each eigenvector's largest-magnitude
/// component made positive, the first such component winning a tie. Applied
/// **after** the back-transform, because `L^-T` can move which component is
/// largest. Negation flips one bit and rounds nothing, so re-applying it costs
/// no accuracy, and it leaves the `B = I` case alone, where the standard
/// solver had already imposed the same rule.
///
/// The convention's degenerate set is the same as [`crate::symmetric_eigen`]'s
/// and is documented there: two components of nearly equal magnitude let an
/// ulp decide the sign of the whole vector. Compare `|<u, v>|`, not components.
///
/// # Panics
///
/// If either matrix is the wrong size or is not symmetric bit for bit.
#[must_use]
pub fn generalized_symmetric_eigen(
    a: &[f64],
    b: &[f64],
    n: usize,
) -> Option<GeneralizedSymmetricEigen> {
    assert_symmetric(a, n, "generalized_symmetric_eigen: A");
    assert_symmetric(b, n, "generalized_symmetric_eigen: B");

    let l = cholesky(b, n)?;
    let mut smallest_pivot = f64::INFINITY;
    let mut largest_pivot = 0.0_f64;
    for i in 0..n {
        smallest_pivot = smallest_pivot.min(l[i * n + i]);
        largest_pivot = largest_pivot.max(l[i * n + i]);
    }

    // `C = L^-1 A L^-T`, in two forward solves. `A` is symmetric, so
    // `(L^-1 A)' = A L^-T` and the second solve finishes the job.
    let left = forward_substitute(&l, n, a, n);
    let reduced = forward_substitute(&l, n, &transpose(&left, n), n);

    // Symmetrise exactly. `(x + x) / 2` is exact, so a `C` that really is
    // symmetric, which is what `B = I` gives, passes through untouched.
    //
    // Taking one triangle and mirroring it would be equally correct and is not
    // distinguishable by any test: the two halves differ by at most an ulp, and
    // the average of two values an ulp apart is one of them. A mutation run
    // confirmed it survives the whole battery. The average is kept because it
    // does not privilege a triangle, not because anything can tell.
    let mut c = vec![0.0; n * n];
    for i in 0..n {
        for j in i..n {
            let value = 0.5 * (reduced[i * n + j] + reduced[j * n + i]);
            c[i * n + j] = value;
            c[j * n + i] = value;
        }
    }

    let eigen = symmetric_eigen(&c, n);

    // Back-transform, then re-impose the sign convention.
    let mut vectors = vec![0.0; n * n];
    for j in 0..n {
        let v = back_substitute_transpose(&l, n, &eigen.eigenvector(j));
        let mut leader = 0;
        let mut largest = 0.0;
        for (i, value) in v.iter().enumerate() {
            let magnitude = libm::fabs(*value);
            if magnitude > largest {
                largest = magnitude;
                leader = i;
            }
        }
        let flip = n > 0 && v[leader] < 0.0;
        for (i, value) in v.iter().enumerate() {
            vectors[i * n + j] = if flip { -value } else { *value };
        }
    }

    // `V' B V - I`, the residual the module docs decline to divide away.
    let mut b_orthonormality = 0.0_f64;
    let mut bv = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..n {
            let mut total = 0.0;
            for k in 0..n {
                total += b[i * n + k] * vectors[k * n + j];
            }
            bv[i * n + j] = total;
        }
    }
    for i in 0..n {
        for j in 0..n {
            let mut total = 0.0;
            for k in 0..n {
                total += vectors[k * n + i] * bv[k * n + j];
            }
            let expected = if i == j { 1.0 } else { 0.0 };
            b_orthonormality = b_orthonormality.max(libm::fabs(total - expected));
        }
    }

    Some(GeneralizedSymmetricEigen {
        n,
        values: eigen.values().to_vec(),
        vectors,
        sweeps: eigen.sweeps(),
        off_diagonal_norm: eigen.off_diagonal_norm(),
        frobenius_norm: eigen.frobenius_norm(),
        smallest_pivot: if n == 0 { 0.0 } else { smallest_pivot },
        largest_pivot,
        b_orthonormality,
    })
}

/// A fitted canonical correlation analysis.
///
/// The canonical correlations, the weight vectors that produce them, and the
/// diagnostics that say whether to believe either.
#[derive(Clone, Debug, PartialEq)]
pub struct Cva {
    samples: usize,
    x_variables: usize,
    y_variables: usize,
    /// The training column means of the `X` block, so [`Cva::variates`] can
    /// centre a new observation the way the fit was centred.
    x_means: Vec<f64>,
    /// The same for the `Y` block.
    y_means: Vec<f64>,
    correlations: Vec<f64>,
    raw_correlations: Vec<f64>,
    /// `x_variables` by `pairs`, row-major. Column `j` is `a_j`.
    x_weights: Vec<f64>,
    /// `y_variables` by `pairs`, row-major. Column `j` is `b_j`.
    y_weights: Vec<f64>,
    eigen: GeneralizedSymmetricEigen,
    spectrum_symmetry: f64,
}

impl Cva {
    /// Fit to two blocks of variables observed together.
    ///
    /// `x` is `samples` by `x_variables` and `y` is `samples` by
    /// `y_variables`, both row-major, row `t` of each being the same
    /// observation. Returns `None` when either block's cross-product matrix is
    /// singular, which is [`cholesky`]'s verdict: a constant column, a
    /// duplicated variable, or fewer observations than variables.
    ///
    /// # Centring, and the scale that cancels
    ///
    /// Both blocks are centred on their column means, by the same compensated
    /// two-pass summation [`crate::pca::Pca::fit`] uses and for the same
    /// reason. They are **not** standardised, and the cross-product matrices
    /// are not divided by `samples - 1`. Neither omission changes the answer:
    /// a canonical correlation is invariant under any non-singular linear
    /// transformation of either block, so column scaling cancels, and dividing
    /// `A` and `B` by the same constant leaves `A v = lambda B v` unchanged.
    /// Leaving the division out keeps the arithmetic exact on integer data,
    /// which is what makes an exact known-answer test possible.
    ///
    /// That invariance is worth stating plainly because it is the difference
    /// between CVA and PCA: PCA on a covariance matrix is dominated by
    /// whichever variable has the largest units, which is why
    /// [`crate::pca`] always standardises. CVA cannot be, so it does not need
    /// to.
    ///
    /// # Panics
    ///
    /// If either matrix's length disagrees with its shape, if the two blocks
    /// have different sample counts, if either block has no columns, or if
    /// there are fewer than two samples.
    #[must_use]
    pub fn fit(
        x: &[f64],
        y: &[f64],
        samples: usize,
        x_variables: usize,
        y_variables: usize,
    ) -> Option<Self> {
        assert!(x_variables > 0 && y_variables > 0, "CVA needs both blocks");
        assert!(
            samples >= 2,
            "a cross-product matrix needs at least two samples, got {samples}"
        );
        assert_eq!(
            x.len(),
            samples * x_variables,
            "the X block is {} entries, not {}",
            x.len(),
            samples * x_variables
        );
        assert_eq!(
            y.len(),
            samples * y_variables,
            "the Y block is {} entries, not {}",
            y.len(),
            samples * y_variables
        );

        let (cx, x_means) = centre(x, samples, x_variables);
        let (cy, y_means) = centre(y, samples, y_variables);

        let n = x_variables + y_variables;
        // `A = [[0, Sxy], [Syx, 0]]`, `B = blkdiag(Sxx, Syy)`. Assembled
        // directly into the two `n`-by-`n` matrices, symmetric by
        // construction: every entry is written into both triangles from one
        // sum, so `assert_symmetric` cannot fail on a rounding difference
        // between two separately accumulated dot products.
        let mut a = vec![0.0; n * n];
        let mut b = vec![0.0; n * n];
        for i in 0..x_variables {
            for j in i..x_variables {
                let value = dot(&cx, samples, x_variables, i, &cx, x_variables, j);
                b[i * n + j] = value;
                b[j * n + i] = value;
            }
        }
        for i in 0..y_variables {
            for j in i..y_variables {
                let value = dot(&cy, samples, y_variables, i, &cy, y_variables, j);
                b[(x_variables + i) * n + x_variables + j] = value;
                b[(x_variables + j) * n + x_variables + i] = value;
            }
        }
        for i in 0..x_variables {
            for j in 0..y_variables {
                let value = dot(&cx, samples, x_variables, i, &cy, y_variables, j);
                a[i * n + x_variables + j] = value;
                a[(x_variables + j) * n + i] = value;
            }
        }

        let eigen = generalized_symmetric_eigen(&a, &b, n)?;

        // The spectrum must be symmetric about zero: if `(a, b)` solves the
        // pencil at `+rho` then `(a, -b)` solves it at `-rho`. Nothing above
        // enforces that, so how nearly it holds is a free check on the whole
        // reduction.
        let mut spectrum_symmetry = 0.0_f64;
        for j in 0..n {
            spectrum_symmetry =
                spectrum_symmetry.max(libm::fabs(eigen.values()[j] + eigen.values()[n - 1 - j]));
        }

        let pairs = x_variables.min(y_variables);
        let raw_correlations: Vec<f64> = eigen.values()[..pairs].to_vec();
        // Clamped, for the reason `pca.rs` clamps a correlation: the quantity
        // is a cosine and cannot leave `[-1, 1]`, so a value outside it is
        // rounding on an exactly collinear pair and nothing else. Both the
        // clamped and the unclamped values are kept, because the size of the
        // excess is the interesting part and a clamp that hides it would be a
        // tuned tolerance with no name.
        //
        // The upper half is reached and tested. The lower half is defensive and
        // no test reaches it: a mutation to `clamp(-1.0, 1.0)` survives the
        // whole battery, because a construction whose smallest canonical
        // correlation is a hair below zero has to have `Sxy` exactly zero, and
        // then the whole pencil is exactly zero and so is the answer. It stays
        // because `pairs` can equal `n / 2`, which puts the last reported value
        // on the boundary of the `+rho, -rho` pair where an ulp decides the
        // sign, and a correlation reported as `-1e-17` would be a lie about a
        // magnitude.
        let correlations: Vec<f64> = raw_correlations.iter().map(|r| r.clamp(0.0, 1.0)).collect();

        let mut x_weights = vec![0.0; x_variables * pairs];
        let mut y_weights = vec![0.0; y_variables * pairs];
        for j in 0..pairs {
            for i in 0..x_variables {
                x_weights[i * pairs + j] = eigen.component(i, j);
            }
            for i in 0..y_variables {
                y_weights[i * pairs + j] = eigen.component(x_variables + i, j);
            }
        }

        Some(Self {
            samples,
            x_variables,
            y_variables,
            x_means,
            y_means,
            correlations,
            raw_correlations,
            x_weights,
            y_weights,
            eigen,
            spectrum_symmetry,
        })
    }

    /// How many observations it was fitted to.
    #[must_use]
    pub const fn samples(&self) -> usize {
        self.samples
    }

    /// Columns in the `X` block.
    #[must_use]
    pub const fn x_variables(&self) -> usize {
        self.x_variables
    }

    /// Columns in the `Y` block.
    #[must_use]
    pub const fn y_variables(&self) -> usize {
        self.y_variables
    }

    /// How many canonical pairs there are: `min(p, q)`.
    #[must_use]
    pub fn pairs(&self) -> usize {
        self.correlations.len()
    }

    /// The canonical correlations, descending, each in `[0, 1]`.
    #[must_use]
    pub fn correlations(&self) -> &[f64] {
        &self.correlations
    }

    /// The same before the clamp to `[0, 1]`.
    ///
    /// Equal to [`correlations`](Self::correlations) except on data where one
    /// block is an exact linear function of the other, where the leading
    /// values are one to within rounding and can land an ulp above it.
    #[must_use]
    pub fn raw_correlations(&self) -> &[f64] {
        &self.raw_correlations
    }

    /// The `X` weights: `x_variables` by `pairs`, row-major, column `j` being
    /// the direction `a_j` in the `X` block.
    #[must_use]
    pub fn x_weights(&self) -> &[f64] {
        &self.x_weights
    }

    /// The `Y` weights: `y_variables` by `pairs`, row-major.
    #[must_use]
    pub fn y_weights(&self) -> &[f64] {
        &self.y_weights
    }

    /// The canonical variates of one observation: `(X a_j, Y b_j)` for every
    /// pair `j`.
    ///
    /// The observation is given as its raw `X` and `Y` rows; centring is
    /// applied here with the training means, so the variates of a training row
    /// have mean zero over the training set.
    ///
    /// # Panics
    ///
    /// If either row is the wrong length.
    #[must_use]
    pub fn variates(&self, x_row: &[f64], y_row: &[f64]) -> Vec<(f64, f64)> {
        assert_eq!(
            x_row.len(),
            self.x_variables,
            "the X row is {} values, not {}",
            x_row.len(),
            self.x_variables
        );
        assert_eq!(
            y_row.len(),
            self.y_variables,
            "the Y row is {} values, not {}",
            y_row.len(),
            self.y_variables
        );
        (0..self.pairs())
            .map(|j| {
                let mut u = 0.0;
                for (i, value) in x_row.iter().enumerate() {
                    u += self.x_weights[i * self.pairs() + j] * (value - self.x_mean(i));
                }
                let mut v = 0.0;
                for (i, value) in y_row.iter().enumerate() {
                    v += self.y_weights[i * self.pairs() + j] * (value - self.y_mean(i));
                }
                (u, v)
            })
            .collect()
    }

    /// The underlying generalised decomposition, for its diagnostics.
    #[must_use]
    pub const fn eigen(&self) -> &GeneralizedSymmetricEigen {
        &self.eigen
    }

    /// The largest `|lambda_j + lambda_{n-1-j}|` over the whole spectrum.
    ///
    /// Zero in exact arithmetic. See the module docs: nothing enforces it, so
    /// it is a free check on the reduction and the back-transform together.
    #[must_use]
    pub const fn spectrum_symmetry(&self) -> f64 {
        self.spectrum_symmetry
    }

    fn x_mean(&self, i: usize) -> f64 {
        self.x_means[i]
    }

    fn y_mean(&self, i: usize) -> f64 {
        self.y_means[i]
    }
}

/// Centre each column of a row-major matrix on its mean.
///
/// Compensated summation for the mean, as [`crate::pca::Pca::fit`] uses: the
/// naive running sum loses digits on data that sits far from zero, which every
/// TEP pressure and temperature channel does.
fn centre(data: &[f64], samples: usize, variables: usize) -> (Vec<f64>, Vec<f64>) {
    let mut means = vec![0.0; variables];
    for (v, m) in means.iter_mut().enumerate() {
        let mut sum = 0.0_f64;
        let mut compensation = 0.0_f64;
        for t in 0..samples {
            let adjusted = data[t * variables + v] - compensation;
            let next = sum + adjusted;
            compensation = (next - sum) - adjusted;
            sum = next;
        }
        *m = sum / samples as f64;
    }
    let mut out = vec![0.0; samples * variables];
    for t in 0..samples {
        for v in 0..variables {
            out[t * variables + v] = data[t * variables + v] - means[v];
        }
    }
    (out, means)
}

/// The dot product of column `i` of one centred matrix with column `j` of
/// another.
fn dot(
    left: &[f64],
    samples: usize,
    left_width: usize,
    i: usize,
    right: &[f64],
    right_width: usize,
    j: usize,
) -> f64 {
    let mut total = 0.0;
    for t in 0..samples {
        total += left[t * left_width + i] * right[t * right_width + j];
    }
    total
}

/// Split a time series into aligned past and future windows.
///
/// This is Larimore's construction, and it is what turns canonical correlation
/// analysis into a *dynamic* process monitor. At each time `t` the past block
/// holds the `past_lags` observations before `t` and the future block holds
/// the `future_leads` observations from `t` on. The canonical correlations
/// between the two blocks measure how much of the plant's future its recent
/// past predicts, which is the state dimension of the process; a fault that
/// changes the dynamics changes them even when it leaves every marginal
/// distribution alone.
///
/// `data` is `samples` by `variables`, row-major. The result is
/// `(past, future, rows)` with
///
/// ```text
/// past  row r, block l = observation (past_lags + r) - 1 - l    l = 0..past_lags
/// future row r, block m = observation (past_lags + r) + m       m = 0..future_leads
/// rows = samples + 1 - past_lags - future_leads
/// ```
///
/// so past block 0 is the most recent history and future block 0 is the
/// present. Most-recent-first on both sides, matching
/// [`crate::dpca::augment_with_lags`]'s present-first convention, so that the
/// two modules cannot disagree about which end of a row is which.
///
/// Returns empty matrices and zero rows when the record is too short to yield
/// a single window, rather than panicking: a horizon that is too short for a
/// chosen lag structure is a fact about the run, and a caller sweeping lag
/// counts should get zero rows rather than a crash.
///
/// # Panics
///
/// If `data.len()` is not `samples * variables`, or if either window length is
/// zero.
#[must_use]
pub fn past_future(
    data: &[f64],
    samples: usize,
    variables: usize,
    past_lags: usize,
    future_leads: usize,
) -> (Vec<f64>, Vec<f64>, usize) {
    assert_eq!(
        data.len(),
        samples * variables,
        "a {samples}-by-{variables} matrix needs {} entries, not {}",
        samples * variables,
        data.len()
    );
    assert!(
        past_lags > 0 && future_leads > 0,
        "both windows need at least one observation, got {past_lags} and {future_leads}"
    );

    let span = past_lags + future_leads;
    // `rows` below is `samples + 1 - span`, so this is the case where it would
    // be zero or negative.
    if samples < span {
        return (Vec::new(), Vec::new(), 0);
    }
    let rows = samples + 1 - span;
    let past_width = variables * past_lags;
    let future_width = variables * future_leads;
    let mut past = vec![0.0; rows * past_width];
    let mut future = vec![0.0; rows * future_width];
    for r in 0..rows {
        let present = past_lags + r;
        for l in 0..past_lags {
            let source = present - 1 - l;
            for v in 0..variables {
                past[r * past_width + l * variables + v] = data[source * variables + v];
            }
        }
        for m in 0..future_leads {
            let source = present + m;
            for v in 0..variables {
                future[r * future_width + m * variables + v] = data[source * variables + v];
            }
        }
    }
    (past, future, rows)
}
