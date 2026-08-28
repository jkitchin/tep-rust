//! Principal component analysis, with Hotelling's T-squared and the squared
//! prediction error.
//!
//! The detector Tier 6 is built around. A PCA model is fitted to fault-free
//! training data, and every later observation is scored twice: T-squared says
//! how far it moved *inside* the subspace the training data lived in, and SPE
//! says how far it moved *out* of that subspace. The two answer different
//! questions and both are needed. A fault that changes the operating point
//! without breaking the correlation structure raises T-squared alone; a fault
//! that breaks a relationship between variables raises SPE alone, often long
//! before it is visible in any single measurement.
//!
//! Jackson, J. E. (1991), *A User's Guide to Principal Components*, Wiley.
//!
//! Chiang, L. H., Russell, E. L. and Braatz, R. D. (2001), *Fault Detection and
//! Diagnosis in Industrial Systems*, Springer. Chapter 3 is the reference for
//! the two statistics and both control limits as the TEP literature uses them.
//!
//! # On the correlation matrix rather than the covariance matrix
//!
//! This always standardises: every column is centred and divided by its sample
//! standard deviation, so the matrix that gets diagonalised is the correlation
//! matrix. That is not a stylistic preference, it is forced by the data. The
//! TEP measurement vector mixes flows in kscmh, pressures in kPa near 2705,
//! temperatures in degrees C and compositions in mole percent. Covariance PCA on
//! that returns the pressure axis as the first principal component and learns
//! nothing else, because pressure's variance in kPa-squared dwarfs a mole
//! fraction's whatever the physics says. The `plot the eigenvalues and see`
//! failure mode is invisible: the model fits, the limits compute, and the
//! detector is a pressure alarm with extra steps.
//!
//! # Reporting
//!
//! [`Pca`] keeps the numbers, not only the model: the eigenvalue spectrum, the
//! variance each component explains, the columns it had to drop, and the
//! eigensolver's own residual. `CLAUDE.md`'s rule about recording numbers
//! rather than verdicts applies to a fitted model as much as to a log entry,
//! because a detector whose false alarm rate moved between two runs is
//! diagnosed from the spectrum and not from the alarm count.

use alloc::vec;
use alloc::vec::Vec;

use crate::distribution::{f_quantile, normal_quantile};
use crate::eigen::{SymmetricEigen, symmetric_eigen};
use crate::special::{not_positive, sqrt};

/// How many principal components to keep.
///
/// Named rather than a bare `k` because the rule is part of the model: two
/// detectors that retain a different number of components are different
/// detectors, and a Tier 6 report that does not say which rule produced its
/// numbers cannot be reproduced.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Retention {
    /// The smallest `k` whose components together explain at least this
    /// fraction of the total variance.
    ///
    /// The usual rule in the TEP literature, usually at 0.9. The fraction is a
    /// parameter because the answer is sensitive to it: on `d00` the difference
    /// between 0.85 and 0.95 is several components, and several components is
    /// the difference between a detector that catches fault 3 and one that does
    /// not.
    CumulativeVariance(f64),
    /// Every component whose eigenvalue exceeds one.
    ///
    /// Kaiser's rule. It means "keep a component only if it explains more than
    /// one original variable's worth of variance", which is a statement that
    /// only makes sense for correlation-matrix PCA, where the eigenvalues
    /// average exactly one. This crate never does covariance PCA, so the rule
    /// is always meaningful here.
    ///
    /// Kaiser, H. F. (1960), "The application of electronic computers to factor
    /// analysis", *Educational and Psychological Measurement* 20(1), 141-151.
    Kaiser,
    /// Exactly this many, clamped to the numerical rank.
    Fixed(usize),
}

/// The T-squared and SPE control limits for one confidence level.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ControlLimits {
    /// The upper control limit for Hotelling's T-squared.
    pub t_squared: f64,
    /// The upper control limit for the squared prediction error.
    pub spe: f64,
    /// The confidence the two were computed at, carried along so a report
    /// cannot quote a limit without its level.
    pub confidence: f64,
    /// How many components the model retained.
    pub components: usize,
    /// How many training samples the model was fitted to.
    pub samples: usize,
}

/// Both monitoring statistics for one observation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Statistics {
    /// Hotelling's T-squared, in the retained subspace.
    pub t_squared: f64,
    /// The squared prediction error, in the residual subspace.
    pub spe: f64,
}

/// A fitted PCA model.
#[derive(Clone, Debug, PartialEq)]
pub struct Pca {
    samples: usize,
    variables: usize,
    mean: Vec<f64>,
    sd: Vec<f64>,
    constant: Vec<usize>,
    eigenvalues: Vec<f64>,
    /// `variables` by `variables`, row-major. Column `j` is the loading vector
    /// of component `j`.
    loadings: Vec<f64>,
    retained: usize,
    rank: usize,
    total_variance: f64,
    eigen_sweeps: usize,
    eigen_off_diagonal_norm: f64,
    eigen_converged: bool,
}

impl Pca {
    /// Fit a model to a training matrix.
    ///
    /// `data` is `samples` by `variables`, row-major: `data[t * variables + v]`
    /// is variable `v` at sample `t`. That is the layout a run recorder
    /// produces and the layout the published `d00` files are in.
    ///
    /// # Standardisation
    ///
    /// Each column is centred on its mean and divided by its **sample**
    /// standard deviation, the `n - 1` one. The mean is formed by compensated
    /// summation and the sum of squares from the already-centred values, which
    /// is the two-pass computation Chan, Golub and LeVeque recommend when the
    /// data can be visited twice. The one-pass `sum(x^2) - n mean^2` form is
    /// what [`crate::summary`] exists to avoid, and it would be worse here than
    /// there: reactor pressure sits near 2705 with a spread of tenths, so the
    /// two terms agree to five digits before they are subtracted.
    ///
    /// # Constant columns
    ///
    /// A column with zero sample variance has no scale to standardise by. It is
    /// set to exactly zero for every observation, contributes nothing to either
    /// statistic, and its index is recorded in
    /// [`constant_columns`](Self::constant_columns). The alternative, keeping
    /// its raw deviation in the residual, would add a quantity in kPa or in
    /// degrees to a sum of standardised squares, which is not a quantity at all.
    /// TEP has such a column in ordinary operation: `XMV(12)`, the agitator
    /// speed, is held fixed by the base-case control scheme.
    ///
    /// A column whose variance is merely tiny rather than exactly zero is a
    /// real hazard and this does not guard against it, because any threshold
    /// would be a tuned number. Read
    /// [`standard_deviations`](Self::standard_deviations) and decide.
    ///
    /// # Panics
    ///
    /// If `data.len()` is not `samples * variables`, if `variables` is zero, or
    /// if there are fewer than two samples. Two is the fewest that has a sample
    /// variance at all.
    #[must_use]
    pub fn fit(data: &[f64], samples: usize, variables: usize, rule: Retention) -> Self {
        assert!(variables > 0, "a PCA model needs at least one variable");
        assert!(
            samples >= 2,
            "a sample standard deviation needs at least two samples, got {samples}"
        );
        assert_eq!(
            data.len(),
            samples * variables,
            "a {samples}-by-{variables} matrix needs {} entries, not {}",
            samples * variables,
            data.len()
        );

        // Pass one: the column means, compensated.
        let mut mean = vec![0.0; variables];
        for (v, m) in mean.iter_mut().enumerate() {
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

        // Pass two: the centred matrix, and each column's sum of squares.
        let mut centred = vec![0.0; samples * variables];
        for t in 0..samples {
            for v in 0..variables {
                centred[t * variables + v] = data[t * variables + v] - mean[v];
            }
        }
        let mut sum_squares = vec![0.0; variables];
        for (v, ss) in sum_squares.iter_mut().enumerate() {
            let mut total = 0.0;
            for t in 0..samples {
                let value = centred[t * variables + v];
                total += value * value;
            }
            *ss = total;
        }

        let divisor = (samples - 1) as f64;
        let sd: Vec<f64> = sum_squares.iter().map(|ss| sqrt(ss / divisor)).collect();
        let constant: Vec<usize> = (0..variables).filter(|&v| not_positive(sd[v])).collect();

        // The covariance matrix of the standardised data, which is the
        // correlation matrix wherever the column is not constant.
        //
        // Off-diagonal entries are `dot / sqrt(ss_i * ss_j)`, one square root of
        // a product rather than a product of two square roots, for the reason
        // spelled out in `correlation.rs`: on exactly collinear columns the
        // first form returns exactly 1 and the second returns one ulp more.
        // A test asserts the two agree bit for bit.
        let mut r = vec![0.0; variables * variables];
        for i in 0..variables {
            for j in i..variables {
                let value = if sum_squares[i] <= 0.0 || sum_squares[j] <= 0.0 {
                    // A constant column has zero variance, so its row and column
                    // of the covariance matrix of the standardised data are
                    // genuinely zero, diagonal included. That gives it an
                    // eigenvalue of exactly zero in its own direction, which the
                    // rank guard then declines to retain. This is the right
                    // answer, not a patch: the standardised column is
                    // identically zero and a zero vector has zero variance.
                    0.0
                } else if i == j {
                    // Exactly one. In exact arithmetic `dot / sqrt(ss * ss)` is
                    // one here; in `f64` it is one to within an ulp, and letting
                    // it be anything else would make the trace, which every
                    // variance fraction below is divided by, drift off the
                    // variable count.
                    1.0
                } else {
                    let mut dot = 0.0;
                    for t in 0..samples {
                        dot += centred[t * variables + i] * centred[t * variables + j];
                    }
                    (dot / sqrt(sum_squares[i] * sum_squares[j])).clamp(-1.0, 1.0)
                };
                r[i * variables + j] = value;
                r[j * variables + i] = value;
            }
        }

        let eigen = symmetric_eigen(&r, variables);
        let eigenvalues = eigen.values().to_vec();
        let loadings = eigen.vectors().to_vec();
        let total_variance: f64 = eigenvalues.iter().sum();
        let rank = numerical_rank(&eigen);
        let retained = retain(&eigenvalues, total_variance, rank, rule);

        Self {
            samples,
            variables,
            mean,
            sd,
            constant,
            eigenvalues,
            loadings,
            retained,
            rank,
            total_variance,
            eigen_sweeps: eigen.sweeps(),
            eigen_off_diagonal_norm: eigen.off_diagonal_norm(),
            eigen_converged: eigen.converged(),
        }
    }

    /// How many training samples.
    #[must_use]
    pub const fn samples(&self) -> usize {
        self.samples
    }

    /// How many variables.
    #[must_use]
    pub const fn variables(&self) -> usize {
        self.variables
    }

    /// How many components are retained.
    #[must_use]
    pub const fn retained(&self) -> usize {
        self.retained
    }

    /// The numerical rank of the correlation matrix.
    ///
    /// Components beyond this have eigenvalues indistinguishable from zero and
    /// are never retained, whatever the rule asked for: dividing a score by one
    /// of them in T-squared would turn rounding noise into an alarm.
    #[must_use]
    pub const fn rank(&self) -> usize {
        self.rank
    }

    /// The full eigenvalue spectrum, descending.
    #[must_use]
    pub fn eigenvalues(&self) -> &[f64] {
        &self.eigenvalues
    }

    /// The loadings: `variables` by `variables`, row-major, column `j` being
    /// component `j`.
    #[must_use]
    pub fn loadings(&self) -> &[f64] {
        &self.loadings
    }

    /// The column means the model centres on.
    #[must_use]
    pub fn means(&self) -> &[f64] {
        &self.mean
    }

    /// The column sample standard deviations the model scales by. Exactly zero
    /// for a constant column.
    #[must_use]
    pub fn standard_deviations(&self) -> &[f64] {
        &self.sd
    }

    /// The indices of the columns that did not move in training.
    #[must_use]
    pub fn constant_columns(&self) -> &[usize] {
        &self.constant
    }

    /// The sum of the eigenvalues, which is the trace of the correlation
    /// matrix: the number of variables that were not constant.
    #[must_use]
    pub const fn total_variance(&self) -> f64 {
        self.total_variance
    }

    /// The fraction of the total variance the retained components explain.
    #[must_use]
    pub fn explained_variance(&self) -> f64 {
        if not_positive(self.total_variance) {
            return f64::NAN;
        }
        let kept: f64 = self.eigenvalues[..self.retained].iter().sum();
        kept / self.total_variance
    }

    /// How many Jacobi sweeps the fit needed, and whether it converged.
    ///
    /// Returns `(sweeps, off_diagonal_norm, converged)`. Carried out of the
    /// eigensolver so a fit can be reported without a second decomposition.
    #[must_use]
    pub const fn eigen_residual(&self) -> (usize, f64, bool) {
        (
            self.eigen_sweeps,
            self.eigen_off_diagonal_norm,
            self.eigen_converged,
        )
    }

    /// Standardise one observation with the training mean and scale.
    ///
    /// # Panics
    ///
    /// If `sample.len()` is not `variables`.
    #[must_use]
    pub fn standardise(&self, sample: &[f64]) -> Vec<f64> {
        assert_eq!(
            sample.len(),
            self.variables,
            "this model takes {} variables, not {}",
            self.variables,
            sample.len()
        );
        (0..self.variables)
            .map(|v| {
                if self.sd[v] > 0.0 {
                    (sample[v] - self.mean[v]) / self.sd[v]
                } else {
                    0.0
                }
            })
            .collect()
    }

    /// The scores of one observation on every component, `P' z`.
    ///
    /// All of them, not only the retained ones: the discarded scores are what
    /// SPE is made of, and a caller diagnosing a fault wants to see where in the
    /// residual subspace it landed.
    ///
    /// # Panics
    ///
    /// If `sample.len()` is not `variables`.
    #[must_use]
    pub fn scores(&self, sample: &[f64]) -> Vec<f64> {
        let z = self.standardise(sample);
        (0..self.variables).map(|j| self.score(&z, j)).collect()
    }

    /// Hotelling's T-squared for one observation.
    ///
    /// ```text
    /// T^2 = sum_{j < k} (p_j' z)^2 / lambda_j
    /// ```
    ///
    /// The squared Mahalanobis distance from the training mean, measured in the
    /// retained subspace only. Each score is divided by its own eigenvalue, so
    /// a component that barely moved in training weighs a movement more heavily
    /// than one that swung widely, which is what makes this a distance rather
    /// than a sum of squares.
    ///
    /// Zero when no component is retained.
    ///
    /// # Panics
    ///
    /// If `sample.len()` is not `variables`.
    #[must_use]
    pub fn t_squared(&self, sample: &[f64]) -> f64 {
        let z = self.standardise(sample);
        self.t_squared_standardised(&z)
    }

    /// The squared prediction error of one observation, also called Q.
    ///
    /// ```text
    /// e   = z - P_k P_k' z
    /// SPE = e' e
    /// ```
    ///
    /// The squared distance from the retained subspace: how much of the
    /// observation the model cannot explain.
    ///
    /// When every component is retained the subspace is the whole space and
    /// this is zero to rounding, not exactly zero. The residual is formed by
    /// subtracting the reconstruction, so what survives is the `O(eps ||z||)`
    /// by which the loadings fall short of exactly orthonormal, squared: about
    /// 1e-31 for a standardised observation. Exactness is available by summing
    /// the discarded scores instead, and it would be exactness about the wrong
    /// quantity. See the note on the implementation below.
    ///
    /// # Panics
    ///
    /// If `sample.len()` is not `variables`.
    #[must_use]
    pub fn spe(&self, sample: &[f64]) -> f64 {
        let z = self.standardise(sample);
        self.spe_standardised(&z)
    }

    /// Both statistics, standardising once instead of twice.
    ///
    /// # Panics
    ///
    /// If `sample.len()` is not `variables`.
    #[must_use]
    pub fn statistics(&self, sample: &[f64]) -> Statistics {
        let z = self.standardise(sample);
        Statistics {
            t_squared: self.t_squared_standardised(&z),
            spe: self.spe_standardised(&z),
        }
    }

    /// Both control limits at one confidence level.
    #[must_use]
    pub fn limits(&self, confidence: f64) -> ControlLimits {
        ControlLimits {
            t_squared: t_squared_limit(self.retained, self.samples, confidence),
            spe: spe_limit(&self.eigenvalues[self.retained..], confidence),
            confidence,
            components: self.retained,
            samples: self.samples,
        }
    }

    fn score(&self, standardised: &[f64], j: usize) -> f64 {
        let mut total = 0.0;
        for (i, z) in standardised.iter().enumerate() {
            total += self.loadings[i * self.variables + j] * z;
        }
        total
    }

    fn t_squared_standardised(&self, z: &[f64]) -> f64 {
        let mut total = 0.0;
        for j in 0..self.retained {
            let score = self.score(z, j);
            total += score * score / self.eigenvalues[j];
        }
        total
    }

    fn spe_standardised(&self, z: &[f64]) -> f64 {
        // Formed by subtracting the retained projection from `z` rather than by
        // summing the discarded scores. The two are equal in exact arithmetic
        // and the first is the one that stays right when the loadings are only
        // orthonormal to rounding: it measures what the model actually failed to
        // reconstruct, which is the definition, rather than what an assumed
        // completeness relation says it should have failed to reconstruct.
        let mut residual = z.to_vec();
        for j in 0..self.retained {
            let score = self.score(z, j);
            for (i, e) in residual.iter_mut().enumerate() {
                *e -= score * self.loadings[i * self.variables + j];
            }
        }
        residual.iter().map(|e| e * e).sum()
    }
}

/// The numerical rank: how many eigenvalues are distinguishable from zero.
///
/// The threshold is `n * eps * lambda_max`, LAPACK's convention for the rank of
/// a matrix from its singular values, and not a tuned number.
fn numerical_rank(eigen: &SymmetricEigen) -> usize {
    let values = eigen.values();
    let Some(&largest) = values.first() else {
        return 0;
    };
    if not_positive(largest) {
        return 0;
    }
    let threshold = values.len() as f64 * f64::EPSILON * largest;
    values.iter().filter(|&&v| v > threshold).count()
}

/// Apply a retention rule to a spectrum.
fn retain(eigenvalues: &[f64], total: f64, rank: usize, rule: Retention) -> usize {
    let wanted = match rule {
        Retention::CumulativeVariance(fraction) => {
            if not_positive(total) {
                0
            } else {
                let mut running = 0.0;
                let mut k = 0;
                while k < eigenvalues.len() && running / total < fraction {
                    running += eigenvalues[k];
                    k += 1;
                }
                k
            }
        }
        Retention::Kaiser => eigenvalues.iter().filter(|&&v| v > 1.0).count(),
        Retention::Fixed(k) => k,
    };
    // Never past the numerical rank, whatever the rule asked for: T-squared
    // divides by the eigenvalue, and a component at the noise floor turns
    // rounding into a fault.
    wanted.min(rank)
}

/// The upper control limit for Hotelling's T-squared on a **new** observation.
///
/// ```text
/// T^2_alpha = k (n + 1)(n - 1) / (n (n - k)) * F_alpha(k, n - k)
/// ```
///
/// `k` components, `n` training samples, `alpha` the confidence.
///
/// # Which limit this is
///
/// There are two, and they are not interchangeable. This is the one for an
/// observation *independent of* the training set, the case a fault-detection
/// experiment is in: the model is fitted to `d00` and then applied to `d01`
/// through `d21`. The other, for an observation that was itself part of the
/// training set, is `(n - 1)^2 / n` times a Beta quantile and is what a phase-I
/// screen of the training data would use. The two differ by a factor of about
/// `(n + 1) / (n - k)` at small `n` and converge as `n` grows.
///
/// Johnson, R. A. and Wichern, D. W. (2007), *Applied Multivariate Statistical
/// Analysis*, 6th ed., section 5.6, for the distribution of the future
/// observation. Chiang, Russell and Braatz (2001) eq. (3.7) states the same
/// limit in PCA notation, which is the form quoted throughout the TEP
/// literature.
///
/// # Behaviour at the edges
///
/// Zero components gives zero, which is consistent: with nothing retained
/// T-squared is identically zero and never exceeds it. `n <= k` gives `NaN`,
/// because the F distribution has no second degrees-of-freedom argument there
/// and a model with as many components as training samples has not been fitted,
/// it has been memorised.
///
/// # Large `n`
///
/// `k F(k, n - k)` tends to a chi-squared with `k` degrees of freedom, and the
/// leading factor tends to one, so the limit tends to the chi-squared quantile.
/// That is the check the tests make against published chi-squared tables.
#[must_use]
pub fn t_squared_limit(components: usize, samples: usize, confidence: f64) -> f64 {
    if components == 0 {
        return 0.0;
    }
    if samples <= components || !(0.0..1.0).contains(&confidence) || confidence <= 0.0 {
        return f64::NAN;
    }
    let k = components as f64;
    let n = samples as f64;
    let f = f_quantile(confidence, k, n - k);
    k * (n + 1.0) * (n - 1.0) / (n * (n - k)) * f
}

/// The upper control limit for SPE, by the Jackson-Mudholkar approximation.
///
/// `residual_eigenvalues` is the discarded tail of the spectrum, the
/// eigenvalues of the components the model did **not** retain.
///
/// ```text
/// theta_i = sum_j lambda_j^i                       i = 1, 2, 3
/// h0      = 1 - 2 theta1 theta3 / (3 theta2^2)
/// SPE_alpha = theta1 [ c_alpha sqrt(2 theta2 h0^2) / theta1
///                      + 1
///                      + theta2 h0 (h0 - 1) / theta1^2 ]^(1/h0)
/// ```
///
/// with `c_alpha` the standard normal deviate at the confidence level.
///
/// Jackson, J. E. and Mudholkar, G. S. (1979), "Control procedures for
/// residuals associated with principal component analysis", *Technometrics*
/// 21(3), 341-349.
///
/// # Why an approximation at all
///
/// SPE is a weighted sum of squared normals, `sum lambda_j chi^2_1`, whose exact
/// distribution has no closed form. Jackson and Mudholkar match the first three
/// moments with a power transformation, which is why only the first three power
/// sums of the residual eigenvalues appear. The approximation is exact in one
/// case worth knowing: when the residual eigenvalues are all equal, `h0` is
/// exactly `1/3` and the formula collapses to the Wilson-Hilferty cube-root
/// approximation of the chi-squared quantile with `p - k` degrees of freedom,
/// which is accurate to about a part in a thousand. The tests check that
/// collapse against published chi-squared tables, which is the only way to check
/// this formula without a second implementation of it.
///
/// Wilson, E. B. and Hilferty, M. M. (1931), "The distribution of chi-square",
/// *PNAS* 17(12), 684-688.
///
/// # Behaviour at the edges
///
/// An empty residual, or one whose eigenvalues are all zero, gives exactly
/// zero: SPE is then identically zero and cannot exceed it. `confidence` at or
/// below one half gives a negative `c_alpha`, which can drive the bracket
/// negative and produce `NaN` from a fractional power; a control limit below the
/// median is not a thing anyone wants, and `NaN` says so rather than returning a
/// number.
///
/// Eigenvalues slightly below zero are clamped to zero. A correlation matrix is
/// positive semi-definite in exact arithmetic, so a negative eigenvalue is
/// rounding noise on a rank-deficient matrix; leaving it in would put a negative
/// cube into `theta3` and bias `h0`.
#[must_use]
pub fn spe_limit(residual_eigenvalues: &[f64], confidence: f64) -> f64 {
    let mut theta1 = 0.0;
    let mut theta2 = 0.0;
    let mut theta3 = 0.0;
    for &lambda in residual_eigenvalues {
        let l = if lambda > 0.0 { lambda } else { 0.0 };
        theta1 += l;
        theta2 += l * l;
        theta3 += l * l * l;
    }
    if not_positive(theta1) {
        return 0.0;
    }
    if !(0.0..1.0).contains(&confidence) {
        return f64::NAN;
    }

    let h0 = 1.0 - 2.0 * theta1 * theta3 / (3.0 * theta2 * theta2);
    if not_positive(libm::fabs(h0)) {
        // The exponent is `1 / h0`. There is no answer here, and returning the
        // untransformed bracket would be a different approximation wearing this
        // one's name.
        return f64::NAN;
    }
    let c = normal_quantile(confidence);
    let bracket = c * sqrt(2.0 * theta2 * h0 * h0) / theta1
        + 1.0
        + theta2 * h0 * (h0 - 1.0) / (theta1 * theta1);
    theta1 * libm::pow(bracket, 1.0 / h0)
}
