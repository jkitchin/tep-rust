//! Known-answer tests for the cyclic Jacobi eigensolver.
//!
//! Same rule as `known_answers.rs`: **nothing is checked against a number this
//! project produced**. Every assertion here is against a closed form, an exact
//! identity, or a published spectrum.
//!
//! The battery is deliberately layered, because an eigensolver has several
//! independent ways to be wrong and no single check catches them all:
//!
//! - **Closed forms.** A diagonal matrix is its own decomposition; a 2-by-2 has
//!   a quadratic formula; a Hadamard rotation of a known diagonal has exactly
//!   representable eigenvectors. These pin the values *and* the vectors.
//! - **A published spectrum.** The 1-D Laplacian's eigenvalues and eigenvectors
//!   are known in closed form for every order, and they are dense, unequally
//!   spaced and clustered at one end, which is the hard case.
//! - **Invariants.** The eigenvalues sum to the trace and multiply to the
//!   determinant, and the determinant is computed here from its own cofactor
//!   formula rather than from the decomposition.
//! - **Reconstruction.** `A = V D V'` and `V' V = I` to machine precision. This
//!   is the strongest single check and it needs no reference at all, which is
//!   why it is applied to matrices that have no closed form: a random symmetric
//!   matrix and an 8-by-8 Hilbert matrix.

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

use tepsim_stats::{SymmetricEigen, symmetric_eigen};

/// A deterministic generator, so these tests are reproducible without pulling
/// in a random-number crate. Not a good generator; good enough to make a
/// scatter that is not a pattern.
struct Lcg(u64);

impl Lcg {
    fn new(seed: u64) -> Self {
        Self(seed | 1)
    }
    /// Uniform on [-0.5, 0.5).
    fn next(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1);
        ((self.0 >> 11) as f64) / ((1_u64 << 53) as f64) - 0.5
    }
}

/// Build a symmetric matrix from the upper triangle of a generator's output.
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

/// The largest absolute difference between the decomposition and its input.
fn reconstruction_error(a: &[f64], eigen: &SymmetricEigen) -> f64 {
    let rebuilt = eigen.reconstruct();
    a.iter()
        .zip(&rebuilt)
        .map(|(x, y)| (x - y).abs())
        .fold(0.0_f64, f64::max)
}

/// The largest absolute deviation of `V' V` from the identity.
fn orthonormality_error(eigen: &SymmetricEigen) -> f64 {
    let n = eigen.n();
    let v = eigen.vectors();
    let mut worst = 0.0_f64;
    for i in 0..n {
        for j in 0..n {
            let mut dot = 0.0;
            for k in 0..n {
                dot += v[k * n + i] * v[k * n + j];
            }
            let expected = if i == j { 1.0 } else { 0.0 };
            worst = worst.max((dot - expected).abs());
        }
    }
    worst
}

// ---------------------------------------------------------------------------
// Closed forms
// ---------------------------------------------------------------------------

/// A diagonal matrix is already diagonal, so the answer is itself.
///
/// Exactly itself: no rotation should fire at all, because every off-diagonal
/// entry starts at zero and the rotation declines to act on a zero pivot. The
/// eigenvalues must therefore be the diagonal entries **bit for bit**, sorted,
/// and the eigenvectors the signed unit basis vectors.
#[test]
fn a_diagonal_matrix_is_its_own_decomposition() {
    let diagonal = [3.0, 1.0, 7.0, -2.0];
    let n = diagonal.len();
    let mut a = vec![0.0; n * n];
    for i in 0..n {
        a[i * n + i] = diagonal[i];
    }

    let eigen = symmetric_eigen(&a, n);
    assert_eq!(eigen.sweeps(), 0, "a diagonal matrix needs no sweeps");
    assert_eq!(eigen.off_diagonal_norm(), 0.0);
    assert!(eigen.converged());
    assert_eq!(eigen.values(), &[7.0, 3.0, 1.0, -2.0]);

    // Column `j` of V is the basis vector for whichever original position the
    // `j`-th largest eigenvalue came from: 7 was at index 2, 3 at index 0, 1 at
    // index 1, -2 at index 3.
    let source = [2, 0, 1, 3];
    for (j, &i) in source.iter().enumerate() {
        for row in 0..n {
            let expected = if row == i { 1.0 } else { 0.0 };
            assert_eq!(
                eigen.component(row, j),
                expected,
                "V[{row}][{j}] for eigenvalue {}",
                eigen.values()[j]
            );
        }
    }
}

/// The 2-by-2 closed form, over a grid that includes the hard cases.
///
/// ```text
/// [ a  b ]  has eigenvalues  (a + d)/2 +- sqrt( ((a - d)/2)^2 + b^2 )
/// [ b  d ]
/// ```
///
/// # The eigenvector reference, and why not the obvious one
///
/// The textbook eigenvector for `lambda` is `(b, lambda - a)`, read off the
/// first row. It is useless as a reference here, and finding that out is worth
/// recording. Take `a = -3`, `d = -0.5`, `b = 1e-8`: the smaller eigenvalue is
/// `-3 - 4e-17`, which rounds to exactly `-3` because an ulp there is `4.4e-16`.
/// So `lambda - a` evaluates to exactly zero and the reference direction comes
/// out `(1, 0)` when the true eigenvector is `(1, -4e-9)`. The second row's form
/// `(lambda - d, b)` fails the same way in the mirror case, and when the two
/// diagonal entries are equal *both* forms lose the small component, because
/// `lambda - a` is then `1e-8` sitting on a magnitude of 3.
///
/// The cure is a form that never subtracts a small number from a large one. The
/// diagonalising rotation angle satisfies
///
/// ```text
/// tan(2 theta) = 2b / (a - d)
/// ```
///
/// and `(cos theta, sin theta)` is the eigenvector of the **larger** eigenvalue,
/// `(-sin theta, cos theta)` that of the smaller. With `theta` from `atan2`
/// there is no cancellation anywhere: the inputs are `2b` and `a - d`, both
/// exact, and `atan2` is well conditioned away from the origin. The origin is
/// `a = d` and `b = 0`, a multiple of the identity, where every vector is an
/// eigenvector and there is nothing to check.
///
/// The comparison is `|<v, reference>| = 1`, which is sign insensitive. The sign
/// convention is pinned by its own test.
#[test]
fn the_two_by_two_matches_the_closed_form() {
    let mut worst_value = 0.0_f64;
    let mut worst_direction = 0.0_f64;
    let mut checked = 0;
    for a in [-3.0_f64, -0.5, 0.0, 1.0, 5.0, 1e6] {
        for d in [-3.0_f64, -0.5, 0.0, 1.0, 5.0, 1e6] {
            for b in [-7.0_f64, -1.0, -1e-8, 1e-8, 0.25, 2.0, 1e6] {
                let matrix = [a, b, b, d];
                let eigen = symmetric_eigen(&matrix, 2);
                assert!(eigen.converged(), "a={a} b={b} d={d}");

                let middle = 0.5 * (a + d);
                let half_gap = 0.5 * (a - d);
                // `hypot` rather than `sqrt(h*h + b*b)`: it is the more
                // accurate of the two, and a reference should be the accurate
                // one. It also overflows nowhere, which the squared form does
                // above 1e154.
                let spread = half_gap.hypot(b);
                let expected = [middle + spread, middle - spread];

                let theta = 0.5 * (2.0 * b).atan2(a - d);
                let (sin, cos) = theta.sin_cos();
                let reference = [(cos, sin), (-sin, cos)];

                for (j, &want) in expected.iter().enumerate() {
                    let got = eigen.values()[j];
                    // Against the scale of the terms the closed form adds up,
                    // not against the eigenvalue: `middle - spread` cancels when
                    // one eigenvalue is near zero, and a relative bound there
                    // would be measuring `f64` subtraction rather than the
                    // solver. That is the same distinction Tier 2 draws.
                    let scale = spread.abs().max(middle.abs()).max(f64::MIN_POSITIVE);
                    worst_value = worst_value.max((got - want).abs() / scale);

                    let dot = eigen.component(0, j) * reference[j].0
                        + eigen.component(1, j) * reference[j].1;
                    worst_direction = worst_direction.max((dot.abs() - 1.0).abs());
                    checked += 1;
                }
            }
        }
    }
    println!(
        "2x2 closed form: {checked} eigenpairs, worst relative eigenvalue error \
         {worst_value:.3e}, worst |<v, closed form>| - 1 = {worst_direction:.3e}"
    );
    // The routine claims about fifteen digits; this asks for thirteen, the same
    // bar `known_answers.rs` sets.
    assert!(worst_value < 1e-13, "{worst_value:.3e}");
    assert!(worst_direction < 1e-13, "{worst_direction:.3e}");
}

/// A known orthogonal rotation of a known diagonal, with every entry exactly
/// representable.
///
/// The normalised 4-by-4 Hadamard matrix `H / 2` has entries `+-1/2`, is exactly
/// orthogonal in `f64`, and `H/2 * diag(8, 4, 2, 1) * (H/2)'` has entries that
/// are quarters, which `f64` also holds exactly. So the input carries no
/// rounding at all and the answer is known in full: eigenvalues 8, 4, 2, 1 and
/// eigenvectors the columns of `H / 2`.
///
/// This is the test the task description calls "a known rotation of a diagonal
/// matrix has known eigenvalues", made exact.
#[test]
fn a_hadamard_rotation_of_a_known_diagonal_returns_both() {
    const N: usize = 4;
    // Columns are mutually orthogonal by construction; this is the Sylvester
    // Hadamard matrix of order 4.
    let h = [
        [1.0, 1.0, 1.0, 1.0],
        [1.0, -1.0, 1.0, -1.0],
        [1.0, 1.0, -1.0, -1.0],
        [1.0, -1.0, -1.0, 1.0],
    ];
    let d = [8.0, 4.0, 2.0, 1.0];

    let mut a = vec![0.0; N * N];
    for i in 0..N {
        for j in 0..N {
            let mut total = 0.0;
            for k in 0..N {
                total += 0.5 * h[i][k] * d[k] * 0.5 * h[j][k];
            }
            a[i * N + j] = total;
        }
    }
    // Every entry is a quarter-integer, so the matrix fed to the solver is
    // exact. Check that rather than assuming it.
    for value in &a {
        let quarters: f64 = value * 4.0;
        assert_eq!(quarters, quarters.round(), "{value} is not a quarter");
    }

    let eigen = symmetric_eigen(&a, N);
    assert!(eigen.converged());
    let mut worst_value = 0.0_f64;
    let mut worst_vector = 0.0_f64;
    for j in 0..N {
        worst_value = worst_value.max((eigen.values()[j] - d[j]).abs() / d[j]);
        // Compared up to sign, deliberately. Every component of a Hadamard
        // column has the same magnitude, so the sign convention -- make the
        // largest-magnitude component positive -- has a four-way tie to break,
        // and which component wins it is decided by whether the rotations left
        // that entry at 0.5 or at 0.49999999999999994. The convention is still
        // reproducible for a given input, which is what it is for, and that is
        // what `the_largest_component_of_every_eigenvector_is_positive` checks.
        // It is not stable under an ulp of perturbation, and a test that
        // pretended otherwise would be testing the rounding.
        let mut error_plus = 0.0_f64;
        let mut error_minus = 0.0_f64;
        for i in 0..N {
            error_plus = error_plus.max((eigen.component(i, j) - 0.5 * h[i][j]).abs());
            error_minus = error_minus.max((eigen.component(i, j) + 0.5 * h[i][j]).abs());
        }
        worst_vector = worst_vector.max(error_plus.min(error_minus));
    }
    println!(
        "Hadamard rotation: worst relative eigenvalue error {worst_value:.3e}, \
         worst eigenvector component error {worst_vector:.3e}, \
         {} sweeps",
        eigen.sweeps()
    );
    assert!(worst_value < 1e-14, "{worst_value:.3e}");
    assert!(worst_vector < 1e-14, "{worst_vector:.3e}");
}

/// The published closed form for the 1-D Laplacian's spectrum.
///
/// The `n`-by-`n` tridiagonal matrix with 2 on the diagonal and -1 beside it
/// has
///
/// ```text
/// lambda_k = 2 - 2 cos(k pi / (n + 1)) = 4 sin^2(k pi / (2 (n + 1)))
/// v_k(j)   = sin(j k pi / (n + 1))                  j = 1..n
/// ```
///
/// for `k = 1..n`. Standard, and in every numerical analysis text; it is the
/// discrete sine transform's diagonalisation. It is a far more demanding test
/// than a 4-by-4: at `n = 12` the eigenvalues run from 0.058 to 3.94 with no
/// gaps, and the ones at the small end are clustered, which is where a solver
/// that is only absolutely accurate loses them.
///
/// Vectors are compared by the absolute dot product with the normalised closed
/// form, which is 1 for a match in either sign.
#[test]
fn the_laplacian_spectrum_matches_the_published_closed_form() {
    for n in [3_usize, 7, 12, 25] {
        let mut a = vec![0.0; n * n];
        for i in 0..n {
            a[i * n + i] = 2.0;
            if i + 1 < n {
                a[i * n + i + 1] = -1.0;
                a[(i + 1) * n + i] = -1.0;
            }
        }

        let eigen = symmetric_eigen(&a, n);
        assert!(eigen.converged(), "n={n}");

        let pi = std::f64::consts::PI;
        let mut worst_value = 0.0_f64;
        let mut worst_vector = 0.0_f64;
        for position in 0..n {
            // Descending order, so position 0 is k = n.
            let k = (n - position) as f64;
            let angle = k * pi / (n + 1) as f64;
            let expected = 2.0 - 2.0 * angle.cos();
            worst_value = worst_value.max((eigen.values()[position] - expected).abs() / expected);

            // The closed-form eigenvector, normalised: the sum of
            // sin^2(j k pi / (n+1)) over j = 1..n is exactly (n+1)/2.
            let scale = (2.0 / (n + 1) as f64).sqrt();
            let mut dot = 0.0;
            for j in 0..n {
                let component = scale * (((j + 1) as f64) * angle).sin();
                dot += component * eigen.component(j, position);
            }
            worst_vector = worst_vector.max((dot.abs() - 1.0).abs());
        }
        println!(
            "Laplacian n={n}: {} sweeps, worst relative eigenvalue error \
             {worst_value:.3e}, worst |<v, closed form>| - 1 = {worst_vector:.3e}, \
             off(A) = {:.3e} against ||A||_F = {:.3e}",
            eigen.sweeps(),
            eigen.off_diagonal_norm(),
            eigen.frobenius_norm()
        );
        assert!(worst_value < 1e-13, "n={n}: {worst_value:.3e}");
        assert!(worst_vector < 1e-13, "n={n}: {worst_vector:.3e}");
    }
}

// ---------------------------------------------------------------------------
// Invariants
// ---------------------------------------------------------------------------

/// The eigenvalues sum to the trace and multiply to the determinant.
///
/// The determinant here is the 3-by-3 cofactor formula written out, not
/// anything this crate computes, so the two sides are genuinely independent.
#[test]
fn the_eigenvalues_reproduce_the_trace_and_the_determinant() {
    let mut rng = Lcg::new(0x1AC0B1);
    let mut worst_trace = 0.0_f64;
    let mut worst_determinant = 0.0_f64;
    for _ in 0..200 {
        let a = random_symmetric(&mut rng, 3, 4.0);
        let eigen = symmetric_eigen(&a, 3);
        assert!(eigen.converged());

        let trace = a[0] + a[4] + a[8];
        let sum: f64 = eigen.values().iter().sum();
        let scale = a.iter().map(|x| x.abs()).fold(0.0_f64, f64::max);
        worst_trace = worst_trace.max((sum - trace).abs() / scale);

        let determinant = a[0] * (a[4] * a[8] - a[5] * a[7]) - a[1] * (a[3] * a[8] - a[5] * a[6])
            + a[2] * (a[3] * a[7] - a[4] * a[6]);
        let product: f64 = eigen.values().iter().product();
        worst_determinant =
            worst_determinant.max((product - determinant).abs() / scale.powi(3).max(1e-300));
    }
    println!(
        "trace and determinant over 200 random 3x3: worst scaled trace error \
         {worst_trace:.3e}, worst scaled determinant error {worst_determinant:.3e}"
    );
    assert!(worst_trace < 1e-14, "{worst_trace:.3e}");
    // Looser than the trace, and not because the solver is worse: the
    // determinant of a random matrix can be many orders below the scale of its
    // entries, and both sides then lose the same absolute precision. The bound
    // is on the error *relative to the matrix scale cubed*, which is the size
    // of the terms the cofactor expansion actually adds up, not the size of
    // their cancelling sum. That is the same distinction `CLAUDE.md` draws for
    // Tier 2 balances.
    assert!(worst_determinant < 1e-14, "{worst_determinant:.3e}");
}

/// `A = V D V'` and `V' V = I`, over matrices with no closed form at all.
///
/// The Hilbert matrix is in here on purpose. `H_ij = 1 / (i + j + 1)` is
/// positive definite and famously ill conditioned: at order 8 the condition
/// number is above 1e10, so the smallest eigenvalue is ten orders below the
/// largest. Jacobi's claim, and the reason this module uses it rather than a
/// tridiagonal QR, is that it computes such an eigenvalue to high relative
/// accuracy and does not return it as a small negative number. That the
/// smallest eigenvalue comes out positive is a mathematical fact about the
/// Hilbert matrix, not a number this project produced.
#[test]
fn the_decomposition_reconstructs_and_stays_orthonormal() {
    let mut rng = Lcg::new(0xE16E);

    let mut cases: Vec<(String, usize, Vec<f64>)> = Vec::new();
    for n in [1_usize, 2, 5, 13, 30] {
        cases.push((
            format!("random {n}x{n}"),
            n,
            random_symmetric(&mut rng, n, 10.0),
        ));
    }
    // A matrix whose entries span twelve orders, so the small eigenvalues are
    // nowhere near the large ones.
    let mut graded = random_symmetric(&mut rng, 8, 1.0);
    for i in 0..8 {
        for j in 0..8 {
            let factor = 10f64.powi(-3 * (i.max(j) as i32));
            graded[i * 8 + j] *= factor;
        }
    }
    cases.push(("graded 8x8".into(), 8, graded));

    let mut hilbert = vec![0.0; 64];
    for i in 0..8 {
        for j in 0..8 {
            hilbert[i * 8 + j] = 1.0 / (i + j + 1) as f64;
        }
    }
    cases.push(("Hilbert 8x8".into(), 8, hilbert.clone()));

    // Degenerate spectra: every eigenvalue equal, and a rank-one matrix.
    let mut identity = vec![0.0; 36];
    for i in 0..6 {
        identity[i * 6 + i] = 1.0;
    }
    cases.push(("identity 6x6".into(), 6, identity));
    cases.push(("zero 4x4".into(), 4, vec![0.0; 16]));
    let u = [0.3, -0.5, 1.2, 0.7, -0.9];
    let mut rank_one = vec![0.0; 25];
    for i in 0..5 {
        for j in 0..5 {
            rank_one[i * 5 + j] = u[i] * u[j];
        }
    }
    cases.push(("rank-one 5x5".into(), 5, rank_one));

    for (name, n, a) in cases {
        let eigen = symmetric_eigen(&a, n);
        let scale = a.iter().map(|x| x.abs()).fold(0.0_f64, f64::max).max(1.0);
        let reconstruction = reconstruction_error(&a, &eigen) / scale;
        let orthonormality = orthonormality_error(&eigen);
        println!(
            "{name}: {} sweeps, ||A - VDV'||_max / scale = {reconstruction:.3e}, \
             ||V'V - I||_max = {orthonormality:.3e}, converged = {}",
            eigen.sweeps(),
            eigen.converged()
        );
        assert!(eigen.converged(), "{name} did not converge");
        // Backward stability for Jacobi is `O(n) eps`; at n = 30 that is 6.7e-15
        // and this asks for 1e-13, which leaves a factor of fifteen and no room
        // for an algorithmic error.
        assert!(reconstruction < 1e-13, "{name}: {reconstruction:.3e}");
        assert!(orthonormality < 1e-13, "{name}: {orthonormality:.3e}");
        // Descending, always.
        for pair in eigen.values().windows(2) {
            assert!(pair[0] >= pair[1], "{name}: {:?}", eigen.values());
        }
    }

    // And the Hilbert matrix's positive definiteness, which is the property
    // this solver was chosen for.
    let eigen = symmetric_eigen(&hilbert, 8);
    let smallest = eigen.values()[7];
    println!(
        "Hilbert 8x8: largest {:.6e}, smallest {:.6e}, condition number {:.3e}",
        eigen.values()[0],
        smallest,
        eigen.values()[0] / smallest
    );
    assert!(
        smallest > 0.0,
        "the Hilbert matrix is positive definite; got {smallest:.6e}"
    );
    assert!(
        eigen.values()[0] / smallest > 1e9,
        "this Hilbert matrix is not ill conditioned enough to demonstrate \
         anything; condition number {:.3e}",
        eigen.values()[0] / smallest
    );
}

// ---------------------------------------------------------------------------
// Contract
// ---------------------------------------------------------------------------

/// The pivot is annihilated **exactly**, not to within rounding.
///
/// The rotation is chosen so that `a[p][q]` becomes zero, so the routine writes
/// zero there rather than evaluating `c s (a_pp - a_qq) + (c^2 - s^2) a_pq` and
/// keeping whatever the arithmetic returns. That is what makes `off(A)` decrease
/// monotonically, which is what makes the sweeps terminate: an entry that is
/// re-dirtied by a few ulp each time it is visited can be re-dirtied for ever.
///
/// A 2-by-2 is the case where the difference is visible, because one rotation
/// finishes it and `off(A)` is then exactly the pivot. Evaluating the expression
/// instead leaves `1e-16` there, which changes no eigenvalue any other test in
/// this file looks at.
///
/// # Provenance
///
/// This test exists because it was missing. A mutation run that flipped the
/// assignment to the evaluated form produced no failures at all across the
/// crate's 136 tests, which is the definition of an untested line.
#[test]
fn a_rotation_annihilates_its_pivot_exactly() {
    let mut rng = Lcg::new(0xA1_1A);
    for case in 0..200 {
        let a = random_symmetric(&mut rng, 2, 5.0);
        let eigen = symmetric_eigen(&a, 2);
        assert_eq!(
            eigen.off_diagonal_norm(),
            0.0,
            "case {case}: one rotation should leave nothing off the diagonal of \
             {a:?}, and it left {:.3e}",
            eigen.off_diagonal_norm()
        );
    }
    // And where the two diagonal entries are decades apart, which is where an
    // evaluated pivot would be largest relative to nothing.
    for (app, aqq, apq) in [
        (1e12_f64, 1e-12_f64, 1.0_f64),
        (1.0, 1e-15, 1e-7),
        (1e8, -1e8, 3.25),
        (-7.5, -7.5, 2.0),
    ] {
        let eigen = symmetric_eigen(&[app, apq, apq, aqq], 2);
        assert_eq!(
            eigen.off_diagonal_norm(),
            0.0,
            "([{app}, {apq}], [{apq}, {aqq}]) left {:.3e} off the diagonal",
            eigen.off_diagonal_norm()
        );
    }
}

/// The sign convention, stated as a test so it cannot drift.
#[test]
fn the_largest_component_of_every_eigenvector_is_positive() {
    let mut rng = Lcg::new(0x516);
    for n in [1_usize, 2, 4, 9] {
        let a = random_symmetric(&mut rng, n, 3.0);
        let eigen = symmetric_eigen(&a, n);
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

/// The same input gives the same bits. The project's determinism rule, applied
/// to the one routine here with an iteration count.
#[test]
fn the_decomposition_is_reproducible_bit_for_bit() {
    let mut rng = Lcg::new(0xDE7);
    let a = random_symmetric(&mut rng, 11, 2.0);
    let first = symmetric_eigen(&a, 11);
    let second = symmetric_eigen(&a, 11);
    assert_eq!(first.sweeps(), second.sweeps());
    for (x, y) in first.values().iter().zip(second.values()) {
        assert_eq!(x.to_bits(), y.to_bits());
    }
    for (x, y) in first.vectors().iter().zip(second.vectors()) {
        assert_eq!(x.to_bits(), y.to_bits());
    }
}

#[test]
#[should_panic(expected = "not symmetric")]
fn an_asymmetric_matrix_is_rejected_rather_than_quietly_symmetrised() {
    // Reading one triangle, as LAPACK does, would return the decomposition of a
    // matrix the caller did not pass.
    let _ = symmetric_eigen(&[1.0, 2.0, 3.0, 4.0], 2);
}

#[test]
#[should_panic(expected = "entries")]
fn a_wrongly_sized_matrix_is_rejected() {
    let _ = symmetric_eigen(&[1.0, 2.0, 3.0], 2);
}

/// An empty problem has an empty answer rather than a panic.
#[test]
fn the_empty_matrix_decomposes_into_nothing() {
    let eigen = symmetric_eigen(&[], 0);
    assert_eq!(eigen.n(), 0);
    assert!(eigen.values().is_empty());
    assert!(eigen.converged());
    assert_eq!(eigen.sweeps(), 0);
}
