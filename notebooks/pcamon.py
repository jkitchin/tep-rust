"""PCA process monitoring, in NumPy and the Python standard library.

This is the detector the notebooks in this directory use: fit a principal
component model to fault-free data, then watch two statistics, Hotelling's
T-squared inside the retained subspace and the squared prediction error (also
written Q, or SPE) outside it. It is the scheme described by Chiang, Russell
and Braatz, and the one whose Tennessee Eastman results are quoted throughout
the fault detection literature.

Nothing here needs SciPy or scikit-learn. The eigendecomposition is
``numpy.linalg.eigh``, the normal quantile is ``statistics.NormalDist``, and
the F quantile is a continued-fraction incomplete beta with a bisection on top,
which is about forty lines and removes the last reason to add a dependency.

Every function is a deliberate match for the corresponding item in
``tepsim-stats``, the Rust crate the repository's own validation ladder runs
on, so the numbers a notebook prints can be checked against the transcript in
``book/src/tutorials/a-detector.md``. Where the two could differ, the
docstrings say which convention was chosen.
"""

import math
from statistics import NormalDist

import numpy as np

__all__ = [
    "fit",
    "statistics",
    "t_squared_limit",
    "spe_limit",
    "alarms_above",
    "detection_delay",
    "report",
    "f_quantile",
    "regularised_beta",
]


# ---------------------------------------------------------------------------
# the model
# ---------------------------------------------------------------------------


def fit(X, retain=0.90):
    """Fit a PCA model to a fault-free training matrix.

    ``X`` is samples by variables. Each column is centred on its mean and
    divided by its **sample** standard deviation (the ``n - 1`` one), so the
    matrix that gets diagonalised is the correlation matrix and every variable
    contributes equally regardless of its engineering units. That matters here:
    reactor pressure is near 2705 kPa and a composition is a percentage, and a
    covariance model would be a model of the pressure alone.

    A column that never moved in training has no scale to divide by. It is set
    to exactly zero, its row and column of the correlation matrix are zeroed,
    and its index is reported in ``constant``. The Tennessee Eastman plant has
    such a column in ordinary closed-loop operation: ``XMV(12)``, the agitator
    speed, is held fixed by the base case control scheme.

    ``retain`` is the fraction of total variance the kept components must
    explain. It is named rather than a bare component count because two
    detectors that keep different numbers of components are different
    detectors, and a reported result that does not say which rule produced it
    cannot be reproduced.

    Returns a dict, which is enough structure for the job and keeps the model
    printable.
    """
    X = np.asarray(X, dtype=float)
    n, p = X.shape
    if n < 2:
        raise ValueError("a sample standard deviation needs at least two samples")

    mean = X.mean(axis=0)
    centred = X - mean
    sd = centred.std(axis=0, ddof=1)
    constant = np.flatnonzero(sd <= 0.0)

    safe = np.where(sd > 0.0, sd, 1.0)
    Z = np.where(sd > 0.0, centred / safe, 0.0)
    R = (Z.T @ Z) / (n - 1)
    R[constant, :] = 0.0
    R[:, constant] = 0.0

    values, vectors = np.linalg.eigh(R)
    order = np.argsort(values)[::-1]
    values = values[order]
    vectors = vectors[:, order]

    total = values.sum()
    # LAPACK's rank convention, not a tuned threshold. Components past this have
    # eigenvalues indistinguishable from zero, and T-squared divides by the
    # eigenvalue, so retaining one would turn rounding noise into an alarm.
    rank = int(np.count_nonzero(values > p * np.finfo(float).eps * values[0]))

    running, k = 0.0, 0
    while k < len(values) and running / total < retain:
        running += values[k]
        k += 1
    k = min(k, rank)

    return {
        "mean": mean,
        "sd": sd,
        "constant": constant,
        "values": values,
        "vectors": vectors,
        "total": total,
        "rank": rank,
        "k": k,
        "samples": n,
        "variables": p,
        "retain": retain,
        "explained": float(values[:k].sum() / total),
    }


def statistics(model, X):
    """Both monitoring statistics for every row of ``X``.

    T-squared is the squared Mahalanobis distance from the training mean,
    measured only in the subspace the model retained::

        T^2 = sum_{j < k} (p_j' z)^2 / lambda_j

    Each score is divided by its own eigenvalue, so a direction that barely
    moved during training weighs a movement far more heavily than one that
    swung widely. That division is what makes it a distance and not just a sum
    of squares.

    SPE is what is left over after the retained subspace has explained all it
    can::

        e = z - P_k P_k' z,   SPE = e'e

    The residual is formed by subtracting the reconstruction rather than by
    summing the discarded scores. The two are equal in exact arithmetic, and
    the first is the one that stays right when the loadings are only
    orthonormal to within rounding, because it measures what the model actually
    failed to reconstruct.

    Returns ``(t_squared, spe)``, each a one-dimensional array with one entry
    per row.
    """
    X = np.asarray(X, dtype=float)
    sd = model["sd"]
    safe = np.where(sd > 0.0, sd, 1.0)
    Z = np.where(sd > 0.0, (X - model["mean"]) / safe, 0.0)

    P = model["vectors"][:, : model["k"]]
    scores = Z @ P
    t2 = np.sum(scores**2 / model["values"][: model["k"]], axis=1)
    residual = Z - scores @ P.T
    spe = np.sum(residual**2, axis=1)
    return t2, spe


# ---------------------------------------------------------------------------
# control limits
# ---------------------------------------------------------------------------


def t_squared_limit(components, samples, confidence):
    """The upper control limit for T-squared on a new observation::

        T2_alpha = k (n + 1)(n - 1) / (n (n - k)) * F_alpha(k, n - k)

    There are two limits in circulation and they are not interchangeable. This
    is the one for an observation independent of the training set, which is the
    case a fault detection experiment is in: the model is fitted to fault-free
    data and then applied to a different, faulted record. The other, for an
    observation that was itself part of the training set, uses a Beta quantile
    and belongs to a screen of the training data.
    """
    if components == 0:
        return 0.0
    if samples <= components or not 0.0 < confidence < 1.0:
        return float("nan")
    k, n = float(components), float(samples)
    return k * (n + 1.0) * (n - 1.0) / (n * (n - k)) * f_quantile(confidence, k, n - k)


def spe_limit(residual_eigenvalues, confidence):
    """The upper control limit for SPE, by the Jackson-Mudholkar approximation::

        theta_i   = sum_j lambda_j^i                    i = 1, 2, 3
        h0        = 1 - 2 theta1 theta3 / (3 theta2^2)
        SPE_alpha = theta1 [ c sqrt(2 theta2 h0^2) / theta1
                             + 1
                             + theta2 h0 (h0 - 1) / theta1^2 ] ^ (1 / h0)

    with ``c`` the standard normal deviate at the confidence level.
    ``residual_eigenvalues`` is the discarded tail of the spectrum.

    SPE is a weighted sum of squared normals whose exact distribution has no
    closed form, so the three moments are matched with a power transformation
    instead, which is why only the first three power sums appear. Eigenvalues
    slightly below zero are rounding noise on a rank deficient matrix and are
    clamped, because a negative cube would bias ``h0``.
    """
    lam = np.clip(np.asarray(residual_eigenvalues, dtype=float), 0.0, None)
    theta1, theta2, theta3 = lam.sum(), (lam**2).sum(), (lam**3).sum()
    if theta1 <= 0.0:
        return 0.0
    if not 0.0 < confidence < 1.0:
        return float("nan")
    h0 = 1.0 - 2.0 * theta1 * theta3 / (3.0 * theta2**2)
    c = NormalDist().inv_cdf(confidence)
    bracket = (
        c * math.sqrt(2.0 * theta2 * h0 * h0) / theta1
        + 1.0
        + theta2 * h0 * (h0 - 1.0) / theta1**2
    )
    return float(theta1 * bracket ** (1.0 / h0))


# ---------------------------------------------------------------------------
# alarms and what they are worth
# ---------------------------------------------------------------------------


def alarms_above(statistic, limit):
    """Threshold a statistic into a boolean alarm series, strictly greater.

    A statistic landing exactly on its limit is not an exceedance. The case is
    unreachable on real data and reachable on constructed data, and choosing
    the side deliberately is cheaper than discovering later which side was
    chosen.
    """
    return np.asarray(statistic, dtype=float) > limit


def detection_delay(alarms, onset, consecutive=3):
    """Samples from the onset to the first run of ``consecutive`` alarms.

    ``None`` if no such run fits inside the record. The persistence
    requirement is what stops the delay being a measurement of luck: at a two
    per cent false alarm rate the first post-onset sample alarms by chance one
    time in fifty, and calling that a delay of zero flatters the detector. The
    run length is a parameter because the literature uses three and six and
    does not agree, so a reported delay has to carry the value it was measured
    with.
    """
    if consecutive < 1:
        raise ValueError("a detection needs at least one alarm")
    post = np.asarray(alarms, dtype=bool)[onset:]
    for start in range(len(post) - consecutive + 1):
        if post[start : start + consecutive].all():
            return start
    return None


def report(alarms, onset, consecutive=3):
    """Detection rate, false alarm rate, delay, and the counts behind them.

    ``onset`` is the index of the first sample **with** the fault present, so
    that sample counts as post-fault. The detection rate is a fraction of
    post-onset samples, not a per-run yes or no: a detector that catches a
    fault and then loses it scores badly, which is the intent, because a
    statistic that drops back inside its limit while the fault is still running
    is flickering rather than detecting. Its complement, the missed detection
    rate, is the quantity the Tennessee Eastman tables in the literature
    report.

    The counts travel with the rates because a rate without its denominator
    cannot be compared against the next run.
    """
    alarms = np.asarray(alarms, dtype=bool)
    boundary = min(onset, len(alarms))
    pre, post = alarms[:boundary], alarms[boundary:]
    return {
        "onset": onset,
        "samples": len(alarms),
        "pre_fault": len(pre),
        "post_fault": len(post),
        "false_alarms": int(pre.sum()),
        "detections": int(post.sum()),
        "fault_detection_rate": float(post.mean()) if len(post) else float("nan"),
        "missed_detection_rate": 1.0 - float(post.mean()) if len(post) else float("nan"),
        "false_alarm_rate": float(pre.mean()) if len(pre) else float("nan"),
        "detection_delay": detection_delay(alarms, onset, consecutive),
        "consecutive": consecutive,
    }


# ---------------------------------------------------------------------------
# quantiles, so that neither SciPy nor a lookup table is needed
# ---------------------------------------------------------------------------


def _beta_continued_fraction(a, b, x):
    """Lentz's method on the continued fraction for the incomplete beta."""
    max_iterations, eps, tiny = 300, 3e-16, 1e-300
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c = 1.0
    d = 1.0 - qab * x / qap
    if abs(d) < tiny:
        d = tiny
    d = 1.0 / d
    h = d
    for m in range(1, max_iterations + 1):
        m2 = 2 * m
        numerator = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + numerator * d
        if abs(d) < tiny:
            d = tiny
        c = 1.0 + numerator / c
        if abs(c) < tiny:
            c = tiny
        d = 1.0 / d
        h *= d * c
        numerator = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + numerator * d
        if abs(d) < tiny:
            d = tiny
        c = 1.0 + numerator / c
        if abs(c) < tiny:
            c = tiny
        d = 1.0 / d
        step = d * c
        h *= step
        if abs(step - 1.0) < eps:
            break
    return h


def regularised_beta(a, b, x):
    """The regularised incomplete beta function I_x(a, b)."""
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    log_beta = math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
    if x < (a + 1.0) / (a + b + 2.0):
        front = math.exp(log_beta + a * math.log(x) + b * math.log1p(-x))
        return front * _beta_continued_fraction(a, b, x) / a
    front = math.exp(log_beta + b * math.log1p(-x) + a * math.log(x))
    return 1.0 - front * _beta_continued_fraction(b, a, 1.0 - x) / b


def f_quantile(p, d1, d2):
    """The p-quantile of the F distribution with ``d1`` and ``d2`` degrees of
    freedom, by bisection on its CDF.

    The CDF is ``I_{d1 x / (d1 x + d2)}(d1/2, d2/2)``, which is monotone in
    ``x``, so bisection converges to the last bit in about sixty iterations and
    needs no derivative and no starting guess.
    """

    def cdf(x):
        return regularised_beta(d1 / 2.0, d2 / 2.0, d1 * x / (d1 * x + d2))

    lo, hi = 0.0, 1.0
    while cdf(hi) < p and hi < 1e12:
        hi *= 2.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if cdf(mid) < p:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)
