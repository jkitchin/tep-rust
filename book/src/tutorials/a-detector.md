# Building a detector, and measuring it

> **The worked example is
> [notebook 2, Fault detection with PCA](../notebooks/02-fault-detection-pca.html),**
> which fits the model, plots both statistics against their limits on fault-free
> and faulted records, scores eight disturbances, and then diagnoses the one
> limit that does not work. Its sequel,
> [notebook 3, The three hard faults](../notebooks/03-hard-faults.html),
> measures the benchmark's detectability floor three separate ways. Sources:
> `notebooks/02-fault-detection-pca.ipynb` and `notebooks/03-hard-faults.ipynb`.

This is the canonical Tennessee Eastman monitoring experiment. Fit a principal
component model to a record of the fault-free plant, watch two statistics on new
data, and raise an alarm when either leaves its control limit. Hotelling's
T-squared measures distance from the training mean inside the subspace the model
retained. The squared prediction error, written SPE or Q, measures how much of
an observation the model could not reconstruct at all.

The detector is `notebooks/pcamon.py`, about four hundred lines of NumPy and the
Python standard library with no SciPy and no scikit-learn: the
eigendecomposition is `numpy.linalg.eigh`, the normal quantile is
`statistics.NormalDist`, and the F quantile the T-squared limit needs is a
continued-fraction incomplete beta with a bisection on top.

That file, rather than a crate, is deliberate. The repository has a Rust
implementation of exactly this detector in `tepsim-stats`, and this page used to
show it. `tepsim-stats` is development-only, and `cargo xtask ci` asserts that
no shipped crate so much as names it, so a reader who installed `tepsim` and
followed the page could not build what the page was teaching. Everything the
notebooks do runs against the package you can install.

The two implementations agree. Notebook 2 opens by reproducing the Rust
transcript, and every printed digit matches: the same 33 components, the same
limits to three decimals, and the same four detection rates and delays, from a
hand-written cyclic Jacobi sweep on one side and LAPACK on the other. What the
notebook measures after that is therefore a property of the plant and the
method, not of the linear algebra.

## The measurement is the point

Any detector can be made to look good by reporting only the faults it catches.
Every detection rate in the notebook is reported next to the false alarm rate on
held-out fault-free data the model never saw, and one of the two statistics
comes out badly.

Three numbers do the work, and each has a trap in it.

The **fault detection rate** is the fraction of post-onset samples that raised
an alarm. It is a rate over samples, not a per-run yes or no, so a detector that
catches a fault and then loses it scores badly. That is the intent: a statistic
that drops back inside its limit while the fault is still running is flickering,
not detecting. Its complement, the missed detection rate, is what the tables in
the literature report.

The **false alarm rate** is the same fraction over the pre-onset samples, and it
is the number that makes a detection rate meaningful. A detector with a 20%
false alarm rate that achieves a 20% detection rate has detected nothing.

The **detection delay** is the number of samples from the onset to the first run
of three consecutive alarms. The persistence requirement is not decoration. With
a run length of one the delay is just the first alarm after the onset, and on a
detector with any false alarm rate at all that is mostly luck: at a 3% false
alarm rate the first post-onset sample alarms by chance one time in thirty, and
calling that a delay of zero flatters the detector. The literature uses three
and six and does not agree, so the run length travels with the number. Delays
are in samples, and one sample is three minutes.

## What the model does with a constant column

`pcamon.fit` standardises each column to zero mean and unit sample variance
before forming the correlation matrix. Standardisation rather than mere centring
is not optional here: reactor pressure lives near 2705 kPa and a composition is
a percentage, so a covariance model would be a model of the pressure and nothing
else.

How many components to keep is passed as a named rule rather than a bare `k`,
because two detectors that retain different numbers of components are different
detectors and a result that does not say which rule produced it cannot be
reproduced. At 90% of the variance the notebook keeps 33 of 52 components. That
is two thirds of them, and it is worth pausing on: the largest eigenvalue is
5.910, only 11.4% of the total of 52, so this plant has no small handful of
dominant directions and a monitoring scheme on it works in a fairly
high-dimensional retained subspace.

`XMV(12)`, the agitator speed, never moves, so its standard deviation is zero
and it cannot be standardised. Rather than divide by zero or quietly drop the
column, `pcamon` records its index, zeroes its row and column of the correlation
matrix, and says so when asked. A silent drop here is how a 53-variable model
becomes a 52-variable model that nobody can reproduce.

## The number the detector would rather you did not see

Both limits are drawn at 99% confidence, so the nominal false alarm rate is 1%.
On a fresh 48-hour fault-free run the model never saw, T-squared alarms on 32 of
960 samples, a rate of 0.0333, which is high but recognisable. SPE alarms on 174
of 960, a rate of 0.1812. That is a factor of eighteen.

It is not an arithmetic error. `pcamon.spe_limit` computes the
Jackson-Mudholkar expression exactly as stated. It is the assumption underneath
the expression failing: the limit is derived for residuals that are normal and
independent, and the plant's are neither. Several feed conditions are driven by
slow random walks that never stop, so a 48-hour record wanders somewhere a
25-hour training record did not go, and when it does, every sample in the
excursion alarms together. The notebook's plot shows exactly that shape, with
the SPE alarms arriving in long blocks rather than as isolated points.

The tempting response is to lower the confidence level until the number looks
right. The notebook does the experiment instead, holding everything else fixed
and lengthening the training record:

| Training record | Samples | SPE false alarm rate |
|---|---|---|
| 25 hours | 500 | 0.1812 |
| 100 hours | 2000 | 0.0219 |
| 200 hours | 4000 | 0.0177 |
| 500 hours | 10000 | 0.0115 |

That settles it. With enough fault-free training data both statistics land on
the nominal 1%, so Jackson-Mudholkar was never the problem. A 25-hour record
simply does not contain the tails of a process whose feed conditions are driven
by random walks that never stop, and the limit it produces is therefore too
tight.

This has a direct consequence for the literature, and not a comfortable one. The
published training file `d00` holds 500 samples, which at a three-minute
interval is exactly 25 hours: the first row of that table. Every SPE false alarm
rate reported for static PCA on the published Tennessee Eastman data inherits
this.

Estimating the limit empirically from the training residuals is the usual
advice, and the notebook checks that too rather than assuming it. It does not
help: the empirical SPE limit gives 0.1448 against the analytic 0.1812, and the
empirical T-squared limit is markedly *worse* than the analytic one at 0.1042
against 0.0333. Both are drawn from the same 25 hours that were too short in the
first place, and neither can know about the excursions it never saw.

The fix is more data, or a method that models the serial correlation rather than
assuming it away. Dynamic PCA, which augments each observation with lagged
copies of itself, and canonical variate analysis are the two the literature
reaches for, and Russell, Chiang and Braatz's 2000 comparison of both against
static PCA is on exactly these files.

## The floor under the benchmark

Notebook 3 is about the best known empirical result on this problem: `IDV(3)`,
`IDV(9)` and `IDV(15)` are effectively undetectable by these methods. It
measures that on the published `d00` through `d21` files, on simulated runs at
the seeds those files record, and over a ten-seed ensemble, and the three
measurements agree.

The plainest form of the result is distributional. The median post-onset
T-squared for those three faults is within one percent of the fault-free median,
and their median SPE within four percent, on quantities whose fault-free values
span two orders of magnitude. No threshold separates them because there is
nothing to separate.

That is worth stating carefully, because it reads like a criticism of PCA and is
not one. The plant really is running normally under those three disturbances. A
detector that alarmed on them would be reporting a fault with no consequence,
and the correct behaviour for a monitoring scheme is what it does. What the
result measures is that this benchmark has a detectability floor set by the
process rather than by the method, which is part of what makes it a good
benchmark: any paper claiming to detect all twenty is claiming something about
its false alarm rate that it has probably not measured.
