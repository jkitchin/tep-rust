# Validation

This is the part of the project that determines whether anyone trusts the
result, so it gets the most care. The strategy is a ladder: prove the pieces
exactly, prove the derivative to near machine precision, prove the stochastic
call order exactly, then prove the long-run behaviour statistically and, finally,
prove it on the downstream task people actually care about.

**Every number on this page was measured, and each one names the `LOG.org`
iteration that measured it.** None of them is a target, a plan, or a rounded
recollection. The project's operating rule is to record numbers rather than
verdicts, precisely so that a degradation inside tolerance is still visible: if
a maximum relative error moves from 3e-14 to 8e-13, both pass a 1e-12 gate and
something has broken, and the only place that is visible is the logged history.

All numbers below were produced with gfortran 15.2.0 and the pinned rustc
1.97.1. The oracle's compiler flags are fixed in `build.rs` and asserted by a
test; changing them invalidates every Tier 1 and Tier 2 number here, so it is a
logged re-baseline rather than a casual edit.

## The oracle harness

`tepsim-oracle` is a development-only crate whose `build.rs` compiles the
unmodified `teprob.f` and `temain_mod.f` with gfortran and links them through a
small C shim exposing `TEINIT`, `TEFUNC`, each `TESUBn`, and read and write
access to every `COMMON` block. A Rust test can therefore set the Fortran into
an arbitrary state, call it, and compare against the Rust implementation in the
same process.

That crate is never a dependency of `tepsim`, `tepsim-py` or `tepsim-wasm`. It
runs on Linux and macOS runners where gfortran is available, and building it
first is what converts the whole port from "read carefully and hope" into a
differential testing exercise.

Alongside it sits a cheaper check that needs no Fortran at all. `cargo xtask
fidelity` runs the port forward 100 steps from the nominal state and diffs
states, derivatives, measurements and the generator word against a golden trace
committed to the repository. It takes about a second and runs at the top of
every session. Since B-0026 it has reported **100 of 100 steps diffed, worst
3.521e-14 at `YP(12)` step 77** against a 1e-12 gate, and that same number has
now been recorded unchanged for eleven consecutive iterations, most recently in
the entry for B-0052.

## Tier 1: the utility routines, exactly

`TESUB1` (enthalpy), `TESUB2` (temperature from enthalpy by Newton), `TESUB3`
(heat capacity) and `TESUB4` (liquid density) are swept over a simplex grid, ten
million random Dirichlet samples, and a boundary pool, at every temperature in
the physical range, for each of the three `ITY` modes. The gate `PLAN.org` sets
is a maximum relative error below 1e-13 with a ULP histogram reported rather
than a pass or fail.

The measured result is not "inside 1e-13". It is zero.

| routine | `ITY` | cases | max relative error | max ULP | histogram | from |
|---|---|---|---|---|---|---|
| `TESUB1` | 0 | 9,987,490 | 0.000e0 | 0 | 0:9987490 | B-0009 |
| `TESUB1` | 1 | 9,987,490 | 0.000e0 | 0 | 0:9987490 | B-0009 |
| `TESUB1` | 2 | 9,987,490 | 0.000e0 | 0 | 0:9987490 | B-0009 |
| `TESUB3` | 0 | 9,987,490 | 0.000e0 | 0 | 0:9987490 | B-0009 |
| `TESUB3` | 1 | 9,987,490 | 0.000e0 | 0 | 0:9987490 | B-0009 |
| `TESUB3` | 2 | 9,987,490 | 0.000e0 | 0 | 0:9987490 | B-0009 |
| `TESUB4` | n/a | 9,987,490 | 0.000e0 | 0 | 0:9987490 | B-0010 |

That is 59,924,940 evaluations across `TESUB1` and `TESUB3` with zero differing
bits, p50 = p90 = p99 = p100 = 0 in all six sweeps. A separate test asserts bit
equality directly, because a 1e-13 threshold would let a drift to 1e-15 pass
unnoticed forever.

The 1e-13 gate is therefore not what is holding the line, and the entry for
B-0009 says why that is the right expectation: both routines are straight-line
arithmetic over constants already proved bit-identical, so once the association
and the literal precisions are right there is nothing left to differ by. Any
future Tier 1 routine that lands at 1e-15 rather than at zero has something
wrong with it that a tolerance would hide.

**The sabotage check.** Substituting a double-precision `273.15` for the widened
`f32` at `teprob.f:1411`, and changing nothing else, gives a maximum relative
error of **1.597e-5 at `grid#704 T=21.875`, 103,098,852,352 ULP** (B-0009). That
is the size of the error a single mis-transcribed literal produces, and it is why
constants in this project are transcribed and asserted rather than retyped.

**`TESUB2`** is the Newton solve, and it is measured from two starting
strategies (B-0011). A warm start, which is what every call site actually does,
converges in one step and recovers the answer exactly: round-trip error 0e0. A
cold start from the far end of the range is what exercises the iteration, and its
worst round trip is **5.12e-13 C**, at `dirichlet#720561 T=171.86`,
`face#1309309 T=146.56` and `face#96 T=144.14`. Over 59,924,940 solves there
were zero differing bits and **zero abandoned** iterations, which is the
measurement that makes delta D-001's effect zero on the physical domain.

**`TESUB7`**, the random generator, must be exact and is (B-0005): 10^7 draws
from the compiled-in seed with the XOR fold and the final state both matching;
10^6 draws each from five dataset seeds; draw-by-draw comparison over 200,000
draws for four seeds with every draw and every intermediate state bit-identical;
and interleaved output modes over 50,000 draws on an irregular pattern. Its
exactness rests on reproducing the rounding rather than removing it: the product
exceeds 2^53 on 0.7716 of draws, against 0.7728 predicted from the arithmetic,
so a "fixed" integer recurrence diverges at draw 0.

**`TESUB5` and `TESUB6`** are exact in both `libm` configurations, since neither
touches a transcendental, and their draw counts were recovered independently by
stepping a port-side generator from the word before each call to the word after:
**3 draws for `TESUB5`, 12 for `TESUB6`**, on every case and for both flag
values (B-0028).

## Tier 2: single-step derivative equivalence

Both implementations are forced into an identical state, all fifty states, the
twelve manipulated variables, the twenty `IDV` flags, the full walk state, the
generator word and the nineteen held analyser readings, evaluated once, and
compared on all fifty components. Sampling is from three pools: states along the
nominal closed-loop trajectory, random perturbations of those states scaled
across several orders of magnitude, and adversarial states placed deliberately
at every discontinuity and clamp in the model.

**The tolerance is relative to the scale of the terms, not to the result**, and
that is a decision with a measurement behind it rather than a relaxation. A
balance is inflow minus outflow, and near steady state those nearly agree.
`YP(2)`, the inert's reactor balance, is a difference of two flows around 660
whose result is a few parts in ten thousand of either. One ULP of difference in
each term, which is all the vendored `libm` costs, is 1e-16 of the terms and
1e-4 of the result. Measured against the result, 28 of the 50 components exceed
1e-12 while the whole right-hand side is *bit-identical* to gfortran under
`libm-system`. Each balance therefore reports the magnitude of its largest term
alongside its value, and the gate is the error over that.

Acceptance, from B-0026: all fifty components, **2,412 running states**, all
three pools, **worst 6.093e-14 at `YP(7)`, `perturbed#300`**, against a 1e-12
gate.

| `YP` | error / scale | error / value | ratio |
|---|---|---|---|
| 2 | 3.811e-14 | 1.393e-4 | 3.7e9 |
| 30 | 4.552e-15 | 2.487e-6 | 5.5e8 |
| 38 | 8.696e-15 | 1.187e-6 | 1.4e8 |
| 14 | 3.113e-14 | 2.618e-6 | 8.4e7 |
| 7 | 6.093e-14 | 7.835e-7 | 1.3e7 |
| 19-27, 37, 39-50 | 0.000e0 | 0.000e0 | 1x |

**28 of 50 components cancel by more than 100x**, and the 22 that do not are
exactly 0.000e0: bit-identical rather than merely inside the gate. The
acceptance test asserts that contrast and fails if the count of heavily
cancelling components drops below twenty, so the reasoning behind the decision
expires loudly if the model ever stops cancelling.

## Tier 3: RNG call-order equivalence

Both sides are instrumented to emit every draw, and the traces are diffed. This
test exists specifically because it is the one the existing Python port would
fail: its documented divergence is attributed to "the exact sequence of calls
differs due to implementation details", which is exactly the defect a trace diff
catches on the first run instead of after a 48-hour statistical comparison.

From B-0029, draw counts at four points in a run, with **exact agreement at
every one**:

| time | draws | what is drawing |
|---|---|---|
| 0 | 0 | noise skipped, walks reset |
| 1e-6 | 264 | noise only |
| 0.15 | 462 | noise, walk advance, gas analysers |
| 0.30 | 522 | and the product analyser |

Scalings in one real evaluation: 30 signed and 432 unit. The 30 is 9 x 3 + 3,
the walk advance; the 432 is 36 compositions at twelve draws each. Worst
evaluation is 462 draws against a trace buffer capacity of 4096. The instrument
is checked against the generator word rather than only against the values, which
covers completeness, ordering and fidelity at once, and it holds over 113,088
draws.

Two independent methods agree on the counts. B-0027 measured them with no
instrumentation at all, by stepping a port-side generator from the word before a
call to the word after, and the trace lengths match that census exactly. Either
method alone would be a claim; the two agreeing is evidence.

## Tier 4: trajectory equivalence, diagnostic

Tier 4 is a diagnostic, not a gate. Long-horizon divergence between two programs
that use different `exp` implementations is expected. `PLAN.org` asks for two
things: that the error stay below the corresponding measurement noise standard
deviation `XNS(i)` for at least the first several hours, and that the onset of
divergence be *explained* by showing that switching `libm` moves it.

**Open loop, nominal, 8 simulated hours (28,800 steps)**, from B-0034:

| `libm` | worst error, as a fraction of `XNS(i)` | ever outside `XNS` |
|---|---|---|
| vendored | 5.149e-5 | never |
| platform | 0.000e0 | never |

All 21 scenarios at 4 hours each stay within `XNS` for the whole run in both
configurations, and every one is exactly 0.00e0 on the platform `libm`. Worst
error by scenario on the vendored `libm` at 4 hours: nominal 2.45e-5, `IDV(10)`
3.66e-8, `IDV(8)` 1.42e-8, `IDV(13)` 1.05e-8, and the rest below 1e-8.

The explanation is better than the one asked for. With identical transcendentals
the two trajectories are not merely close, they are **bit-identical for 28,800
steps on all 41 measurements**. So the divergence is transcendental rounding and
nothing else, and the residual under the vendored `libm` is twenty thousand
times below what the instruments could resolve after eight hours.

**Closed loop, 48 hours, 172,800 steps**, nominal scenario with the driver's
forced `IDV(12)` from hour eight, from B-0041:

| `libm` | worst `XMEAS` | worst `XMV` (fraction of range) | within `XNS` | at the end |
|---|---|---|---|---|
| platform | 0 | 0 | 48.000 h | 0 |
| vendored | 1.705e-10 at `XMEAS(1)` | 1.246e-11 at `XMV(10)` | 48.000 h | 7.882e-11 x `XNS(13)` |

Bit-identical for all 172,800 closed-loop steps under the platform `libm`. Not
one measurement, not one valve, not one step. Under the vendored `libm` the
error never reaches a tenth of a billionth of any instrument's noise over two
simulated days.

Compare the open-loop figure above: nominal at 8 hours ends at 5.149e-5 of
`XNS`, six orders of magnitude worse. Closing the loop suppresses the
amplification, which is what a controller pulling toward a setpoint should do
and is worth having measured rather than assumed.

The same run measures what the control layer is actually buying. Open loop from
the same start, with the valves held, **trips at step 11,017 (3.060 h) on
reactor pressure high**. So "48 hours without tripping" is a statement about the
control layer, not about a placid plant. B-0052 reproduced that trip from the
public API and the CLI, at exactly the same step, without either being told
about the other.

## Tier 5: statistical equivalence, the real gate

Tier 5 tests *equivalence* rather than difference, because a failure to reject
the null of no difference is not evidence of equivalence and two one-sided tests
are. Every statistic ships in Rust, in `tepsim-stats`, with known-answer tests:
Welch's t and the TOST wrapper, the two-sample Kolmogorov-Smirnov statistic,
energy distance, autocorrelation, Welch power spectra and the Pearson
correlation matrix. Nothing calls out to numpy or scipy.

### The harness

From B-0047a: a 48-hour run produces **960 samples of 53 variables**, a port run
takes **766 ms** in release, and a full battery of 2100 runs takes about **27
minutes per source**. Over three scenarios the two sources agree **bit-identically
(0)** under the platform `libm` and to **2.056e-13** under the vendored one, while
five nominal seeds spread by 1.366e1, which is the scale that makes those first
two numbers meaningful.

The same entry checked that the scenarios are not vacuous. All twenty
disturbances move the plant within two hours; the smallest departures are
`IDV(3)` at 4.076e-3, `IDV(9)` at 3.531e-3 and `IDV(15)` at 9.839e-3, and the
largest are `IDV(6)` at 3.467e0 and `IDV(1)` at 2.456e0.

### The battery

Four of the six statistics have margins that are *measured* rather than chosen
(B-0047b). The reference's seeds are split in half, the statistic is computed
against itself over twenty deterministic splits, and the cross-source value is
treated as one more draw from that null, giving `p = (1 + #{within >= cross}) /
(K + 1)` gated at 0.05. With twenty splits the cross-source value fails exactly
when it is the strict maximum of the twenty-one.

The smoke battery in the CI gate is 3 scenarios, 4 seeds, 2 hours, 12 runs per
source, 14 seconds. The full battery is **partial: 3 of 21 scenarios** at 100
seeds by 48 hours, stopped by direction after about 25 minutes.

| scenario | worst mean power | Frobenius, cross | within, max | p |
|---|---|---|---|---|
| nominal | 1.59 | 1.186e-2 | 7.595 | 1.0000 |
| `IDV(1)` | 1.61 | 7.080e-3 | 5.703 | 1.0000 |
| `IDV(2)` | 0.23 | 2.827e-5 | 1.924 | 1.0000 |

Every calibrated statistic passed at p = 1.0000 on all three: the cross-source
value was never the maximum of its null. The correlation matrix's Frobenius
distance is three to five orders of magnitude inside the within-source spread,
which matters because that matrix is exactly what a PCA-based detector consumes.

TOST power against battery size, measured in the same entry:

| battery | worst power |
|---|---|
| 4 seeds, 2 h (smoke) | 11.93 |
| 8 seeds, 12 h | 3.38 |
| 8 seeds, 48 h | 1.04 |
| 100 seeds, 48 h, nominal | 1.59 |
| 100 seeds, 48 h, `IDV(2)` | 0.23 |

A power above 1 means the battery is **underpowered for the margin**, not that
the sources differ. `PLAN.org` sets the mean margin at a tenth of the pooled
standard deviation, and a disturbed plant has a much larger pooled spread than a
quiet one, so the same absolute run-to-run wander sits comfortably inside the
margin under `IDV(2)` and outside it at the nominal operating point. At the
nominal plant the margin needs about **255 seeds**, not 100. The variables that
run out of power are the manipulated ones: an integrating controller's output
has a random-walk component, so its mean over 48 hours varies between seeds by
much more than a tenth of its own within-run spread. That is a fact about the
plant, not about the port, and the battery reports it as undecided rather than
as a difference while still asserting on the measured gap.

### Physics invariants

These are the only tests in the whole ladder that can catch an error the port
*faithfully inherited*, because every other tier compares against the Fortran
and an error in the Fortran is invisible to all of them. From B-0046:

| invariant | Fortran | port | gate |
|---|---|---|---|
| I-1 reaction mass, 200 states | 4.664e-16 | 3.498e-16 | 1e-14 |
| I-2 plant mass balance, nominal | 2.079e-16 | 2.344e-16 | 1e-13 |
| I-2 along 2 h open loop | 6.556e-16 | not run | 1e-13 |
| I-3 inert reaction term | exactly 0 | exactly 0 | exact |
| I-4 per-component moles, 1 h | 7.436e-15 | 2.892e-15 | 1e-13 |

All four reactions balance exactly with the published molecular weights, `2 + 28
+ 32 = 62`, `2 + 28 + 46 = 76`, `2 + 46 = 48` and `3 x 32 = 2 x 48`, which is
why I-1 is an equality rather than a tolerance.

The invariants' teeth were checked by mutation, and one result is worth
repeating because it is not obvious. Deleting the reaction term from the
reactor's component balance leaves the **total mass balance passing at 2e-16**.
That is not a weak test, it is a true fact about the invariant: I-2 is the
molecular-weight-weighted sum of I-4, and I-1 says the reaction is mass-neutral,
so the reaction term cancels out of I-2 exactly. Total mass conservation is
blind to stoichiometry by construction. The unweighted per-component balance,
I-4, sees it immediately. An invariant that is a sum of other invariants is
weaker than the set, and how much weaker is not obvious from reading it.

## Tiers 6 through 10

These have not run. They are listed here so that the shape of the remaining
claim is visible rather than implied.

**Tier 6, downstream-task equivalence**, is the operational definition of
"practically equivalent" for this project's research audience: train a detector
suite on Fortran `d00` data, evaluate on Fortran and on Rust test data, then
reverse it, and compare detection rate, false alarm rate and detection delay.
The claim to be able to make is that cross-source performance matches
within-source performance within its own run-to-run variability. Backlog item
B-0050.

**Tier 7, published dataset reproduction**, attempts direct reproduction of the
bundled `d00` through `d21` files under the documented generation protocol,
reporting per-file agreement with the Tier 5 machinery. The groundwork exists:
`teprob.f:1187-1256` carries fifty-four generator words in comments, one per
published dataset, and B-0047a transcribed and asserted them. Three facts about
that table are worth knowing in advance. Thirty-four of the fifty-four exceed
2^32, which is not a transcription error because `TESUB7` reduces any seed on
the first draw. Twenty-seven are even, and a multiplicative generator modulo a
power of two keeps the factors of two its seed has, so those runs have a shorter
period and low bits that never move. Reproducing the published files means doing
the same rather than fixing it. Backlog item B-0051.

**Tier 8** is differential fuzzing, **Tier 9** is cross-platform determinism by
golden BLAKE3 digest including wasm in a real browser, and **Tier 10** requires
that every quirk fix ship with a measured delta from the full Tier 5 battery
with the fix on and off. None has started.

## Two bugs that only long runs found

Worth recording because they are the argument for running the battery long
rather than often (B-0047b).

`TRCN`, the Tier 3 trace counter, is a Fortran `INTEGER` that nothing clears
between evaluations. At about 264 draws per step over 172,800 steps it passes
2^31 after roughly fifty runs, goes negative, passes the capacity guard, and
writes outside the array. The symptom was a SIGSEGV deep inside the Fortran,
tens of millions of steps into the battery, with nothing nearby to suggest why:
the same seed ran perfectly in isolation, and starting at seed 40 moved the
crash to seed 90, which is what identified it as cumulative rather than
data-dependent. Nothing in Tiers 1 to 4 was ever affected, because no shorter run
approached the overflow.

`welch_t` computed `(n - 1)` on a summary with no observations. In debug that
panics; in release it wraps and returns a degrees-of-freedom figure that looks
like a number. The case is `XMV(12)`, the agitator, which no controller ever
writes, so every run has zero variance and the log-variance ensemble is empty.

Both are the same lesson. A validation harness is code, it has bugs, and the
bugs it has are the ones only its longest runs reach.
