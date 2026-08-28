# Status

This port is not finished, and the useful thing a status page can do is say
exactly where the edge is. What follows is the state of the tree at the time
this chapter was written, taken from the current-state block of `BACKLOG.org`
and from the closing entries of `LOG.org`.

## What works

**The plant model is complete and validated.** `tepsim-core` ports the whole of
`TEFUNC` (`teprob.f:196-816`), from unpacking the fifty states through the
vapour-liquid equilibrium, the kinetics, the stream table, the flows and the
compressor, the stripper, the heat transfer, the measurements and the fifty
balances. It closed with B-0025 for the balances and B-0024b for the noise and
the analysers.

**The control layer is complete.** `tepsim-control` ports the twenty control
loops of `temain_mod.f`, plus the driver's own scheduling (see
`temain_mod.f:246-317` for the velocity-form PI loop and `temain_mod.f:365-412`
for the driver). Nineteen of the twenty are ones the driver calls; the
twentieth, `CONTRL22`, is defined, tuned and never called, which is delta D-008.
Phase 4 closed with B-0041 and its full 48-hour differential.

**The public API exists.** `tepsim` gives you `Scenario`, `Simulation` and
`Run`, and `tep` drives them from a terminal. B-0052 established that
`Simulation` reproduces the previously validated run loop bit for bit, so every
claim Tiers 4 and 5 made now attaches to the public surface.

**The differential harness runs.** `tepsim-oracle` compiles and links the
unmodified Fortran and carries the Tier 1, 2, 3, 4 and 5 machinery.
`tepsim-stats` implements every statistic Tier 5 needs in Rust, with
known-answer tests, so no part of the validation depends on numpy or scipy.

The baseline recorded in `BACKLOG.org` at the close of B-0052 is 387 tests by
default and 577 with the `oracle` feature, `cargo xtask validate --tiers
1,2,3,4` green, and provenance at 1123 of 1589 claimed lines of `teprob.f` and
183 of 1408 of `temain_mod.f`, 43.6% combined.

## What does not exist yet

**There is no browser application.** Phase 8 has not started. It additionally
needs `trunk`, which is not installed on the development machine.

**There are no Python bindings.** `crates/tepsim-py` is eleven lines of doc
comment. Phase 7 was never blocked on PyO3; it was blocked on there being an
API to bind, which B-0052 has now supplied.

**There is no scenario engine.** `crates/tepsim-scenario` is sixteen lines. A
`Scenario` today is the struct described in [Getting started](getting-started.md),
which covers seed, duration, step, sampling cadence, the twenty disturbance
flags and the quirk switches, and nothing more. Arbitrary user-defined faults,
scheduled onsets and scenario files are backlog item B-0054.

**Recorder sinks exist** as of B-0055: `Columnar` for analysis, `Csv` for
reading by eye, `Ring` for a bounded live display, and `Decimating` and
`Selecting` as composable wrappers. `Simulation::run_into` streams into one
rather than collecting, which matters because a 48-hour run sampled at every
step is 172,800 rows of 53 channels. Apache Arrow and Parquet are not
implemented; whether their dependency cost is acceptable is still an open
question.

**Three integrators exist** as of B-0053: fixed-step explicit Euler, classical
RK4, and Dormand-Prince 5(4) with an embedded error estimate. **Only Euler
reproduces the original**, and it is the default; the validation ladder's every
claim is a claim about Euler. The step size is not yet adapted, because a
variable step changes when the discrete phases run and that is a decision about
fidelity rather than about numerics.

The three-phase split of the right-hand side described in [The right-hand
side](process/rhs.md) exists precisely so a multi-stage method can be correct,
and B-0053 is where that was finally demonstrated: all three methods advance
the disturbance walks once per outer step and update the gas analyser the same
number of times, which a naive four-stage method built on `TEFUNC` would not.

That comparison also produced a result about the original rather than about the
port. RK4 and Dormand-Prince agree with each other to 1.5e-6 while both differ
from Euler by about 1.1e-2, so the published Tennessee Eastman data carries
roughly one percent of integration error against an accurate solution of the
same equations. Reproducing that is the point, but it means "the TEP" names a
particular discretisation and not only a set of differential equations.

**Tier 5 is partial and Tiers 6 through 10 have not run.** The Tier 5 battery
exists and passes, but the full run covers **3 of the 21 scenarios**, stopped by
direction after about 25 minutes to redirect effort to the user-facing surface
(LOG entry B-0047b). Tier 6, the downstream detector experiment, is B-0050.
Tier 7, reproducing the published `d00` through `d21` files, is B-0051. Tiers 8,
9 and 10 have not been started.

**The validation report is written by hand.** `PLAN.org` calls for `cargo xtask
validate` to generate a chapter of this book with tables, histograms and plots.
It does not do that yet: `xtask` writes nothing into `book/`. The numbers in
[Validation](validation.md) were transcribed from `LOG.org`, and each one names
the iteration that measured it so the transcription can be checked.

**The delta register is written by hand too.** The `@delta` annotation
convention exists in the source from the first entry, but the generator that
would collect them into [the register](deltas.md) is a Phase 9 item.

## Two items awaiting sign-off

Two Class C quirks are implemented, measured, and deliberately **not** the
default, because `CLAUDE.md` requires an explicit sign-off before a
behaviour-defining fix ships. B-0025b covers the shutdown freeze
(`teprob.f:807-811`, delta D-007) and B-0040a covers the driver's hard-coded
`IDV(12)` (`temain_mod.f:366-368`, delta D-011). Neither blocks other work. The
faithful behaviour is what a default `Scenario` gives you.
