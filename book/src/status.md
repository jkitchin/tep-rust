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

**There are no recorder sinks beyond the CLI's CSV writer.** The `Recorder`
trait the plan calls for is B-0055.

**There is one integrator, and it is the original's.** Fixed-step explicit
Euler at one second, as `temain_mod.f` uses. RK4 and Dormand-Prince are B-0053.
The three-phase split of the right-hand side described in [The right-hand
side](process/rhs.md) exists precisely so that a multi-stage integrator can be
correct, and nothing has yet demonstrated one.

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
