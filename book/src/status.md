# Status

**Everything planned is built.** Every phase has landed, every validation tier
has a harness, and the four decisions that were once open have been taken. The
backlog is 85 items done, none open and none blocked.

That makes this page shorter than it used to be, and its job is now to say
where the *limits* are rather than where the edge of construction is. Numbers
here come from the current-state block of `BACKLOG.org` and the closing entries
of `LOG.org`, and each names the iteration that measured it.

## What exists

The plant model, the control layer, the public API, the command line, the
Python wheel and the browser app are all complete and validated.
`tepsim-core` ports the whole of `TEFUNC` (`teprob.f:196-816`);
`tepsim-control` ports the twenty control loops and the driver's own
scheduling; `tepsim` is the API most callers want; `tep` drives it from a
terminal, including `tep dataset` for generating `d00`-`d21` shaped files; the
wheel carries no C dependency; and [TEP Studio](studio/) runs the whole
simulator in a browser tab with no server.

The differential harness is `tepsim-oracle`, which compiles and links the
unmodified Fortran, and `tepsim-stats`, which implements every statistic the
ladder needs in Rust with known-answer tests, so no part of the validation
depends on numpy or scipy.

**All ten tiers have a harness and all ten have run.** Tier 8, differential
fuzzing with shrinking, and Tier 9, cross-platform determinism, were the last
two and landed together. The [validation chapters](validation/index.md) are
generated from the suite's own output for the tiers that have a generator; the
narrative in [Validation](validation.md) is hand-written and transcribed from
`LOG.org`, and each number there names the iteration that measured it so the
transcription can be checked.

## What is genuinely still missing

These are limits, not unfinished work, and most of them are limits of the
machine this was built on rather than of the code.

**Tier 9 has no x86-64 leg and no real-browser leg.** Six committed digests are
identical on aarch64 and on wasm32 under Node, across three build profiles, and
the browser app's transport path reproduces them independently. Nobody has run
them on x86-64, Windows or aarch64 Linux, and Node is not a browser. The table
is committed constants rather than a value computed twice in one process,
precisely so that running `cargo xtask tier9` on another machine completes the
claim with no code change.

**Tier 8 has one open counterexample.** In five million tuples, `fuzz#863105`
misses the 1e-12 gate at 4.607e-12 of the scale of its terms. It is recorded
and attributed rather than fixed, because there is nothing to fix: it is an
accepted one-ULP `exp` difference amplified by 973 simulated hours of
`IDV(13)`'s kinetic drift, and it is bit-identical under `libm-system`.

**Two `IDV` faults miss the Tier 5 margin under the vendored libm.** `IDV(14)`
and `IDV(19)`, on exactly the valves those faults stick, and bit-identical
under the platform libm. The cause is `teprob.f:801`, a discontinuous branch on
a floating-point comparison. Those valves are judged on their distribution
rather than on a mean, because a series of plateaux has no meaningful centre.
No margin was widened.

**Apache Arrow and Parquet sinks are not implemented.** The recorder sinks that
exist are `Columnar`, `Csv`, `Ring`, `Decimating` and `Selecting`, and
`Simulation::run_into` streams into one rather than collecting. Arrow's
dependency cost was judged not to be worth paying until someone has a dataset
large enough that CSV is the bottleneck; the right home would be `tepsim-cli`,
never `tepsim`, which is `no_std` and compiles to wasm32 under a size budget.

**A historian export does not exist.** Everything is wide tabular. A long-format
`tag,timestamp,value,quality` record would need a caller-supplied epoch, since
the core may not read a clock, a decision about what quality code a dead-time
analyser and a frozen channel deserve, and units promoted from prose in
`measurements.rs` to data on a channel.

**The step size is not adapted.** Three integrators exist as of B-0053: fixed
step explicit Euler, classical RK4, and Dormand-Prince 5(4) with an embedded
error estimate. **Only Euler reproduces the original**, and it is the default.
A variable step changes when the discrete phases run, which is a decision about
fidelity rather than about numerics, so it has not been taken.

That comparison produced a result about the original rather than about the
port. RK4 and Dormand-Prince agree with each other to 1.5e-6 while both differ
from Euler by about 1.1e-2, so the published Tennessee Eastman data carries
roughly one percent of integration error against an accurate solution of the
same equations. Reproducing that is the point, but it means "the TEP" names a
particular discretisation and not only a set of differential equations.

## The two Class C quirks, and what a default `Scenario` now does

Both were signed off on 2026-08-28, so a default `Scenario` **fixes** them and
`Scenario::faithful()` reproduces them.

A trip ends the run (delta D-007, `teprob.f:807-811`). The original freezes the
plant and keeps reporting, which is where four of the forty-four published files
get their frozen tails: 1,832 rows in total, and 363 of `d06.dat`'s 480. The fix
is pure truncation, since every sample before the trip is bit-identical either
way, and the argument for keeping the freeze was that it preserved an option it
does not preserve, because the plant cannot be restarted in either case.

The driver does not force `IDV(12)` at hour eight (delta D-011,
`temain_mod.f:366-368`). Tier 7 established that the published files were
generated with that line replaced rather than kept: every `dNN_te` except
`d12_te` sits at the nominal operating point straight across row 160.

Every comparison against the Fortran or against published data runs
`Scenario::faithful()`, and `tier5::run_port` pins it so no differential can
lose it. From the command line the flags are `tep run --faithful`, or
`--force-idv12` and `--freeze-on-trip` individually.
