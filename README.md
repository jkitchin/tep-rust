# TEP-Rust

A pure-Rust port of the Tennessee Eastman Process, taken from Downs and Vogel's
original Fortran rather than from any later reimplementation, and shipped three
ways: a Rust crate, a Python package with no C dependencies, and a
self-contained WebAssembly app that runs in the browser with no server.

TEP is the standard benchmark for plant-wide process control and for
fault-detection research: a reactor, condenser, vapor-liquid separator, recycle
compressor and stripper, with eight chemical species, four reactions, 50 states,
41 measurements, 12 manipulated variables and 20 canonical disturbances.

## The headline result

Given the same `exp` and `pow`, a complete 48-hour closed-loop run of **172,800
integrator steps is bit-identical to the Fortran** in all 41 measurements and
all 12 manipulated variables. Not one measurement, not one valve, not one step.

That comparison is run in CI against the original Fortran as a live oracle, not
against a recorded trace. Under the vendored `libm` the two disagree only where
`exp` and `pow` disagree, and over those same two simulated days the worst error
is 1.705e-10 on `XMEAS(1)` and 1.246e-11 on `XMV(10)`, which is under a tenth of
a billionth of any instrument's noise band.

## Status

**Complete.** Every phase has landed and every design decision has been taken.
The backlog is 85 items done, none open and none blocked.

## Quick start

The command line:

```console
$ cargo run -p tepsim-cli -- run --fault 4 --hours 24 --labels > idv4.csv
$ cargo run -p tepsim-cli -- dataset --out ./data
$ cargo run -p tepsim-cli -- faults
```

The library:

```rust
use tepsim::{Scenario, Simulation};

let run = Simulation::new(Scenario::fault(4).with_hours(24.0)).run();
println!("{:?} after {} samples", run.outcome, run.samples.len());
```

Python, from a wheel that carries no C dependency:

```python
import tepsim as tep

run = tep.Simulation(tep.Scenario.baseline(seed=42, hours=48)).run()
run.to_numpy()        # (960, 53) float64, XMEAS(1..41) then XMV(1..12)
run.measurement(7)    # XMEAS(7), reactor pressure
run.outcome           # 'completed'
```

The browser app is `apps/studio`: simulation on a Web Worker, charts and an SVG
flowsheet on the main thread, scenarios shareable as a URL fragment. It is
74,010 bytes gzipped, under 5% of the project's 1.5 MB budget, and runs between
40,662 and 85,529 times real time.

## How it is validated

A ten-tier ladder, all of it implemented, running against the Fortran as a live
oracle. `cargo xtask ci` is the gate and `cargo xtask validate` is the periodic
deep check.

| Tier | What it proves | Where it stands |
|---|---|---|
| 1 | `TESUB1`-`TESUB8` match the oracle | relative error under 1e-13, and `TESUB7` bit-exact |
| 2 | Single-step derivatives match, measured against the scale of the terms | all 50 components over three sampling pools |
| 3 | RNG call *order* matches | trace diff empty |
| 4 | Trajectories (diagnostic, not a gate) | 48 h bit-identical under the platform libm |
| 5 | Statistical equivalence | 21 scenarios by 100 seeds by 48 h, 4,200 runs; 19 of 21 pass everything, and the other two are the stuck-valve case below |
| 6 | Downstream detectors cannot tell the sources apart | cross-source within the reference's own split-half variation |
| 7 | Reproduces the published `d00`-`d21` | measured per file, with every unknown protocol detail named rather than fitted |
| 8 | Differential fuzzing with shrinking | 5,000,000 tuples, 2,551,618 compared, one counterexample, attributed to the vendored `libm` |
| 9 | Cross-platform determinism | six committed digests, identical on aarch64 and wasm32; suite digest `953892e861fd5e68` |
| 10 | Every quirk fix has a measured delta | 11 register entries, each with a marker in the source that CI cross-checks |

Two `IDV` faults miss the Tier 5 margin under the vendored `libm` and are
bit-identical under the platform one. The cause is established rather than
assumed: `teprob.f:801` moves a valve only when the command has travelled past a
dead band, which is a discontinuous branch, so one ULP of `exp` decides which
side it lands on and the valve then holds a different position for the rest of
the run. Those valves are judged on their distribution rather than on a mean,
because a series of plateaux has no meaningful centre. No margin was widened.

## Faithful by default, and honest about the rest

Every deliberate deviation from the Fortran is a numbered entry in
[`book/src/deltas.md`](book/src/deltas.md) with a marker in the source, and
`cargo xtask deltas` fails if either half goes missing or the two disagree about
the class. There are eleven.

Two of them are behaviour-defining and were decided rather than inherited. A
trip ends the run, where the original freezes the plant and keeps reporting; and
the driver does not force `IDV(12)` on at hour eight, because the published
files turn out to have been generated with that line replaced. `Scenario::faithful()`
reproduces both, `tep run --faithful` is the same thing from the shell, and
every comparison against the Fortran or against published data pins it, so no
differential can quietly lose one.

Beyond fidelity, faults are data rather than code: the 20 disturbances are
`(injection point, profile)` pairs, so they can be scheduled, composed, given
continuous magnitudes, or defined from scratch, and every sample carries
machine-readable ground truth so detection delay stops being guesswork.

## Layout

| Crate | What it is |
|---|---|
| `tepsim` | The API most callers want: scenarios, runs, recorders, dataset output |
| `tepsim-core` | The plant model. `no_std`, `forbid(unsafe_code)`, deterministic |
| `tepsim-control` | Controllers and the multi-rate scheduler |
| `tepsim-scenario` | Disturbances, schedules and reproducible scenario digests |
| `tepsim-cli` | `tep`, the command line |
| `tepsim-py` | PyO3 bindings and the wheel |
| `tepsim-wasm` | WebAssembly bindings behind the browser app |
| `tepsim-oracle` | Development only. FFI harness driving the original Fortran |
| `tepsim-stats` | Development only. The statistics behind the ladder, in Rust |

Determinism is a hard invariant: no `f32`, no SIMD or rayon in the core, no
reordered reductions, and no clock or randomness outside the seeded generator.
`tepsim-oracle` is never a dependency of anything shipped, and CI asserts it.

## Documents

| File | What it is |
|---|---|
| [`PLAN.org`](PLAN.org) | Design of record: architecture, validation strategy, roadmap |
| [`BACKLOG.org`](BACKLOG.org) | Ordered work queue and current state |
| [`LOG.org`](LOG.org) | Iteration history with the measured validation numbers |
| [`book/`](book/) | mdBook: theory, unit operations, tutorials, delta register, generated validation report |
| [`CLAUDE.md`](CLAUDE.md) | Development protocol for this repository |
| [`NOTICE.md`](NOTICE.md) | Attribution, upstream license, and citation requirements |

## Building

Rust 1.97.1, pinned in `rust-toolchain.toml`. `cargo xtask ci` is the full gate;
`cargo xtask ci --fast` skips the oracle differential for local iteration. The
oracle needs `gfortran`, the Python step needs `python3` and `maturin`, and
Tier 9's wasm leg needs the `wasm32-unknown-unknown` target and Node. Each step
skips with a clear message when its toolchain is absent.

## License

MIT ([`LICENSE`](LICENSE)). Portions are derived from Fortran licensed under the
University of Illinois/NCSA Open Source License ([`LICENSE-NCSA`](LICENSE-NCSA)),
whose attribution conditions apply to source and binary redistribution alike and
cannot be dropped. The combined work is `MIT AND NCSA` in SPDX terms. See
[`NOTICE.md`](NOTICE.md) before redistributing.

The original process and control problem are due to J. J. Downs and E. F. Vogel
(Tennessee Eastman Company); the widely used modified code is due to the Large
Scale Systems Research Laboratory at the University of Illinois under Prof.
Richard D. Braatz.
