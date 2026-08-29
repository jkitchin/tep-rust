# Getting started

This page is the Rust API and the command line tool. If you want to drive the
simulator from Python, [The Python package](python.md) is the equivalent page
for the `tepsim` wheel, and the [tutorials](tutorials/first-run.md) and the
notebooks they link to are all Python.

## Building

The workspace pins its toolchain in `rust-toolchain.toml`, so a `rustup`
install picks the right compiler on its own. Nothing in the shipping crates
needs Fortran; `gfortran` is required only to build `tepsim-oracle`, the
development-only differential harness, which sits behind the `oracle` feature
and is never a dependency of `tepsim`.

```console
$ cargo build --release
$ cargo xtask ci          # fmt, clippy, tests, docs, cargo-deny
```

## The library

The whole public surface is three types. A `Scenario` describes a run, a
`Simulation` performs one, and a `Run` holds the output as columns.

```rust,ignore
use tepsim::{Scenario, Simulation};

// Twenty-four hours with disturbance IDV(4), the reactor cooling water
// inlet temperature step.
let run = Simulation::new(Scenario::fault(4).with_hours(24.0)).run();

// One measurement across the run, one-based as XMEAS(n) is.
let reactor_pressure = run.measurement(7);

// One manipulated variable, one-based as XMV(n) is.
let coolant_valve = run.manipulated(10);

assert!(run.outcome.is_completed());
```

`Scenario` is a plain `Copy` struct, so a caller can build one, tweak it, and
keep both. Its builders are `with_seed`, `with_hours`, `sampling_every`,
`with_fault` and `open_loop`, and the two constructors are
`Scenario::baseline()` and `Scenario::fault(n)`.

| Field | Default | Meaning |
|---|---|---|
| `seed` | 4651207995 | the generator word compiled into `teprob.f:1187` |
| `hours` | 48.0 | `NPTS = 172800` at a one-second step, the run `temain_mod.f` was written to do |
| `step_hours` | 1/3600 | the step the original's `INTGTR` uses |
| `sample_every` | 180 | `temain_mod.f:401` writes every 180 seconds, and `d00` through `d21` are at that spacing |
| `disturbances` | all off | twenty flags, one-based as `IDV(n)` is |
| `controlled` | `true` | closed loop under the published control scheme, or open loop with the valves held |
| `quirks` | all fixed | which Class C quirks are fixed rather than reproduced |
| `driver_forces_idv12` | `false` | whether the driver switches `IDV(12)` on at hour eight, as `temain_mod.f:366-368` does |

Three of those defaults deserve a sentence. Open loop is not a useful operating
mode, it is a diagnostic one: the plant trips on reactor pressure after about
three hours, and the difference between the two settings is the clearest single
statement of what the control layer does.

The other two are the Class C quirks, signed off on 2026-08-28 and fixed by
default: a trip ends the run (delta D-007) and the driver does not force
`IDV(12)` (delta D-011). Both are one call away.
`Scenario::faithful()` gives the original's behaviour on both, and it is what
every comparison against the Fortran or against published data runs. See
[the delta register](deltas.md) for the evidence behind each.

A `Run` holds `samples`, each carrying `step`, `hours`, `measurements`
(`XMEAS(1..41)`), `manipulated` (`XMV(1..12)`) and `labels`. `Run::column`,
`Run::columns`, `Run::measurement` and `Run::manipulated` reshape it into the
columns a statistic or a detector wants, and `Sample::row` gives all 53 channels
in one array, measurements first, which is the layout every downstream consumer
uses. `channel_names()` returns matching names in the same order.

`Run::outcome` is how a run ended: `Completed`, `Tripped` with the step, the
hour and the first shutdown condition that fired, or `SolveFailed` with the step
at which a temperature solve failed to converge. A trip **ends the run** by
default, so a tripped run is shorter than its scenario planned.
`teprob.f:807-811` instead freezes the plant and keeps reporting, which is what
`Scenario::faithful()` gives and what the constant tails in `d06` and `d18`
are.

`Labels` is ground truth: which disturbances were active at that instant, and
how long each had been. The original records nothing of the sort, so every
detection-delay figure in the literature is computed against an onset the author
assumed.

## Stepping a run by hand

`Simulation::run()` is the whole scenario in one call. For an online detector, a
live display or a reinforcement learning loop, `Simulation::step()` advances one
integrator step and hands back a `Sample` on the steps where one is due.

```rust,ignore
use tepsim::{Scenario, Simulation};

let mut sim = Simulation::new(Scenario::baseline().with_hours(2.0));
while let Some(sample) = sim.step() {
    println!("{:.3} h  pressure {:.1}", sample.hours, sample.measurements[6]);
}
```

The order inside a step is `temain_mod.f`'s: force `IDV(12)` if it is due, run
the controllers on the *previous* step's measurements, integrate, then clamp.
That one plant step of dead time in every loop is delta D-010, and getting it
backwards leaves `XMEAS(14)` 23% out after four hours.

## The command line

```console
$ tep run --fault 4 --hours 24 --labels
```

`tep` has four subcommands, `run`, `dataset`, `faults` and `help`, and `run`
takes these flags:

| Flag | Default | Effect |
|---|---|---|
| `--fault <1-20>` | none | inject a disturbance |
| `--hours <h>` | 48 | simulated duration |
| `--seed <n>` | 4651207995 | generator word |
| `--every <steps>` | 180 | sample every N steps, so three minutes at the default step |
| `--open-loop` | off | hold the valves instead of controlling |
| `--force-idv12` | off | switch `IDV(12)` on at hour eight whatever was asked for, as `temain_mod.f:367` does |
| `--freeze-on-trip` | off | on a trip, freeze the plant and keep reporting rather than ending the run, as `teprob.f:807-811` does |
| `--faithful` | off | both of the above: reproduce every Class C quirk |
| `--labels` | off | include the ground-truth columns |

`tep faults` prints the twenty disturbances with the shape of each and the
`teprob.f` line it acts on; the same table appears in
[The twenty disturbances](disturbances.md).

Output is CSV on stdout, with progress and the outcome on stderr, so that
redirecting stdout to a file gives a clean file. The header is `step,hours`
followed by the 53 channel names, plus `fault` and `hours_since_onset` when
`--labels` is given. Values are written with seventeen significant digits, which
round-trips an `f64` exactly, so a CSV written here is reproducible rather than
approximately reproducible.

A trip is reported on stderr and exits successfully, because a trip is a result
rather than a malfunction; the message says whether the run ended there or the
plant froze and carried on. Only a failed temperature solve is an error exit.

### Generating a dataset

```console
$ tep dataset --out ./data
$ tep dataset --set all --faults 0,4,6 --format csv --out ./small
$ tep dataset --list
```

`tep dataset` writes files with the same geometry as the published `d00`
through `d21`: 960 rows of 52 columns for the `_te` half, 480 (500 for `d00`)
for the training half, at the seeds `teprob.f:1187-1256` records, in the same
fixed sixteen-column layout the shipped files use. `--set` chooses `te` (the
default), `training` or `all`; `--faults` takes `all` or a comma list;
`--format` is `dat` or `csv`, the latter at full precision with a header.

It is emphatically **not** the published data, and it says so on stderr every
run. Four things stand in the way and none is a matter of trying harder: the
toolchain that made the shipped files is unrecorded and its `exp` was not this
one's, the output was rounded to five significant figures before anyone saw it,
`IDV(21)` is not in this revision of the model so `d21` cannot be generated at
all, and the protocol behind the twenty-two training files is written down
nowhere. [Tier 7](validation/index.md) measures the gap rather than closing it.

## Reproducibility

A run is a pure function of its `Scenario`. There is no clock, no thread-local
state and no global, so the same scenario gives the same bits on x86-64,
aarch64 and wasm32. That is what makes a recorded dataset reproducible from its
description rather than from a file, and it is the property Tier 9 will assert
across the CI matrix.
