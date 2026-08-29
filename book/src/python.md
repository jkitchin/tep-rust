# The Python package

The simulator ships as a Python package called `tepsim`. It is not a
reimplementation and not a wrapper around a subprocess: it is the same Rust core
the rest of this book documents, compiled into an extension module with PyO3.
There is no C in the build, no Fortran, and no runtime dependency except NumPy.

A run started from Python is bit-identical to the same run started from Rust or
from the browser, because all three call the same code with the same scenario.
That is what makes the worked examples in the notebooks usable as evidence
rather than as illustrations.

## Installing it

It is **not on PyPI**. The version is still `0.0.0` and no release tag has been
pushed, so nothing has ever been published. The `publish` job in
`.github/workflows/wheels.yml` is wired up and gated on a `v*` tag, and the name
is free, so `pip install tepsim` will be the line once there is a release. It is
not the line today.

Install it from the repository instead:

```console
$ pip install "git+https://github.com/jkitchin/tep-rust#subdirectory=crates/tepsim-py"
```

Two things about that command, both of which will bite otherwise.

**It builds from source, so you need a Rust toolchain.** There is no published
wheel to fall back on, so `pip` clones the repository, invokes maturin, and
compiles the whole core in release mode. That takes a couple of minutes on a
warm machine and it fails partway through if `cargo` is not on your `PATH`.
Install Rust from [rustup.rs](https://rustup.rs/) first. NumPy is pulled in
automatically as a declared dependency, so nothing else is needed.

**The quotes are not optional.** Most shells treat `#` as the start of a
comment, and without the quotes the `#subdirectory=` fragment is silently
dropped, at which point `pip` tries to build the workspace root and fails with
an error that does not mention the fragment at all.

The result is a `cp39-abi3` wheel, so one build covers every GIL-enabled CPython
from 3.9 upwards. Free-threaded interpreters have no stable ABI to target and
need a version-specific build; 3.14t is the earliest that works, because PyO3
0.29 does not support the free-threaded build of anything below 3.14.

### From a checkout

If you have the repository already, `cargo xtask python` is the supported path.
It builds a release wheel, creates a throwaway virtualenv at
`.xtask-python/venv`, installs the wheel and NumPy into it, proves that
`import tepsim` resolves inside that virtualenv rather than somewhere else on
the machine, and runs the binding's pytest suite against it.

```console
$ cargo xtask python
$ .xtask-python/venv/bin/python -c "import tepsim; print(tepsim.__version__)"
```

That interpreter is the one the book's own tests use to check the quickstart
transcript below. Note that `cargo xtask python` deletes and rebuilds the
virtualenv every time it runs, so anything else installed into it, `jupyter` and
`matplotlib` for the notebooks in particular, has to be reinstalled afterwards.

To build a wheel without installing it, `maturin build --release -m
crates/tepsim-py/Cargo.toml` writes one and prints the path.

## Quickstart

```python
from concurrent.futures import ThreadPoolExecutor

import numpy as np

import tepsim as tep

print("tepsim %s: XMEAS(1..%d), XMV(1..%d), %d channels, IDV(1..%d)"
      % (tep.__version__, tep.MEASUREMENTS, tep.MANIPULATED, tep.CHANNELS,
         tep.DISTURBANCES))

# A Scenario says what to simulate, a Simulation does it, a Run holds what
# came out. That is the whole API.
run = tep.Simulation(tep.Scenario.baseline(seed=42, hours=48)).run()
print()
print("run:      %r" % run)
print("matrix:   %s %s" % (run.to_numpy().shape, run.to_numpy().dtype))
print("outcome:  %s" % run.outcome)
print("XMEAS(7): mean %.2f kPa over %.0f h"
      % (run.measurement(7).mean(), run.hours[-1]))

# The twenty disturbances say what they do, not only what the original header
# called them. Five of the published descriptions are the word "Unknown".
print()
print("the five the original leaves unexplained")
for fault in tep.faults():
    if fault.published == "Unknown":
        print("  IDV(%2d) %-9s %s" % (fault.index, fault.shape, fault.effect))

# Ground truth travels with the data, which the original records nowhere.
faulted = tep.Simulation(tep.Scenario.fault(1, hours=8)).run()
labels = faulted.labels()
print()
print("IDV(1) over 8 h: active at the last sample %s, %.2f h since onset"
      % (labels["active"][-1, 0], labels["since_onset"][-1, 0]))

# run() releases the GIL, so an ensemble is a thread pool and nothing has to be
# pickled.
sims = [tep.Simulation(tep.Scenario.fault(n, hours=8)) for n in range(1, 21)]
with ThreadPoolExecutor() as pool:
    runs = list(pool.map(tep.Simulation.run, sims))
print()
print("twenty 8-hour faulted runs: %d completed, %d tripped"
      % (sum(r.outcome == "completed" for r in runs),
         sum(r.outcome == "tripped" for r in runs)))

# A run is a pure function of its scenario, and a scenario is one line of text
# that parses back to an equal scenario.
scenario = tep.Scenario.fault(4, hours=8)
print()
print("digest:   %s" % scenario.digest)
print("text:     %s" % scenario.to_text())
print("parses back equal:      %s"
      % (tep.Scenario.from_text(scenario.to_text()) == scenario))
print("two runs bit-identical: %s"
      % np.array_equal(tep.Simulation(scenario).run().to_numpy(),
                       tep.Simulation(scenario).run().to_numpy()))
```

```text
tepsim 0.0.0: XMEAS(1..41), XMV(1..12), 53 channels, IDV(1..20)

run:      <tepsim.Run 960 samples x 53 channels, 48.0 h, completed>
matrix:   (960, 53) float64
outcome:  completed
XMEAS(7): mean 2706.16 kPa over 48 h

the five the original leaves unexplained
  IDV(16) random    enables walk channel 9, the stripper steam valve capacity
  IDV(17) random    enables spike channel 10, the reactor coolant duty
  IDV(18) random    enables spike channel 11, the condenser coolant duty
  IDV(19) sticking  sticks valves 5, 7, 8 and 9; touches no equation in the model
  IDV(20) random    enables spike channel 12, the reactor outlet flow

IDV(1) over 8 h: active at the last sample True, 8.00 h since onset

twenty 8-hour faulted runs: 19 completed, 1 tripped

digest:   7f9accc04bb61e7f
text:     tepsim.scenario.v1;seed=4651207995;hours=8;step=2.777777777777778e-4;every=180;faults=4;controlled=1;idv12=0;trip=1;continuous=0;integrator=euler;events=
parses back equal:      True
two runs bit-identical: True
```

## The API

Three classes carry the whole surface, and they divide the way the problem
does. A `Scenario` is a description and holds no state. A `Simulation` is the
machinery. A `Run` is the result.

### Module level

| Name | Value | Meaning |
|---|---|---|
| `MEASUREMENTS` | 41 | `XMEAS(1..41)` |
| `MANIPULATED` | 12 | `XMV(1..12)` |
| `CHANNELS` | 53 | a row: measurements then manipulated variables |
| `DISTURBANCES` | 20 | `IDV(1..20)` |
| `DEFAULT_SEED` | 4651207995.0 | the generator word compiled into `teprob.f:1187` |
| `DEFAULT_STEP_HOURS` | 1/3600 | one simulated second, the step `INTGTR` uses |
| `DEFAULT_SAMPLE_EVERY` | 180 | the three-minute spacing of `d00` through `d21` |
| `FORCED_DISTURBANCE_STEP` | 28800 | the step the driver forces `IDV(12)` on at |

`channel_names()` returns the 53 names in row order, so a CSV header and a
matrix column cannot disagree about which is which. `faults()` returns the
twenty `Fault` records.

### `Scenario`

Immutable and cheap to copy, so a caller can build one, derive variants from it,
and keep all of them.

There are four constructors: the bare `Scenario(...)`, `Scenario.baseline(...)`,
`Scenario.fault(n, ...)` and `Scenario.from_text(text)`. All but the last take
the same keyword arguments, which are `seed`, `hours`, `step_hours`,
`sample_every`, `controlled`, `driver_forces_idv12` and `trip_ends_the_run`,
plus `faults` on the bare one.

| Member | What it is |
|---|---|
| `seed`, `hours`, `step_hours`, `sample_every` | the four numbers a run is measured in |
| `controlled` | closed loop under the published control scheme, or open loop with the valves held |
| `driver_forces_idv12` | whether the driver switches `IDV(12)` on at hour eight, delta D-011 |
| `trip_ends_the_run` | whether a shutdown stops the run, delta D-007 |
| `faults` | the active disturbances, one-based and ascending |
| `steps`, `samples` | how many integrator steps and how many recorded rows |
| `digest` | a content hash over everything affecting the run, sixteen hex characters |
| `to_text()`, `from_text()` | the canonical one-line form, and its strict parser |
| `with_seed`, `with_hours`, `with_fault`, `sampling_every`, `open_loop` | derive a variant |

The two defaults worth knowing are that a trip ends the run and that the driver
does not force `IDV(12)`. Both are Class C quirks of the original that this port
fixes by default and can reproduce on request; see
[the delta register](deltas.md).

### `Simulation`

`Simulation(scenario)` holds a plant, a controller stack and an integrator
state, all of it owned. The original keeps its whole working set in six Fortran
`COMMON` blocks, which allows exactly one simulation per process and no
reentrancy at all. This allows as many as there are threads.

`run()` runs the whole scenario and returns a `Run`. It runs a *copy*, so
calling it twice gives two equal runs rather than one run and one empty one, and
it never raises for a plant that misbehaves: a trip or a failed temperature
solve is reported through `Run.outcome`, because a run that ended early is data.

`run()` releases the GIL for the entire integration, which is where all the time
goes, so an ensemble is a `ThreadPoolExecutor` and not a `ProcessPoolExecutor`.
Nothing has to be pickled, and the plant does no I/O and touches no Python
object while the GIL is down.

### `Run`

| Member | What it is |
|---|---|
| `to_numpy()` | the whole run, one `(n_samples, 53)` float64 array, C-contiguous |
| `hours`, `steps` | simulated time and integrator step per row |
| `measurement(n)`, `manipulated(n)` | one channel, one-based as `XMEAS(n)` and `XMV(n)` are |
| `column(i)`, `columns()` | one channel by zero-based position, or all 53 keyed by name |
| `labels()` | ground truth: `'active'` and `'since_onset'`, both `(n_samples, 20)` |
| `outcome` | `'completed'`, `'tripped'` or `'solve_failed'` |
| `tripped_at`, `tripped_hours`, `trip_cause` | where and why the plant shut down, or `None` |
| `solve_failed_at` | the step a temperature solve failed at, or `None`; delta D-001 |
| `scenario` | what was asked for |

Every array is a read-only view over one buffer, which is filled once when the
run finishes and moved into NumPy rather than copied. `to_numpy()` returns the
same object every call, and `measurement`, `manipulated`, `column` and `columns`
are strided views into it. Call `.copy()` for something writable, or
`numpy.ascontiguousarray` if contiguity matters. `columns()` is a dict in
channel order, which makes it a `pandas.DataFrame` constructor argument as it
stands.

`labels()` is the part with no counterpart in the original. A published
Tennessee Eastman dataset is a matrix and a filename, so every detection-delay
figure in the literature is computed against whatever onset its author assumed.
Here the onset is in the data.

### `Fault`

`faults()` returns twenty of these. `published` is the description at
`teprob.f:172-191` verbatim, five of which say only "Unknown". `effect` is what
the source actually does, which for those five is perfectly explicit: only the
physical interpretation was withheld. `shape` is `'step'`, `'random'` or
`'sticking'`, `line` names the `teprob.f` line the fault acts on, and
`affects_the_plant` is `False` for the three sticking faults, which touch no
equation in the model at all. The same table is in
[The twenty disturbances](disturbances.md).

## What `Scenario(...)` cannot say

Three of a scenario's eleven fields are out of the constructors' reach: the
event schedule, the `continuous` extension that allows a fault at a fraction of
its full strength, and the choice of integrator. `Scenario.from_text` is the way
to them from Python, and since the text is the scenario's real serialisation,
the thing the digest is taken over and the browser app puts in a URL fragment,
building through it is building in the form the run will be described by.

```python
text = tep.Scenario.baseline(hours=18).to_text()
scheduled = tep.Scenario.from_text(text.replace("events=", "events=6:start:4,12:stop:4"))
```

`repr()` knows about this. A scenario the constructor can express is printed as
a constructor call, because that is what is readable at a prompt, and anything
else is printed as `Scenario.from_text(...)`. Both round-trip, so
`eval(repr(s)) == s` holds for either shape. The alternative, printing only the
constructor's arguments, produced valid Python that evaluated to a *different*
scenario with a different digest and said nothing about it.

`notebooks/04-custom-scenarios.ipynb` works through all three fields, with a
helper that rebuilds the text line from a dictionary rather than patching it, so
that a mistyped field name is an error rather than a silently different run.

## What is not in the package

`tepsim-stats`, the crate the repository's own validation ladder runs its
statistics on, is development-only. `cargo xtask ci` asserts that no shipped
crate so much as names it, and it is therefore not reachable from Python. That
is deliberate: the ladder's statistics have to be free of any dependency the
port itself does not have, so they can never be the reason a validation number
changed.

For monitoring work the notebooks use `notebooks/pcamon.py` instead, which is
PCA, both monitoring statistics and their control limits in NumPy and the Python
standard library alone, with no SciPy and no scikit-learn. Where the two could
differ they agree: `02-fault-detection-pca.ipynb` checks every printed digit of
its detector against the Rust one.

## Types

The package ships `py.typed`, and `_tepsim.pyi` beside `__init__.py` describes
the compiled half, so `mypy` and an editor see the real signatures rather than
`Any`.

## Worked examples

The four notebooks in `notebooks/` are where the package is actually used. They
are executed and committed with their outputs and plots, and rendered copies are
published beside this book:

- [Getting started](notebooks/01-getting-started.html)
- [Fault detection with PCA](notebooks/02-fault-detection-pca.html)
- [The three hard faults](notebooks/03-hard-faults.html)
- [Custom scenarios](notebooks/04-custom-scenarios.html)

Those four links resolve on the published site, where
`.github/workflows/pages.yml` renders the notebooks beside the book with
`jupyter nbconvert`. In a local `mdbook build` they point at files that are not
there; the sources are `notebooks/*.ipynb` in the repository, and
`notebooks/README.md` says how to run them.

The [tutorials](tutorials/first-run.md) in this book are the narrative around
them.
