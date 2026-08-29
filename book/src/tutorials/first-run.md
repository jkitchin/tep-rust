# A first run, and the 53 channels

> **The worked example is
> [notebook 1, Getting started](../notebooks/01-getting-started.html).** It runs
> the plant, plots what normal operation looks like, walks the twenty faults,
> injects one, reads the ground truth back, checks reproducibility from a seed,
> and finishes with a trip. Source: `notebooks/01-getting-started.ipynb`.

The whole library is three objects. A `Scenario` says what to simulate, a
`Simulation` does it, and a `Run` holds what came out. The notebook uses all
three in its first cell; this page is the part of the story that is worth
stating once in prose rather than re-deriving from a plot.

The examples are Python, and [The Python package](../python.md) is how to
install it. The simulator itself is Rust, and its Rust API is
[Getting started](../getting-started.md), but a run started from either binding
is bit-identical, because both call the same core.

## Why forty samples in two hours, and not 7200

The integrator takes one step per simulated second, so two hours is 7200 steps.
A sample is written every 180 of them, which is the three-minute spacing
`temain_mod.f:401` writes at and the spacing of the published `d00` through
`d21` files. Two hours is therefore forty rows, and a 48-hour run is 960.

The first sample of a run is at 0.0497 hours rather than at 0.05, and that is
not an off-by-one. `run.steps[0]` on that row is 180, and `run.hours[0]` is the
time at which step 180 *began*, which is 179 seconds. The simulated clock is
advanced at the end of a step, after the row has been written, because that is
the order `temain_mod.f` writes in, and a row that carried the post-step time
would be labelled with a clock the plant had not reached when it was measured.
The last row of a two-hour run is at 1.9997 hours for the same reason. If you
want the initial condition itself, it is the model's nominal state and not a row
of the run.

## What the 53 channels are

A row is the plant as an operator sees it: 41 measurements and then the 12 valve
positions the controllers are holding. `run.to_numpy()` returns exactly that as
one `(n_samples, 53)` array, measurements first, and `channel_names()` returns
names in the same order, so a CSV header and a matrix column can never disagree
about which is which.

The 41 measurements split into two groups that behave quite differently.

`XMEAS(1)` through `XMEAS(22)` are continuous instruments: flows, pressures,
levels, temperatures, the compressor work. They are read every step and carry
Gaussian measurement noise whose standard deviation is the `XNS` table in
`teprob.f`.

`XMEAS(23)` through `XMEAS(41)` are the three gas chromatographs, and they are
not instruments in the same sense at all. They sample on a schedule, take time
to run, and then hold the answer until the next result arrives. The notebook
plots this, and it is the single most surprising thing about the data the first
time you see it: a composition channel is a staircase where a flow meter is a
curve. The reactor feed and purge analysers run every 0.1 hours and the product
analyser every 0.25, so at the three-minute output cadence each answer appears
twice or five times. A detector that treats those repeats as independent
observations is counting the same measurement several times, and that is worth
knowing before it produces a p-value.

Both `0.1` literals in `teprob.f` are single precision, so the gas interval is
really 0.10000000149011612 and a step landing on exactly 0.1 does not sample.
That is faithfully reproduced here; see [the delta register](../deltas.md).

`XMV(1)` through `XMV(12)` are the manipulated variables. Eleven of them move.
`XMV(12)`, the agitator speed, sits at 50 forever, because the published control
scheme never touches it. A statistical model fitted to all 53 channels has to
cope with that constant column, and [the detector tutorial](a-detector.md) shows
what a PCA model has to do about it.

## Getting at the data

`run.to_numpy()` is the record itself, one read-only `(n_samples, 53)` array of
`float64`. Three views on top of it cover most uses. `run.measurement(n)` and
`run.manipulated(n)` take the one-based indices of the original, so
`run.measurement(7)` is `XMEAS(7)` and needs no mental arithmetic.
`run.column(i)` takes the zero-based row position, and `run.columns()` returns
all 53 at once keyed by name, which is the shape a covariance matrix, a chart or
a `pandas.DataFrame` wants. Alongside them, `run.hours` and `run.steps` are the
simulated time and the integrator step each row was written at.

Every one of those is a view into the same buffer rather than a copy, so
`.copy()` is needed before writing into one.

## The digest

`scenario.digest` is a content hash over everything the run's output depends on:
the seed, the duration, the step, the cadence, the disturbances, the control
mode, the quirk flags, the schedule and the integrator. Two scenarios that
describe the same experiment produce the same sixteen characters, and two that
differ in any respect do not.

It is worth putting in a filename or a file header, because it turns "this is
the fault 4 data, I think" into something checkable. `scenario.to_text()` writes
the whole scenario out in one line that `Scenario.from_text` reads back, so a
dataset can carry its own description rather than a memory of one. The notebook
checks both properties rather than asserting them: the same scenario run twice
is bit-identical, a different seed is not, and the text round-trips.

## What a trip looks like

The plant has eight shutdown conditions, on reactor pressure, reactor level,
reactor temperature, separator level and stripper level. Cross one and the
simulation is over. `run()` never raises for a plant that misbehaves, because a
run that ended early is data, and throwing it away would hide the difference
between a port that trips where the original does and one that does not. The
outcome, the step, the hour and the condition that fired are all on the `Run`.

The notebook's survey of all twenty faults at 48 hours is the useful version of
this: at the default seed exactly one of them trips, `IDV(6)`, the total loss of
the A feed. That is a property of the seed as much as of the fault, which is the
kind of thing worth measuring before designing an experiment around it.

## The same run from the command line

```console
$ tep run --hours 2
```

writes the same 40 rows as CSV on stdout, with seventeen significant digits, so
the file round-trips an `f64` exactly.
