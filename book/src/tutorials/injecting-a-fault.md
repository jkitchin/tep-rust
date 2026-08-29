# Injecting a fault, and finding it in the data

> **The worked example is the second half of
> [notebook 1, Getting started](../notebooks/01-getting-started.html),** which
> switches `IDV(4)` on eight hours into a 24-hour run, plots the reactor
> temperature and the cooling water valve side by side, and reads the ground
> truth back. Source: `notebooks/01-getting-started.ipynb`.

A disturbance is one line: `Scenario.fault(n)` is the baseline with `IDV(n)`
switched on for the whole run, one-based exactly as the Fortran's `IDV(n)` is.
That is what the original does, and it is the right thing when what you want is
a faulted record. It is the wrong thing when what you want is to watch the fault
*arrive*, and arrival is what a detector is judged on, so the notebook schedules
it instead: `IDV(4)` off for eight hours and on afterwards, which is the layout
the published test sets use.

`IDV(4)` is a five degree step in the reactor cooling water inlet temperature.
The interesting part is that it is nearly invisible where you would first look.

## The fault is in the valve, not in the temperature

Over the fourteen hours after the fault settles, the mean reactor temperature is
120.399542 degrees on the fault-free run and 120.399514 degrees on the faulted
one. The difference is 28 millionths of a degree, against a fault-free standard
deviation of 0.018728 degrees: two thousandths of one standard deviation. To any
instrument, and to any detector watching that channel, nothing happened.

The cooling water valve tells the other half of the story. It sits at 41.103%
open without the fault and 44.868% with it, a shift of 3.766 percentage points
of valve travel. That is the whole of `IDV(4)`: the cooling water arrives
hotter, the temperature controller notices immediately, and it opens the valve
until the temperature comes back. The controller has converted a disturbance in
a measured variable into a disturbance in a manipulated one, which is what a
controller is for.

The lesson generalises well past this fault, and it is the single most important
thing to understand before building a monitor for a plant under closed-loop
control. The evidence of a disturbance often sits in the manipulated variables
rather than in the measurements, because the controller has been busy cleaning
up the measurements. Both halves are in the array `run.to_numpy()` returns, and
the detector in [the next tutorial](a-detector.md) uses all 53 channels for
exactly this reason.

## Ground truth

`run.labels()` records what was actually true at each instant. It returns two
`(n_samples, 20)` arrays indexed by `IDV(n) - 1`: `active` says whether that
disturbance was on at that sample, and `since_onset` says how many hours it had
been on, `nan` where it was not on at all.

The original records nothing of the sort. A published Tennessee Eastman dataset
is a matrix and a filename, so every detection delay in the literature is
measured against an onset its author knew from the experimental protocol rather
than from the data. That is fine until two papers disagree about where sample
160 falls. Here the onset is in the data, and finding it is
`int(np.argmax(labels["active"][:, 3]))` rather than a constant somebody has to
remember.

The precision is deliberate too. A fault live from the first step reports
0.04972222222222225 hours at the first sample, not "about 0.05": the label
carries the same simulated clock the row does, so it is 179 seconds. That is
what lets a delay of one sample be distinguished from a delay of zero.

## The disturbance you did not ask for, and no longer get

Pass `driver_forces_idv12=True` and run for longer than eight hours, and a
second fault appears at hour eight without being asked for. That is not a bug in
the scheduler. It is `temain_mod.f:366-368`, which switches `IDV(12)` on at
eight hours whatever the scenario said, and both `IDV(4)` and `IDV(12)` act on
cooling water, so a run nominally labelled `IDV(4)` was really the two together.

It was the default here until 2026-08-28, on the belief that reproducing `d01`
through `d21` required it. Tier 7 showed the opposite: every `dNN_te` file
except `d12_te` sits at the nominal operating point straight across row 160,
which is hour eight, so the published files were made with that line replaced.
The default now follows the files.

It is delta D-011 in [the register](../deltas.md), and it used to be the single
most common way a comparison against the published files went quietly wrong,
because it was the default. It is not any more: a request for `IDV(4)` gets
`IDV(4)`, and `driver_forces_idv12=True` is how to ask for the driver's version.
The Rust equivalent is `Scenario::faithful()`, which sets that quirk and the
freeze-on-trip one together, and the command line spelling is
`tep run --force-idv12`. Say which one you used when you report numbers. The
labels make the difference visible either way, which is the point of recording
ground truth rather than assuming it.

## From the command line

```console
$ tep run --fault 4 --hours 8 --labels
```

`--labels` adds the `fault` and `hours_since_onset` columns to the CSV, so the
ground truth travels with the data.
