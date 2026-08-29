# A fault that arrives, and clears

> **The worked example is
> [notebook 4, Custom scenarios](../notebooks/04-custom-scenarios.html),** which
> schedules `IDV(4)` between hours six and twelve and plots the valve moving and
> coming back, composes two faults, sweeps a fault's magnitude from zero to one,
> and compares the three integrators. Source:
> `notebooks/04-custom-scenarios.ipynb`.

Every published Tennessee Eastman dataset is the same shape of experiment: set
some `IDV` flags before the run, then leave them. That is all the original
admits, which is why the literature's fault onsets are all in the same place and
why almost nobody studies a fault that arrives, persists, and then goes away.

A `Scenario` here carries a schedule of events, each at a time and each doing
one thing, so a disturbance can arrive at a stated hour, several can arrive
independently of each other, and any of them can clear. Two further things the
original cannot express sit beside it: a fault applied at a fraction of its full
strength, and a choice of integrator. All three live in the scenario's canonical
text form, which is where Python reaches them, and none of it would be worth
much if the description of a run could not travel with the run.

## The schedule is text, and that is the point

The Python constructors cover eight of a scenario's eleven fields. The other
three, `events`, `continuous` and `integrator`, are reached by rendering a
scenario to text, editing a field, and parsing it back. That sounds like a
workaround and is closer to the opposite. The text is the scenario's real
serialisation: it is what the digest is taken over, what the browser app puts in
a URL fragment, and what the Rust and wasm sides read and write, so building a
schedule through it is building it in the form the run will be described by.

An event is `time:verb` followed by the verb's own fields, so an events field of
`6:start:4,12:stop:4` switches `IDV(4)` on at hour six and off at hour twelve.
The notebook's helper rebuilds the whole line from a dictionary rather than
patching it with a string replacement, so a mistyped field name is a `KeyError`
and not a silently different run.

An event is applied on the first step whose time is at or past its own. The
window is half-open at the start and closed at the end, so an event is applied
exactly once however the step size divides its time, and an event at time zero
is applied on the first step rather than never. Two events at the same instant
keep the order they were written, because `stop` then `start` on one fault is a
different scenario from the reverse and the digest has to be able to tell them
apart.

## What the trace shows

The reactor cooling water valve sits at 41.063% open before the fault. Half an
hour after it arrives the valve is at 45.492%, it holds near 44.5% for the six
hours the fault is live, and half an hour after it clears it is back at 41.045%,
within a fiftieth of a percentage point of where it started. The controller does
not know a fault happened in either direction; it is rejecting a disturbance
both times.

That return is worth noticing rather than assuming. It is this loop doing its
job on a step disturbance that was removed, and it is not guaranteed in general:
a fault that shifts an inventory leaves the plant somewhere else even after it
clears, and the plant after a fault is not the plant before it, because the
controllers have moved. The recovery half of the problem is expressible here and
barely studied, which is most of the reason the schedule exists.

The labels follow the schedule exactly. Once the fault clears, `since_onset`
goes back to `nan`, which means the ground truth describes the plant's *current*
condition rather than its history. A study of recovery has to look at
transitions in `active`, not at that column.

## Composition

Several disturbances can be described independently, each with its own arrival
time. The original allows more than one `IDV` flag to be set, but not when each
arrives, and interactions are usually where the interesting behaviour is.

`IDV(1)` steps the A/C feed ratio and `IDV(8)` enables random walks on the A and
B feed compositions. Individually the plant absorbs both: reactor pressure means
of 2709.58 and 2702.19 kPa against a fault-free 2706.89. Together the mean is
2731.92 with a standard deviation of 113.07, where adding the two individual
shifts would predict 2704.87. The effect of the pair is not the sum of the
effects of the parts, which is the whole reason to be able to compose them.

## Half a fault

`teprob.f:341-346` opens `TEFUNC` by forcing every `IDV` to zero or one, so the
original has exactly two states per disturbance. The disturbances are then used
multiplicatively, `XST(1,4) = TESUB8(1,TIME) - IDV(1)*0.03` at `teprob.f:407`
being the pattern, so a magnitude between zero and one scales a fault smoothly
and needs no other change to the model.

That is the `continuous` extension, and it is off by default. With it on, a
sweep of `IDV(4)` from magnitude 0 to 1 moves the mean cooling water valve from
41.023% to 44.779% in steps of about 0.9 percentage points, linear in the
magnitude to within 0.042 percentage points of valve travel, while the reactor
temperature it is defending stays at 120.400 degrees at every magnitude. A
detection threshold study is then a sweep over one number rather than an
argument about what "harder" means.

Two things are worth saying plainly. A run with the extension on is not
comparable to any published dataset, and none of the validation ladder applies
to it. And a fractional magnitude without the extension is *refused* rather than
rounded, with an error naming the line of Fortran that makes it impossible.
Silently turning a request for half a fault into a whole one would produce a run
that does not match its own description, which is exactly what the content hash
exists to prevent. A magnitude of exactly 1.0 is bit-identical to the faithful
path, which is asserted, so the extension can be left on for a study that mixes
full and partial faults.

## The trap in refining the step

The notebook also compares Euler, RK4 and Dormand-Prince, and the result is
about the original rather than about the port: RK4 and Dormand-Prince agree with
each other to about 2e-6 relative while both differ from Euler by about 1.5e-2
on the worst channel. Two independent methods agreeing that closely and
disagreeing with the third is what convergence looks like. Euler is not the
accurate choice, it is the *faithful* one, and everything the validation ladder
claims is a claim about Euler.

The obvious next move, keeping Euler and halving the step, does not do what you
expect. Halving the step moves the answer by more than changing the method does,
and it does not converge as the step shrinks. That is neither a bug nor
stiffness: the disturbance walks and the measurement noise advance once per
*step*, so a run at half the step draws twice as many random numbers and is a
different realisation of the stochastic forcing. The step size is part of the
disturbance model in this plant and not only a numerical parameter. Change the
method to integrate the same realisation more accurately; change the step only
when you mean to change the noise.

## The description travels with the run

A scenario, schedule included, writes out as one line of text, that line parses
back to an equal scenario, and the digest is unchanged across the round trip. A
dataset generated from a scenario can carry the line in its header, and a reader
can reconstruct the exact run rather than the description of a run. The same
line is what TEP Studio puts in a URL fragment, so a scheduled experiment can be
handed to somebody as a link.

The format is versioned and its parser is strict: a missing field, an unknown
field, a value out of range or a version this build does not know is an error
that says what was wrong, rather than a default quietly substituted.

`repr()` handles the two shapes separately, and knowing why is useful. A
scenario the constructors can express prints as a constructor call, because that
is what is readable at a prompt. Anything else prints as
`Scenario.from_text(...)`, which round-trips by construction. Both satisfy
`eval(repr(s)) == s`. The obvious implementation, printing only the
constructor's arguments, produced valid Python that evaluated to a *different*
scenario with a different digest and said nothing about it.

A schedule holds up to 32 events, which is far more than any experiment in the
literature uses. The published datasets have exactly one event each: the fault,
switched on before the run.
