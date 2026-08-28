# The twenty disturbances

`IDV(1..20)` is a bare integer array in the original, and the header at
`teprob.f:172-191` names each entry in prose. Nothing there connects a name to
the line that implements it, and five of the twenty are called only "Unknown".
This page makes the connection: every fault names the `teprob.f` line it acts
on, where it is injected, and what shape it has.

There are twenty, not the twenty-one of the later literature. `teprob.f:340` is
`DO 500 I=1,20`. The extra one comes from later versions of the model, and this
port follows the vendored source.

## Three shapes, and the third is not a plant disturbance at all

**Step** faults change a feed condition the moment they are switched on and hold
it. Seven of the twenty, at `teprob.f:407-414`, `567` and `568`.

**Random** faults enable a walk channel through `IDVWLK` (`teprob.f:347-358`).
Ten of the twenty, of which the last three drive *spike trains* rather than
walks.

**Sticking** faults do not touch the model. They set `IVST`
(`teprob.f:793-798`), which widens the dead band a valve command must cross
before the valve follows it. Three of the twenty.

That third kind matters more than its size. A sticking fault is not a
disturbance to the plant; it is a disturbance to the controller's authority over
the plant. In an open-loop run, where the command never moves, it does nothing
whatever, and a scenario engine that treated it as a plant fault would report an
injected disturbance with no effect and look broken. Tier 4 confirmed this from
a direction that could not have been arranged: over a four-hour run, `IDV(14)`,
`IDV(15)` and `IDV(19)` report worst errors identical to the nominal case to
every digit, because with the command held still their trajectory *is* the
nominal trajectory.

## The table

| `IDV` | Published description (`teprob.f:172-191`, verbatim) | What the source does | Shape | Line |
|---|---|---|---|---|
| 1 | A/C Feed Ratio, B Composition Constant (Stream 4) | steps the mixed feed's A fraction down by 0.03 | step | `teprob.f:407` |
| 2 | B Composition, A/C Ratio Constant (Stream 4) | steps B up by 0.005 and A down by 2.43719e-3, on two lines | step | `teprob.f:408-409` |
| 3 | D Feed Temperature (Stream 2) | steps the D feed temperature up by 5 C | step | `teprob.f:411` |
| 4 | Reactor Cooling Water Inlet Temperature | steps the reactor coolant inlet up by 5 C | step | `teprob.f:413` |
| 5 | Condenser Cooling Water Inlet Temperature | steps the condenser coolant inlet up by 5 C | step | `teprob.f:414` |
| 6 | A Feed Loss (Stream 1) | shuts the A feed off entirely, not partially | step | `teprob.f:567` |
| 7 | C Header Pressure Loss - Reduced Availability (Stream 4) | reduces the mixed feed's capacity by 20% | step | `teprob.f:568` |
| 8 | A, B, C Feed Composition (Stream 4) | enables two walk channels, on A and on B | random, channels 1 and 2 | `teprob.f:347-348` |
| 9 | D Feed Temperature (Stream 2) | enables the D feed temperature walk | random, channel 3 | `teprob.f:349` |
| 10 | C Feed Temperature (Stream 4) | enables the mixed feed temperature walk | random, channel 4 | `teprob.f:350` |
| 11 | Reactor Cooling Water Inlet Temperature | enables the reactor coolant inlet walk | random, channel 5 | `teprob.f:351` |
| 12 | Condenser Cooling Water Inlet Temperature | enables the condenser coolant inlet walk | random, channel 6 | `teprob.f:352` |
| 13 | Reaction Kinetics | enables two walks, one per rate constant of reactions 1 and 2 | random, channels 7 and 8 | `teprob.f:353-354` |
| 14 | Reactor Cooling Water Valve | sticks valve 10; touches no equation in the model | sticking | `teprob.f:793` |
| 15 | Condenser Cooling Water Valve | sticks valve 11; touches no equation in the model | sticking | `teprob.f:794` |
| 16 | Unknown | enables walk channel 9, the stripper steam valve capacity | random, channel 9 | `teprob.f:355` |
| 17 | Unknown | enables spike channel 10, the reactor coolant duty | random, spiking, channel 10 | `teprob.f:356` |
| 18 | Unknown | enables spike channel 11, the condenser coolant duty | random, spiking, channel 11 | `teprob.f:357` |
| 19 | Unknown | sticks valves 5, 7, 8 and 9; touches no equation in the model | sticking | `teprob.f:795-798` |
| 20 | Unknown | enables spike channel 12, the reactor outlet flow | random, spiking, channel 12 | `teprob.f:358` |

The channel column is not decoration. A test asserts that this table agrees with
the code that maps `IDV` flags to channel flags, that every one of the twelve
channels is driven by exactly one fault, that the spiking flag is set exactly for
channels 10 and above, and that exactly `IDV(14)`, `IDV(15)` and `IDV(19)` fail
to reach the plant. Two statements of one fact, which is the point.

## The five "Unknown" faults are not unknown

The header calls `IDV(16)` through `IDV(20)` unknown, and every paper on TEP
repeats it. The *source* is perfectly explicit about what they do; only their
physical interpretation was withheld. They enter the model at these points:

| Fault | Where it lands |
|---|---|
| `IDV(16)` | the stripper steam valve capacity, `UAC` at `teprob.f:572` |
| `IDV(17)` | the reactor coil duty, the drift factor at `teprob.f:673` |
| `IDV(18)` | the condenser duty, the drift factor at `teprob.f:676` |
| `IDV(19)` | sticks valves 5, 7, 8 and 9 |
| `IDV(20)` | the reactor outlet flow resistance, at `teprob.f:582-583` |

So `IDV(19)` is a sticking fault and the other four are not, which the shared
label hides. Three of the four are the *spike* channels, which is why they are
reported in the literature as the hardest to detect: they are intermittent
rather than sustained.

## Nine walks and three spike trains

The twelve channels are not the same kind of object, which the shared array
names hide (`teprob.f:340-406`).

**Channels 1 to 9 are random walks.** When one runs out, `teprob.f:359-371`
evaluates the old segment at its endpoint, takes the value and the slope there,
and builds a segment that continues smoothly from them. `TESUB5`
(`teprob.f:1506-1537`) chooses the next knot value and slope from the uniform
generator and fits a Hermite cubic; `TESUB8` (`teprob.f:1300-1359`) evaluates
the cubic at the current time.

**Channels 10 to 12 are spike trains**, and `teprob.f:372-396` gives them their
own rule. They alternate between two states. Dwelling: the channel sits at zero
for a randomly drawn interval, its segment being a parabola rising from zero,
and the dwell ends when the value reaches 0.1. Spiking: once the value exceeds
0.1, the channel is given a cubic through the current value and slope that lasts
exactly 0.1 hours, with coefficients that drive it hard and then back down. So a
spike channel is off, off, off, then briefly on.

**The flag scales the dwell, not the schedule.** `CDIST(I) = IDVWLK(I) / h^2` at
`teprob.f:391` is the only place the flag enters a spike channel. With the
disturbance off it is zero, the parabola is flat at zero, the channel never
reaches 0.1, and it dwells forever, drawing a fresh interval each time. With it
on, the parabola climbs and the channel eventually spikes. Either way it keeps
drawing at the same rate.

That last property is load-bearing and is easy to optimise away. `IDVWLK`
multiplies the two *endpoint* draws at `teprob.f:1529-1530` but not the
*duration* draw at `teprob.f:1528`, so an inactive walk channel still consumes
all three draws. Skipping the two endpoint draws when the flag is zero produces
*identical* segment values, because an inactive channel lands on `SZERO` with
zero slope either way, and leaves the generator two steps behind. Every
subsequent draw in the run then differs: the noise, the other channels,
everything. That mutation was implemented deliberately to check the tests have
teeth, and it was caught on the *stream position* rather than on any value. It
is the entire argument for [Tier 3](validation.md) demonstrated at Tier 1.

## Using them

```console
$ tep faults                       # the table above, from the source
$ tep run --fault 4 --hours 24 --labels
```

or, from Rust, `Scenario::fault(4)` and `Scenario::baseline().with_fault(4)`.
Faults can be combined; the flags are independent.

One caveat that is not the caller's doing. The driver switches `IDV(12)` on at
eight hours whatever the scenario asked for (`temain_mod.f:366-368`), because
that is what the original driver does and what every published dataset longer
than eight hours contains. It is delta D-011, the `Scenario` field
`driver_forces_idv12` turns it off, and the ground-truth labels record it, so a
"fault-free" 48-hour run is visibly fault-free for eight hours and then not.
