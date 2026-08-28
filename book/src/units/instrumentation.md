# Instrumentation

Forty-one measurements come out of the plant, and they are not all the same kind
of thing. `XMEAS(1..22)` are continuous instruments, read every step, computed at
`teprob.f:679-701` and given additive noise at `teprob.f:711-735`.
`XMEAS(23..41)` are composition analysers, read on a schedule and reporting the
composition from their *previous* sample (`teprob.f:736-761`).

Nothing in the continuous block is physics. Every one of the twenty-two lines
takes a quantity the model has already computed and converts it into the unit an
operator's instrument would read. Getting a conversion wrong changes no state and
no derivative; it changes only what the controller sees, which is worse, because
the plant then runs correctly and is controlled wrongly.

## The conversion factors

| Factor | Where | Meaning |
|---|---|---|
| `0.359` | `679`, `682`-`684`, `688` | standard cubic feet per lbmol |
| `35.3145` | the same lines, and `692`, `695`, `704`-`710` | cubic feet per cubic metre |
| `0.454` | `680`, `681`, `697` | kilograms per pound |
| `760` | `685`, `691`, `694` | mmHg per atmosphere |
| `101.325` | the same three | kPa per atmosphere |

So `FTM * 0.359 / 35.3145` is lbmol/h to standard cubic metres per hour, and
`(P - 760)/760 * 101.325` is mmHg absolute to kPa gauge. Levels are reported as
a percentage of a span, and the three vessels do not do it the same way: the
reactor and the separator carry hard-coded ranges (`teprob.f:686`, `690`) while
the stripper's span is the vessel volume `VTC` itself (`teprob.f:693`).

## `XMEAS(20)` is assigned twice

```fortran
      XMEAS(20)=CPDH*0.0003927D6
      XMEAS(20)=CPDH*0.29307D3
```

at `teprob.f:698-699`. The first is dead, and the two factors are not equal:
392.7 against 293.07, a third apart. So this is not a harmless duplicate but a
superseded conversion, and a port that took the first line would report
compressor work 34% high. That is delta D-006.

## The shutdown detector

Eight limits, checked at `teprob.f:702-710`:

| Condition | Limit | Line |
|---|---|---|
| reactor pressure high | above 3000 kPa gauge | `703` |
| reactor level high | above 24 cubic metres | `704` |
| reactor level low | below 2 cubic metres | `705` |
| reactor temperature high | above 175 C | `706` |
| separator level high | above 12 cubic metres | `707` |
| separator level low | below 1 cubic metre | `708` |
| stripper level high | above 8 cubic metres | `709` |
| stripper level low | below 1 cubic metre | `710` |

The original records only that *something* tripped, in a single integer `ISD`.
This port reports which, because "the plant tripped" without a reason is nearly
useless to a caller and the information is free.

All eight comparisons are strict, so a state exactly on a limit does not trip.
That matters for how the adversarial sampling pool is built: states placed *on*
the limits exercise the not-tripped side, and the tripping side needs states
past them. Both were built.

Two of the eight are phrased in terms of `XMEAS` rather than the underlying
quantity. `teprob.f:703` tests the *converted* reactor pressure against 3000 kPa
gauge, and `teprob.f:706` tests `XMEAS(9)`, which is `TCR` unconverted. Testing
`PTR` against an equivalent mmHg threshold instead would be arithmetically
different in the last bits.

What a trip *does* is described in [The right-hand side](../process/rhs.md): it
freezes all fifty derivatives (`teprob.f:807-811`), which is delta D-007.

## Noise and dead time

Noise is drawn by `TESUB6` (`teprob.f:1538-1546`), twelve uniform draws summed
and scaled, and is skipped entirely at `TIME = 0` and on a tripped plant
(`teprob.f:711`). Only the *continuous* noise is skipped, though. The analyser
blocks at `teprob.f:744-761` have no such guard, so a tripped plant still draws:
258 draws in a tripped evaluation against 522 in a healthy one, measured in
B-0027. A port that silenced everything on a trip would leave the generator 264
steps behind and desynchronise every later draw.

The dead time is a latch, and the order of two lines makes it:

```fortran
      XMEAS(I)=XDEL(I)
      CALL TESUB6(XNS(I),XMNS)
      XMEAS(I)=XMEAS(I)+XMNS
      XDEL(I)=XCMP(I)
```

The reported value is taken from the store *before* the store is updated.
Swapping those two lines gives an analyser with no dead time at all, which
produces entirely plausible numbers and a plant that is much easier to control
than the real one.

Three further details are about *when* rather than about *what*. The schedules
advance from their own previous value rather than from the current time (`TGAS =
TGAS + 0.1` at `teprob.f:751`), so a step arriving late does not shift the
schedule. Both `0.1` literals, at `teprob.f:741` and `751`, are single
precision, so the gas interval is 0.10000000149011612 and a step landing on
exactly 0.1 does *not* sample; `0.25` is exactly representable, so the product
analyser is unaffected, which is precisely the kind of inconsistency that has to
be read off the line rather than inferred from its neighbour. And at `TIME = 0`
the analysers are primed rather than sampled (`teprob.f:736-743`): the store and
the reported value are both set to the current composition, with no noise and no
draw.

## The 53 channels

Measurements and manipulated variables together are the 53 columns every
downstream consumer sees, in this order. The measurement names follow Downs and
Vogel's Table 4 and the manipulated ones their Table 3.

| # | `XMEAS` | # | `XMEAS` |
|---|---|---|---|
| 1 | A feed | 22 | condenser cooling water outlet |
| 2 | D feed | 23 | reactor feed, A |
| 3 | E feed | 24 | reactor feed, B |
| 4 | total feed | 25 | reactor feed, C |
| 5 | recycle flow | 26 | reactor feed, D |
| 6 | reactor feed rate | 27 | reactor feed, E |
| 7 | reactor pressure | 28 | reactor feed, F |
| 8 | reactor level | 29 | purge, A |
| 9 | reactor temperature | 30 | purge, B |
| 10 | purge rate | 31 | purge, C |
| 11 | separator temperature | 32 | purge, D |
| 12 | separator level | 33 | purge, E |
| 13 | separator pressure | 34 | purge, F |
| 14 | separator underflow | 35 | purge, G |
| 15 | stripper level | 36 | purge, H |
| 16 | stripper pressure | 37 | product, D |
| 17 | stripper underflow | 38 | product, E |
| 18 | stripper temperature | 39 | product, F |
| 19 | stripper steam flow | 40 | product, G |
| 20 | compressor work | 41 | product, H |
| 21 | reactor cooling water outlet | | |

`XMEAS(23..28)` and `XMEAS(29..36)` are the two gas analysers, on a 0.1 hour
schedule; `XMEAS(37..41)` is the product analyser, on 0.25 hours.

| # | `XMV` | # | `XMV` |
|---|---|---|---|
| 1 | D feed flow | 7 | separator underflow |
| 2 | E feed flow | 8 | stripper underflow |
| 3 | A feed flow | 9 | stripper steam |
| 4 | total feed flow | 10 | reactor cooling water flow |
| 5 | compressor recycle | 11 | condenser cooling water flow |
| 6 | purge valve | 12 | agitator speed |

`XMV(12)`, the agitator, is never written by any controller the driver calls, so
in a closed-loop run it has zero variance. That is a fact about the control
scheme rather than about the plant, and it is the reason the Tier 5 harness had
to learn to report a degenerate ensemble as `NaN` rather than as a number.
