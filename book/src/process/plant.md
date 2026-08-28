# The plant

This page is the vocabulary the rest of the process chapters use: the eight
components, the fifty states, and the thirteen internal streams. All three are
recovered from `teprob.f` rather than from the 1993 paper, because on two of the
three the paper and the source disagree.

## Components

Eight, `A` through `H`. `A`, `B` and `C` are non-condensible and are treated as
ideal gases throughout (`teprob.f:478`). `D` through `H` are condensible and get
an Antoine vapour pressure (`teprob.f:484`). `B` is the inert: it appears in
none of the four reactions, arrives with the mixed feed, and leaves only through
the purge.

`A`, `B` and `C` have no real liquid density correlation. `AD` is 1.0 with `BD`
and `CD` zero for all three (`teprob.f:973`, `983`, `993`), so they contribute a
flat, temperature-independent term. That is a placeholder keeping the mixing
rule finite rather than a fitted number, because the model never puts them in a
liquid phase in quantity.

## The fifty states

The original carries the state in a bare `YY(50)` and unpacks it by index
arithmetic inside `TEFUNC` (`teprob.f:417-440`). Recovering that mapping is the
first act of the port, and here it becomes typed structure, pinned by a test
that reads the corresponding `COMMON/TEPROC/` variables back out of the Fortran
rather than trusting a comment.

| `YY` (1-based) | Fortran | Meaning | Count |
|---|---|---|---|
| 1-3 | `UCVR(1:3)` | reactor vapour holdup, A, B, C | 3 |
| 4-8 | `UCLR(4:8)` | reactor liquid holdup, D through H | 5 |
| 9 | `ETR` | reactor internal energy | 1 |
| 10-12 | `UCVS(1:3)` | separator vapour holdup, A, B, C | 3 |
| 13-17 | `UCLS(4:8)` | separator liquid holdup, D through H | 5 |
| 18 | `ETS` | separator internal energy | 1 |
| 19-26 | `UCLC(1:8)` | stripper liquid holdup, all eight | 8 |
| 27 | `ETC` | stripper internal energy | 1 |
| 28-35 | `UCVV(1:8)` | mixing zone vapour holdup, all eight | 8 |
| 36 | `ETV` | mixing zone internal energy | 1 |
| 37 | `TWR` | reactor cooling water outlet temperature | 1 |
| 38 | `TWS` | condenser cooling water outlet temperature | 1 |
| 39-50 | `VPOS(1:12)` | valve positions, one first-order lag each | 12 |

**The eight slots do not mean the same thing in every vessel**, and this is the
part that is easy to get wrong. For the reactor and the separator the array is
split by phase: slots 1 to 3 are the *vapour* holdups of A, B and C, and slots 4
to 8 are the *liquid* holdups of D through H. `UCLR(1..3)` is set to zero at
`teprob.f:420-421` because the non-condensibles never form a liquid, and
`UCVR(4..8)` does not come from the state at all: it is derived from the
vapour-liquid equilibrium later in the same call (`teprob.f:500-501`). For the
stripper all eight slots are liquid, and for the mixing zone all eight are
vapour.

### The four temperatures are state, not derived quantities

`TESUB2` takes its temperature argument as both the initial guess and the result
(`teprob.f:1432`, `1438`), and the four call sites at `teprob.f:460-465` pass
`TCR`, `TCS`, `TCC` and `TCV` straight out of `COMMON`. Every evaluation
therefore starts its Newton solves from the previous evaluation's answers, and
since the iteration stops on a step below 1e-12 the converged value depends on
where it started.

That is not a detail. B-0015 measured the cost of getting it wrong: seeding the
solves from a different point on the nominal trajectory moves up to 21 of the 50
derivatives. A port that solved from a fixed guess would be tidier and would not
be bit-exact. The warm-start temperatures are carried explicitly here, and B-0034
found the same thing again from the other end, where a trajectory started from
the nominal literals instead of from the values `TEINIT`'s own evaluation leaves
behind is a *different trajectory* rather than a rounding of the same one.

| vessel | after `TEINIT` | nominal literal |
|---|---|---|
| reactor | 120.3999996050374 | 120.4 |
| separator | 80.1094039945582 | 80.109 |
| stripper | 65.7310297718018 | 65.731 |
| mixing zone | 86.1201119771066 | 86.120 |

Those four values are asserted against the oracle bit for bit; they are from the
`LOG.org` entry for B-0052.

## The thirteen internal streams

**The Fortran's stream indices are not the stream numbers in the paper.**
`FTM(1)` is the D feed, which Downs and Vogel call stream 2. `FTM(3)` is the A
feed, which they call stream 1. Nothing in the source says so, and every
reimplementation of TEP has to rediscover it; getting it wrong produces a plant
that runs, looks plausible, and is wired up incorrectly.

| Internal | Paper | Stream |
|---|---|---|
| 1 | 2 | D feed |
| 2 | 3 | E feed |
| 3 | 1 | A feed |
| 4 | 4 | A and C feed |
| 5 | 5 | stripper overhead vapour to the mixing zone |
| 6 | 6 | mixing zone outlet to the reactor |
| 7 | 6 | reactor inlet, an alias of 6 |
| 8 | 7 | reactor outlet to the condenser and separator |
| 9 | 8 | separator vapour through the compressor, the recycle |
| 10 | 9 | purge |
| 11 | 10 | separator liquid underflow to the stripper |
| 12 | none | stripper liquid downflow, internal only |
| 13 | 11 | product |

The mapping was established from the source, not from the paper: `teprob.f:565`
drives `FTM(1)` from valve 1 and `XMV(1)` is documented as "D Feed Flow (stream
2)"; `teprob.f:567` gates `FTM(3)` on `IDV(6)`, documented as "A Feed Loss
(Stream 1)"; `teprob.f:688` reports `FTM(10)` as `XMEAS(10)`, "Purge Rate
(stream 9)"; `teprob.f:683` reports `FTM(9)` as `XMEAS(5)`, "Recycle Flow
(stream 8)"; and so on for the remaining six.

Streams 6 and 7 are the same fluid. `teprob.f:656-661` copies flow, enthalpy,
temperature, composition and component flows from 6 to 7 wholesale, with no
mixing, no pressure drop and no heat loss. Stream 7 exists so that the reactor's
balance at `teprob.f:763-772` can name its own inlet.

## Vessel volumes

Four vessels, four fixed total volumes, all four written in the original as
single-precision literals (`teprob.f:1118-1121`):

| Vessel | Fortran | Value, cubic feet |
|---|---|---|
| reactor | `VTR` | 1300 |
| separator | `VTS` | 3500 |
| stripper | `VTC` | 156.5 |
| mixing zone | `VTV` | 5000 |

The reactor and the separator hold two phases, so their vapour space is
whatever the liquid does not occupy. The stripper is treated as liquid only and
the mixing zone as vapour only.

## Precision is a property of each literal

The line above is not pedantry. 182 of the assignments in `TEINIT` are
single-precision literals, and a literal written without a `D` suffix is stored
by gfortran as a single-precision value widened to double, which differs from
the decimal number by up to about 6e-8 relative. Since the port must reproduce
the original's arithmetic bit for bit, every constant has to be transcribed
according to the suffix on its own line.

The original is not consistent, so the precision cannot be inferred from
elsewhere in the file. `teprob.f:1411` writes `273.15` and `teprob.f:594` writes
`273.15D0`. The `1.8` at `teprob.f:790` and `792` is single while every other
occurrence in the file (`1396`, `1404`, `1464`, `1471`) is `1.8D0`. The gas
constant at `teprob.f:475` is `RG=998.9`, single, and it multiplies six of the
eight partial pressures in every vessel.

The canary the constants table was built around is `XMW(2)`, which must come out
25.399999618530273 and not 25.4.
