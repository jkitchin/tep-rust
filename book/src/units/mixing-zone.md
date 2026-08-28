# The mixing zone

The feed mixing zone is a single vapour volume in which the three pure feeds,
the compressor recycle and the stripper overhead combine before entering the
reactor. It is the simplest of the four vessels: one phase, a fixed volume, no
reaction, no heat duty.

Source: `teprob.f:428` and `434` for its states, `teprob.f:465-466` for its
temperature, `teprob.f:492` for its pressure, `teprob.f:576-579` for its outlet
flow, and `teprob.f:783-788` for its energy balance.

## Equations

All eight components are vapour and the volume is fixed at `VTV`, so the
equilibrium block needs only the ideal gas law applied to the mixture
(`teprob.f:492`):

\\[ P_v = \\frac{N_v R T_K}{V_{tv}} \\]

The temperature comes from the specific internal energy, not from a state
directly. The block totals the holdups, forms the mole fractions and the energy
per mole, and then solves for whatever temperature makes the mixture's specific
internal energy equal to it (`teprob.f:465-466`):

\\[ N = \\sum_i n_i, \\qquad x_i = \\frac{n_i}{N}, \\qquad e = \\frac{E}{N} \\]

That solve is Newton's method in `TESUB2` (`teprob.f:1415-1442`), warm-started
from the previous evaluation's answer, which is why `TCV` is state rather than a
derived quantity. See [The plant](../process/plant.md).

Flow out of the mixing zone is not valve-driven. It is a square-root resistance
across the pressure difference to the reactor, converted from mass to moles by
the stream's mean molecular weight (`teprob.f:576-579`):

\\[ F_6 = \\frac{1937.6 \\, \\sqrt{\\max(P_v - P_r,\\, 0)}}{\\overline{M}_6} \\]

The clamp at zero is what stops a reversed pressure gradient from producing a
`NaN` out of the square root, and it is reachable from an adversarial state
though not from the nominal trajectory.

The component and energy balances have five inlets and one outlet
(`teprob.f:762-770` and `783-788`):

\\[
  \\frac{dn_i}{dt} = \\dot n_{i,1} + \\dot n_{i,2} + \\dot n_{i,3}
  + \\dot n_{i,5} + \\dot n_{i,9} - \\dot n_{i,6}
\\]

\\[
  \\frac{dE}{dt} = h_1 F_1 + h_2 F_2 + h_3 F_3 + h_5 F_5 + h_9 F_9 - h_6 F_6
\\]

There is no `Q` term: the mixing zone is adiabatic.

## Variables

| Fortran | Meaning | Where |
|---|---|---|
| `UCVV(1:8)` | vapour component holdup, `YY(28..35)` | `teprob.f:428` |
| `ETV` | internal energy, `YY(36)` | `teprob.f:434` |
| `UTVV` | total vapour moles | `teprob.f:443-449` |
| `XVV` | vapour mole fractions | `teprob.f:450-455` |
| `ESV` | specific internal energy | `teprob.f:459` |
| `TCV` | temperature, degrees Celsius | `teprob.f:465-466` |
| `PTV` | total pressure, mmHg | `teprob.f:492` |
| `VTV` | vessel volume, 5000 cubic feet, single precision | `teprob.f:1121` |
| `FTM(6)` | outlet molar flow | `teprob.f:576-579` |
| `YP(28..35)`, `YP(36)` | the nine derivatives | `teprob.f:762-770`, `783-788` |

## Three things the source settles

**The mixed A/C feed does not pass through here.** Stream 4 goes directly to the
stripper, and it appears in the stripper's energy balance at `teprob.f:778-782`
and in the stripper's feed at `teprob.f:614-662`, never in the mixing zone's.
The three feeds that do enter are streams 1, 2 and 3, the D, E and A feeds.

**Stream 7 is an alias of stream 6, made in the stripper block.**
`teprob.f:656-661` copies flow, enthalpy, temperature, composition and component
flows from 6 to 7 wholesale. There is no mixing, no pressure drop and no heat
loss between them; stream 7 exists so that the reactor's balance at
`teprob.f:763-772` can name its own inlet.

**The pressure published as stripper pressure is this vessel's.** `teprob.f:694`
reports `PTV` as `XMEAS(16)`, which Downs and Vogel's measurement table names
stripper pressure. The model carries no separate stripper vapour space: the
stripper's overhead discharges into the mixing zone as stream 5, and `PTV` is
the pressure of that shared vapour node. The mixing zone's own outlet flow is
reported as `XMEAS(6)`, the reactor feed rate (`teprob.f:684`).
