# The stripper

Steam strips light components out of the separator liquid. What leaves overhead
rejoins the mixing zone as stream 5; what does not becomes the stripper's own
liquid, stream 12, which leaves as the product, stream 13. The vessel is treated
as liquid only: it has no vapour space of its own in the model.

Source: `teprob.f:614-662` for the column, `teprob.f:677-678` for the reboiler,
and `teprob.f:778-782` for the energy balance.

## Equations

The feed is the mixed A/C feed plus the separator underflow (`teprob.f:635-639`).
Note that stream 4 arrives here directly rather than through the mixing zone:

\\[ f_i = \\dot n_{i,4} + \\dot n_{i,11} \\]

A vapour-to-liquid ratio sets how hard the column strips, scaled by a
temperature factor (`teprob.f:622`):

\\[ \\Lambda = \\frac{F_4}{F_{11}} \\, \\tau(T_c) \\]

Each condensible then strips according to a Langmuir-shaped saturating function
of that ratio (`teprob.f:623-627`):

\\[ s_i = \\frac{k_i \\Lambda}{1 + k_i \\Lambda}, \\qquad i \\in \\{D \\ldots H\\} \\]

| species | \\(k_i\\) |
|---|---|
| D | 8.5010 |
| E | 11.402 |
| F | 11.795 |
| G | 0.0480 |
| H | 0.0242 |

and the split is simply (`teprob.f:643-644`)

\\[ \\dot n_{i,5} = s_i f_i, \\qquad \\dot n_{i,12} = f_i - \\dot n_{i,5} \\]

Both product streams leave at the stripper's own temperature (`teprob.f:652-653`),
and their enthalpies are taken on different bases: stream 5 with `ITY = 1`, the
vapour basis, and stream 12 with `ITY = 0`, the liquid one (`teprob.f:654-655`).

## The temperature factor has a pole at 177 C

\\[
  \\tau(T) = \\begin{cases}
    T - 120.262 & T > 170 \\\\
    0.1 & T < 5.292 \\\\
    \\dfrac{363.744}{177 - T} - 2.22579488 & \\text{otherwise}
  \\end{cases}
\\]

from `teprob.f:615-621`. The middle branch diverges at 177 C, which is inside
the range the two outer branches leave for it only if `TCC` exceeds 170, and it
does not: the `T > 170` branch takes over first. So the pole is unreachable by
seven degrees, and the two branches are continuous to within 0.1% at 170.

The adversarial state catalogue built for Tier 2 places a state at 176 C anyway,
to sit near the pole and confirm that it stays on the linear branch. That state
is coverage of the *guard*, not of the pole.

## `FTM(11) > 0.1` switches the whole block

Below that threshold the column is not really running, and `teprob.f:629-633`
substitutes five fixed stripping factors (0.9999, 0.999, 0.999, 0.99, 0.98)
rather than evaluating the correlation. The reason is visible in the arithmetic:
\\(\\Lambda = F_4 / F_{11}\\) diverges as \\(F_{11} \\to 0\\).

Both sides are covered by the adversarial catalogue, which places a state
exactly on `FTM(11) = 0.1`. Since the test at `teprob.f:614` is `.GT.`, that
state takes the fixed-factor branch.

## `SFR(1..3)` are never recomputed

`teprob.f:623-627` and `629-633` both write slots 4 through 8 only. Slots 1, 2
and 3 are set once in `TEINIT` (`teprob.f:1126-1128`) and are read at
`teprob.f:643` on every evaluation, so A, B and C strip at a fixed 99.5%, 99.1%
and 99.0% no matter what the column is doing.

That is the intended physics rather than an oversight: the non-condensibles are
gases, they leave overhead essentially completely, and no temperature or flow
ratio in the plant's range would change that. It is worth stating because the
loop at `teprob.f:643` runs `I=1,8` and looks as though all eight factors come
from the branch above it.

## The reboiler

Steam is at 100 C, so above that there is nothing to transfer and the original
sets the duty to zero rather than letting it go negative (`teprob.f:677-678`):

\\[ Q_c = \\begin{cases} U A_c \\, (100 - T_c) & T_c < 100 \\\\ 0 & \\text{otherwise} \\end{cases} \\]

The coefficient `UAC` is not computed here. It is a valve-lagged capacity with a
disturbance drift factor, set at `teprob.f:572`:

\\[ U A_c = \\frac{v_9 R_9 (1 + d_9)}{100} \\]

which is the point at which `IDV(16)`, published as "Unknown", enters the model.

The nominal trajectory sits near 65 C, so the *cutoff* is the branch at risk of
never being exercised. B-0021 measured 300 of 300 nominal states below 100.

## The reactor inlet is an alias, and this is where it is made

`teprob.f:656-661` copies flow, enthalpy, temperature, composition and component
flows from stream 6 to stream 7 wholesale. There is no mixing, no pressure drop
and no heat loss between them: stream 7 exists so that the reactor's balance at
`teprob.f:763-772` can name its own inlet. It lives in the stripper block for no
reason other than that is where the original put it.

## Balances, and an asymmetry between them

The component balances are a straight pass-through of the two streams the column
produced (`teprob.f:762-770`), while the energy balance names the column's
*inputs* instead (`teprob.f:778-782`):

\\[
  \\frac{dn_i}{dt} = \\dot n_{i,12} - \\dot n_{i,13}
\\]

\\[
  \\frac{dE}{dt} = h_4 F_4 + h_{11} F_{11} - h_5 F_5 - h_{13} F_{13} + Q_c
\\]

That is what the source does, and the two forms describe the same vessel: the
component split at `teprob.f:643-644` conserves moles by construction, so
streams 4 and 11 in, less stream 5 out, is stream 12.

## Variables

| Fortran | Meaning | Where |
|---|---|---|
| `VTC`, `VLC` | vessel and liquid volume | `teprob.f:1120`, `472` |
| `UCLC(1:8)` | liquid component holdup, `YY(19..26)` | `teprob.f:427` |
| `ETC` | internal energy, `YY(27)` | `teprob.f:433` |
| `XLC` | liquid mole fractions | `teprob.f:453` |
| `TCC` | temperature, degrees Celsius | `teprob.f:464` |
| `DLC` | liquid molar density | `teprob.f:469` |
| `TMPFAC` | temperature scaling | `teprob.f:615-621` |
| `VOVRL` | vapour-to-liquid ratio | `teprob.f:622` |
| `SFR(1:8)` | fraction of each species stripped | `teprob.f:623-633`, `1126-1128` |
| `FIN(1:8)` | combined feed to the column | `teprob.f:635-639` |
| `UAC`, `QUC` | reboiler coefficient and duty | `teprob.f:572`, `677-678` |
| `YP(19..26)`, `YP(27)` | component and energy derivatives | `teprob.f:762-770`, `778-782` |

Two of the eight shutdown conditions belong to this vessel: stripper liquid
volume above 8 or below 1 cubic metre (`teprob.f:709-710`). Its level
measurement is also the odd one out among the three, because its span is the
vessel volume `VTC` itself rather than a separately hard-coded range
(`teprob.f:693`).
