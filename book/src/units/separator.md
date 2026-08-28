# The condenser and separator

The condenser cools the reactor effluent and the separator splits it into a
vapour and a liquid. The vapour leaves through two paths, the compressor recycle
and the purge; the liquid goes to the stripper. The compressor and its recycle
valve belong here too, because their whole job is deciding how much of the
separator's vapour goes back around the loop.

Source: `teprob.f:473-502` for the equilibrium, `teprob.f:585-601` for the flow
network and the compressor, `teprob.f:674-676` for the condenser duty, and
`teprob.f:773-777` for the balances.

## Vapour-liquid equilibrium

The separator's equilibrium has exactly the shape of [the reactor's](reactor.md)
and shares its code path. The vapour space is `VTS` less the liquid volume
(`teprob.f:474`); A, B and C get ideal gas partial pressures (`teprob.f:481`);
D through H get Raoult's law with an Antoine vapour pressure (`teprob.f:488-490`);
the composition is \\(y_i = p_i / P\\) (`teprob.f:495`); and the vapour holdup
comes back out of the ideal gas law at `teprob.f:498` and `501`.

The one number worth keeping in mind is that `PTS` floors around 811 mmHg over
the whole sampled domain. That matters for the purge clamp below.

## The condenser

A smooth saturating function of reactor outlet flow, approaching 0.404655 as the
flow grows (`teprob.f:674`), with the duty taken against the *stream 8*
temperature rather than the separator's own, and scaled by a disturbance drift
factor (`teprob.f:675-676`):

\\[ U A_s = 0.404655 \\left(1 - \\frac{1}{1 + (F_8/3528.73)^4}\\right) \\]

\\[ Q_s = U A_s \\, (T_{ws} - T_{st,8}) \\, (1 - 0.25 \\, d_{11}) \\]

### `**2` and `**4` are integer powers and must not go through `pow`

gfortran expands an integer exponent into multiplications rather than calling
libm, and the shape of that expansion is load-bearing. Measured over 200,000
values with this project's pinned flags:

| candidate for `X**4` | matches gfortran |
|---|---|
| `(x*x)*(x*x)` | 200,000 of 200,000 |
| `((x*x)*x)*x` | 132,040 |
| `pow(x, 4.0)` | 99,523 |

So it is binary exponentiation, squaring twice, and the two plausible
alternatives are each wrong about a third and a half of the time. `X**2` is
`x*x` on all 200,000, which is the only thing it could be.

## Three ways out

**The underflow to the stripper** is valve-lagged, and is the simple case
(`teprob.f:570`):

\\[ F_{11} = \\frac{v_7 R_7}{100} \\]

**The purge** is pressure-driven to atmosphere through valve 6
(`teprob.f:585-588`):

\\[ F_{10} = \\frac{v_6 \\times 0.151169 \\times \\sqrt{\\max(P_s - 760,\\, 0)}}{\\overline{M}_{10}} \\]

That clamp at zero is the one clamp in the whole flow network that **cannot be
reached**. `PTS` is a sum of eight partial pressures in a vessel that always
holds material, and it floors around 811 mmHg against a threshold of 760, so no
trajectory state, no random perturbation and no adversarial boundary takes that
branch. It is therefore covered by a unit test at a composition chosen to reach
it, rather than by the differential, on the principle that a branch no test
enters is indistinguishable from a branch that is wrong.

**The recycle** goes through the compressor. The operating point on the curve is
the pressure ratio between the mixing zone and the separator, clamped at both
ends (`teprob.f:589-591`), and the machine is fixed-speed with a cubic
pressure-ratio curve (`teprob.f:592-593`). The recycle valve then bleeds flow
back and the result has a floor (`teprob.f:596-599`):

\\[
  \\dot m = \\max\\!\\left(
    F_{\\max}\\left(1 + \\frac{1 - r^3}{1.197}\\right)
    - v_5 \\times 53.349 \\times \\sqrt{\\max(P_v - P_s,\\, 0)},
    \\; 10^{-3}
  \\right)
\\]

with \\(r = \\mathrm{clamp}(P_v / P_s,\\, 1,\\, 1.3)\\), \\(F_{\\max} = 280275\\)
(`teprob.f:1170`) and the ratio ceiling `CPPRMX` at 1.3 (`teprob.f:1171`). The
floor at \\(10^{-3}\\) exists so that the division at `teprob.f:600-601` cannot
blow up.

`PR**3` at `teprob.f:593` is an integer power, for the same reason `**4` is
above, and is written out as three multiplications rather than routed through
`pow`.

The compressor work appears twice: as an enthalpy bump on the recycle stream,
and as measurement 20 (`teprob.f:594-595`, `601`, `699`):

\\[
  W = \\dot m \\, (T_{cs} + 273.15) \\times 1.8 \\times 10^{-6} \\times 1.9872
      \\times \\frac{P_v - P_s}{\\overline{M}_9 P_s}
\\]

Note that the `273.15` on that line is written `273.15D0`, double precision,
unlike the one in `TESUB2` at `teprob.f:1411`. The original is not consistent
about that constant and each occurrence has to be read off its own line.

### `HST(10) = HST(9)` is a snapshot, not an alias

`teprob.f:562` copies the separator vapour enthalpy into the purge. Both streams
leave the separator vapour space, so at that moment they are the same fluid at
the same temperature and the copy is exact.

Then `teprob.f:601` adds the compressor work to `HST(9)`. The recycle gains the
work; the purge does not, because it was copied first. Reading line 562 as an
alias rather than as a copy would give the purge a share of compressor work it
never receives, and the two lines are seventy apart, so the ordering is easy to
miss. Both energy balances that read stream 9 come after line 601 and therefore
see the bumped value.

## Balances

One inlet, three outlets, and the condenser duty (`teprob.f:762-770` and
`773-777`):

\\[
  \\frac{dn_i}{dt} = \\dot n_{i,8} - \\dot n_{i,9} - \\dot n_{i,10} - \\dot n_{i,11}
\\]

\\[
  \\frac{dE}{dt} = h_8 F_8 - h_9 F_9 - h_{10} F_{10} - h_{11} F_{11} + Q_s
\\]

The condenser cooling water wall temperature is `YP(38)` (`teprob.f:791-792`).

## Variables

| Fortran | Meaning | Where |
|---|---|---|
| `VTS`, `VLS`, `VVS` | total, liquid and vapour volume | `teprob.f:1119`, `471`, `474` |
| `PPS(1:8)`, `PTS` | partial and total pressure, mmHg | `teprob.f:481`, `489`, `490` |
| `XVS`, `XLS` | vapour and liquid mole fractions | `teprob.f:495`, `452` |
| `UTVS`, `UCVS` | total and per-component vapour moles | `teprob.f:498`, `501` |
| `TCS`, `TKS` | temperature, Celsius and kelvin | `teprob.f:462-463` |
| `DLS` | liquid molar density | `teprob.f:468` |
| `UAS`, `QUS` | condenser coefficient and duty | `teprob.f:674-676` |
| `PR`, `CPPRMX` | compressor pressure ratio and its ceiling | `teprob.f:589-591`, `1171` |
| `CPFLMX` | maximum compressor flow | `teprob.f:1170` |
| `CPDH` | compressor enthalpy bump | `teprob.f:594-595` |
| `FTM(9)`, `FTM(10)`, `FTM(11)` | recycle, purge, underflow | `teprob.f:600`, `588`, `570` |
| `TWS` | cooling water outlet temperature, `YY(38)` | `teprob.f:436` |
| `YP(10..17)`, `YP(18)` | component and energy derivatives | `teprob.f:762-770`, `773-777` |

Two of the eight shutdown conditions belong to this vessel: separator liquid
volume above 12 or below 1 cubic metre (`teprob.f:707-708`).
