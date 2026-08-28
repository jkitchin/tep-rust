# The reactor

The reactor is a two-phase vessel with an internal cooling coil and an
agitator. Four exothermic gas-phase reactions run in its vapour space, its
liquid holdup sets how much of the coil is wetted, and it is the only vessel in
the plant with a reaction term in its balances.

Source: `teprob.f:473-502` for the vapour-liquid equilibrium, `teprob.f:503-528`
for the kinetics, `teprob.f:663-673` for the coil, and `teprob.f:762-772` for
the balances.

## Vapour-liquid equilibrium

The vapour space is what the liquid does not occupy (`teprob.f:473`):

\\[ V_{vr} = V_{tr} - V_{lr} \\]

A, B and C are non-condensible and are treated as ideal gases, so their partial
pressures come from the holdup directly (`teprob.f:478-483`):

\\[ p_i = \\frac{n_i R T_K}{V_{vr}}, \\qquad i \\in \\{A, B, C\\} \\]

D through H are condensible, so their partial pressures come from Raoult's law
with an Antoine vapour pressure in degrees Celsius (`teprob.f:484-491`):

\\[
  p_i = x_i \\exp\\!\\left(A_i + \\frac{B_i}{T_c + C_i}\\right),
  \\qquad i \\in \\{D \\ldots H\\}
\\]

The total is the sum of all eight, the vapour composition is \\(y_i = p_i / P\\)
(`teprob.f:493-496`), and the vapour holdup follows from the ideal gas law
applied to the mixture (`teprob.f:497` and `500`):

\\[ N_v = \\frac{P V_{vr}}{R T_K}, \\qquad n_i = N_v y_i \\]

### `UCVR` is both an input and an output

`teprob.f:418` fills `UCVR(1..3)` from the state and `teprob.f:500` fills
`UCVR(4..8)` from the equilibrium computed here. It is one Fortran array written
by two different mechanisms, and the halves are not interchangeable: the
non-condensibles are integrated, the condensibles are derived. Read in that
order it also explains why the A, B and C partial pressures are computed first,
at `teprob.f:478-483`: they need `UCVR` before it is overwritten.

### This is where bit equality with gfortran ends

`DEXP` at `teprob.f:485` and `teprob.f:488` is the model's first transcendental
call. The port answers it from the vendored pure-Rust `libm` rather than the
platform's, for the determinism reason set out in [The right-hand
side](../process/rhs.md). Measured over the whole Antoine range this model
reaches, the two disagree on 9.945% of arguments, by exactly one ULP.

Everything downstream of a condensible partial pressure therefore carries about
1.1e-16 of relative difference from the Fortran that no care in the algebra
removes. The bit-exactness claim does not disappear, it moves: under the
`libm-system` feature the transcendental is the one gfortran calls, and the port
is bit-identical again, so the algebra is still held to zero ULP rather than to
a tolerance.

## Kinetics

```text
1:  A + C + D -> G          2:  A + C + E -> H
3:  A + E     -> F          4:  3 D       -> 2 F
```

Rates 1 and 2 are Arrhenius in reactor temperature with fractional pressure
orders on A and C, multiplied by a disturbance drift factor
(`teprob.f:503-504`, `508-511`):

\\[
  r_1 = f_1 \\, e^{\\,a_1 - E_1/T_K} \\; p_A^{1.1544} \\, p_C^{0.3735} \\, p_D \\, V_{vr}
\\]

\\[
  r_2 = f_2 \\, e^{\\,a_2 - E_2/T_K} \\; p_A^{1.1544} \\, p_C^{0.3735} \\, p_E \\, V_{vr}
\\]

Both are guarded: `teprob.f:507` requires \\(p_A > 0\\) and \\(p_C > 0\\), and
sets \\(r_1 = r_2 = 0\\) otherwise. Rates 3 and 4 are first order in each
reactant, and rate 4 shares rate 3's exponential rather than having one of its
own (`teprob.f:505-506`, `516-517`):

\\[
  r_3 = e^{\\,a_3 - E_3/T_K} \\, p_A \\, p_E \\, V_{vr},
  \\qquad
  r_4 = 0.767488334 \\; e^{\\,a_3 - E_3/T_K} \\, p_A \\, p_D \\, V_{vr}
\\]

All four are multiplied by the vapour volume at `teprob.f:518-520`, so a rate is
an extent in moles per hour rather than a volumetric rate. Net production per
species follows the stoichiometry (`teprob.f:521-527`), and the heat release
comes from reactions 1 and 2 only (`teprob.f:528`):

\\[ Q_{rxn} = r_1 h_1 + r_2 h_2 \\]

with \\(h_1 = 0.06899381054\\) and \\(h_2 = 0.05\\) (`teprob.f:1122-1123`).
Reactions 3 and 4 contribute no heat: the original simply does not include them
in `RH`, and `HTR(3)` is declared but never assigned or read.

### `R1F` and `R2F` are two different quantities under one name

They arrive from `TESUB8(7)` and `TESUB8(8)` at `teprob.f:415-416` as the
`IDV(13)` kinetics-drift multipliers, are consumed at `teprob.f:503-504`, and
are then **reassigned in place** at `teprob.f:508-509` to hold the fractional
pressure powers. The two meanings share nothing but the storage.

Reading `teprob.f:510` as though `R1F` were still the drift factor gives a
plausible and completely wrong rate law, so the port gives the two roles
separate names. That is delta D-002: no numerical effect, and the only defence
against a misreading no test would catch, because a wrong-but-consistent reading
still reproduces itself.

### `CRXR(2)` is never assigned

Seven of the eight slots are written at `teprob.f:521-527`. `CRXR(2)`, the
inert, is not, and it is read anyway at `teprob.f:763`. It works because
`COMMON` is zero-initialised and nothing ever writes it, so B's net production is
zero by static initialisation rather than by statement. That is delta D-003,
class A: the value is right, the mechanism is an accident. Here the slot is
explicitly zero and a test asserts the oracle agrees.

### Precision hazards in this range

Every literal except `0.767488334D0` and `1.5D0` is single precision: the three
pre-exponentials, the three activation energies, the gas constant `1.987`, and
both fractional exponents.

Worse, `40000.0/1.987` at `teprob.f:503` is a quotient of *two*
single-precision literals, so Fortran evaluates the division itself in single
precision. Widening the operands first and dividing in double is wrong by 4e-9
relative, inside a `DEXP` argument.

## The cooling coil

The coil's effective area ramps with liquid level, because the coil is only
wetted over part of its height. `VLR/7.8` is the level as a percentage, and the
ramp is piecewise linear between 10% and 50% (`teprob.f:663-669`):

\\[
  \\lambda = \\begin{cases}
    1 & \\ell > 50 \\\\
    0 & \\ell < 10 \\\\
    0.025\\,\\ell - 0.25 & \\text{otherwise}
  \\end{cases}
  \\qquad \\ell = V_{lr} / 7.8
\\]

The overall coefficient is quadratic in agitator speed (`teprob.f:670-671`), and
the duty is the coefficient times the driving temperature difference, scaled by
a disturbance drift factor (`teprob.f:672-673`):

\\[ U A_r = \\lambda \\left(-0.5\\,\\omega^2 + 2.75\\,\\omega - 2.5\\right) \\times 855490 \\times 10^{-6} \\]

\\[ Q_r = U A_r \\, (T_w - T_c) \\, (1 - 0.35 \\, d_{10}) \\]

That parabola peaks at \\(\\omega = 2.75\\), which is *above* the agitator's whole
range: `teprob.f:575` puts it between 1.5 and 2.5. So the coefficient rises
monotonically with speed everywhere the plant can go, from 0.5 to 1.25 times the
scale, and the falling half of the parabola is unreachable. Its roots are at
1.149 and 4.351, both outside that range too, so the coefficient never reaches
zero from the agitator alone. The model is a fit, not a mechanism.

### The ramp does not quite meet its flat sections

`0.025` at `teprob.f:668` is single precision, so it is stored as
0.02500000037252903 and the ramp misses both of its endpoints:

| level | ramp gives | flat section gives | gap |
|---|---|---|---|
| 10 | 3.725290298461914e-9 | 0 | 3.7e-9 |
| 50 | 1.0000000186264515 | 1 | 1.9e-8 |

Both breakpoint comparisons are strict, so a level of exactly 10 or exactly 50
takes the *ramp*, and `UARLEV` is discontinuous by those amounts as the level
crosses either one. This is faithful reproduction rather than a delta: it is
what the original computes, the gaps are eight orders below the quantity itself,
and the coefficient they scale is an empirical fit in the first place. It is
written down because "the ramp meets the flat sections" is the obvious
assumption, it is false, and a test asserting it would fail for a reason that
looks like a porting error.

## Balances

Inlet is stream 7, outlet is stream 8, and the reaction term is the only one of
its kind in the plant (`teprob.f:762-772`):

\\[
  \\frac{dn_i}{dt} = \\dot n_{i,7} - \\dot n_{i,8} + r_i,
  \\qquad
  \\frac{dE}{dt} = h_7 F_7 - h_8 F_8 + Q_{rxn} + Q_r
\\]

`QUR` is heat *removed*, so it enters positive here because \\(U A_r (T_w -
T_c)\\) is already negative when the coil is cooling. The reactor's cooling
water wall temperature is `YP(37)` (`teprob.f:789-790`).

## Variables

| Fortran | Meaning | Where |
|---|---|---|
| `VTR`, `VLR`, `VVR` | total, liquid and vapour volume | `teprob.f:1118`, `470`, `473` |
| `PPR(1:8)`, `PTR` | partial and total pressure, mmHg | `teprob.f:479`, `486`, `487` |
| `XVR`, `XLR` | vapour and liquid mole fractions | `teprob.f:494`, `451` |
| `UTVR`, `UCVR` | total and per-component vapour moles | `teprob.f:497`, `500` |
| `TCR`, `TKR` | temperature, Celsius and kelvin | `teprob.f:460-461` |
| `RR(1:4)` | extent of each reaction | `teprob.f:503-520` |
| `CRXR(1:8)` | net production per species | `teprob.f:521-527` |
| `RH` | heat of reaction | `teprob.f:528` |
| `HTR(1:2)` | heats of reactions 1 and 2 | `teprob.f:1122-1123` |
| `AGSP` | agitator speed, fraction of nominal | `teprob.f:575` |
| `UARLEV`, `UAR`, `QUR` | wetted fraction, coefficient, duty | `teprob.f:663-673` |
| `TWR` | cooling water outlet temperature, `YY(37)` | `teprob.f:435` |
| `YP(1..8)`, `YP(9)` | component and energy derivatives | `teprob.f:762-772` |

Three of the eight shutdown conditions belong to this vessel: reactor pressure
above 3000 kPa gauge (`teprob.f:703`), liquid volume outside 2 to 24 cubic
metres (`teprob.f:704-705`), and temperature above 175 degrees Celsius
(`teprob.f:706`). See [Instrumentation](instrumentation.md).
