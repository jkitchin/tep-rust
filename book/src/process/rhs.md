# The right-hand side

`TEFUNC` occupies `teprob.f:196-816`, six hundred lines that present themselves
as a derivative evaluation. They are not one. Understanding why is the single
most consequential structural decision in this port, so it comes before the unit
operations rather than after them.

## Three phases, because `TEFUNC` is not a pure function of `(t, y)`

Inside what looks like a right-hand side, the original also advances the
disturbance random walks, draws measurement noise, ticks the three sampled
analysers, and latches the valve commands. That is harmless for the fixed-step
Euler integrator the original uses, which evaluates the right-hand side exactly
once per step. It is wrong for anything else: an RK4 step would advance the
walks four times and draw four sets of noise.

The impure work does not sit in one place, which is the part that is not obvious
until you look. Some of it must happen *before* the derivative and some can only
happen *after* it. The walks are read at `teprob.f:407-416`, so they must be
advanced first. The measurement vector is assembled at `teprob.f:679-701` out of
flows the evaluation computes, so noise cannot be added until afterwards. One
impure call cannot sit on both sides.

The port therefore splits the routine in three:

| Phase | `teprob.f` | What it does |
|---|---|---|
| `advance_discrete` | 341-406, 793-804 | the `IDV` clamp, the `IDVWLK` mapping, the walk advance and spike draws, the `TIME = 0` initialisation, and the valve-command latch |
| `derivatives` | 407-710, 762-792, 805 | the entire physical model, the noise-free measurements, the shutdown test, and the fifty balances |
| `sample_measurements` | 711-761 | additive measurement noise, and the three sampled analysers with their dead time |

The pure phase hands back the signals it computed alongside the derivative, so
the post-phase does not have to re-run the model to find out what to add noise
to.

### The valve latch is hoisted, and the hoist is proved mechanically

`teprob.f:793-798` sets `IVST` from `IDV` and `799-804` latches `VCV` from
`XMV`, at the very end of the routine. That block reads only `XMV`, `VST`,
`IVST`, `IDV` and `TIME`, and nothing in `345-792` writes any of them, so moving
it into the pre-phase changes no number.

That is a claim about four hundred and fifty lines of Fortran, which is too
large to check by eye, so it is checked by machine instead: a dedicated oracle
test drives both orderings and asserts they agree. The latch shares the `DO
9020` loop with the valve derivative at `teprob.f:805`, so the port splits that
loop, sending the latch to the pre-phase and leaving `YP(I+38)` in the pure one.

## What each stage of the pure phase does

The pure phase runs in the original's order, and each block is a module in
`tepsim-core` with its own page or section in this chapter.

| `teprob.f` | Stage | Where it is documented |
|---|---|---|
| 417-472 | unpack the state into per-vessel inventories, fractions, temperatures, densities and volumes | [The plant](plant.md) |
| 473-502 | vapour-liquid equilibrium and the three vessel pressures | [reactor](../units/reactor.md), [separator](../units/separator.md), [mixing zone](../units/mixing-zone.md) |
| 503-528 | the four reactions, their rates and the heat of reaction | [The reactor](../units/reactor.md) |
| 529-564 | the stream table: compositions, molecular weights, temperatures, enthalpies | [The plant](plant.md) |
| 565-613 | valve-lagged flows, pressure-driven flows, the compressor | [The condenser and separator](../units/separator.md) |
| 614-662 | the stripper, and the reactor-inlet alias | [The stripper](../units/stripper.md) |
| 663-678 | the reactor coil, the condenser, the stripper reboiler | the three vessel pages |
| 679-710 | the twenty-two continuous measurements and the shutdown detector | [Instrumentation](../units/instrumentation.md) |
| 762-811 | the fifty balances | below |

## The fifty balances

Everything above exists to feed `teprob.f:762-811`. Four vessels, each with
eight component balances and one energy balance, plus two cooling-water wall
temperatures and twelve valve lags.

\\[
  \\frac{dn_i}{dt} = \\sum_{\\text{in}} \\dot n_i - \\sum_{\\text{out}} \\dot n_i + r_i
\\]

\\[
  \\frac{dE}{dt} = \\sum_{\\text{in}} h F - \\sum_{\\text{out}} h F + Q
\\]

The reactor is the only vessel with a reaction term. The cooling water walls
follow

\\[
  \\frac{dT_w}{dt} = \\frac{F_w \\times 500.53 \\times (T_{in} - T_w) - Q \\times 10^6 / 1.8}{H_w}
\\]

where `500.53` converts a cooling water flow to a heat capacity rate and the
factor \\(10^6 / 1.8\\) undoes the scaling the enthalpy correlations carry
(`teprob.f:789-792`). The `1.8` on those two lines is single precision and is
the only such occurrence in the file, which is why it has to be read off the
line rather than inferred from `teprob.f:1396`, `1404`, `1464` or `1471`, where
it is written `1.8D0`.

Each valve is a first-order lag toward its latched command (`teprob.f:805`):

\\[ \\frac{dv_i}{dt} = \\frac{c_i - v_i}{\\tau_i} \\]

| `YP` | Vessel |
|---|---|
| 1-8, 9 | reactor: component balances, then energy |
| 10-17, 18 | separator |
| 19-26, 27 | stripper |
| 28-35, 36 | mixing zone |
| 37, 38 | cooling water outlet temperatures |
| 39-50 | valve lags |

### A shutdown freezes the plant

`teprob.f:807-811` zeroes all fifty derivatives whenever any shutdown condition
holds. That does not stop the plant, it *freezes* it: the state stops moving,
the clock keeps running, and nothing in the original says so.

`PLAN.org` classes this as a Class C quirk, "behaviour-defining and
benchmark-relevant", so the *fix* needed a measured delta and a sign-off. Both
arrived on 2026-08-28, and a default `Scenario` now **ends the run** at the
trip. It is delta D-007. The port reports the trip and its cause either way,
rather than leaving a caller to infer a freeze from a vector of zeros.

`Scenario::faithful()` reproduces the freeze, and Tier 2 needs it to: the
adversarial sampling pool contains states that trip, and a port that did not
freeze would disagree with the oracle on all fifty components for every one of
them. So does any comparison against `d06` or `d18`, whose published files are
between 45% and 76% frozen tail.

## Determinism, and the two `libm` builds

`tepsim-core` is `no_std`, forbids `unsafe`, and contains no `f32`, no SIMD, no
rayon, no reordered reductions and no source of time or randomness outside the
model's own generator. `exp`, `pow` and `ln` come from a vendored pure-Rust
`libm` so that the answer does not depend on the host's C library.

That choice costs something measurable, and the project measures it rather than
hoping. Over the range of Antoine arguments this model reaches, the vendored
`libm` and gfortran's disagree on 9.945% of them, by exactly one ULP (LOG entry
B-0018). So a differential against the default build can only assert a
tolerance. Every Tier 2 comparison therefore runs a second time under a
`libm-system` feature, where `exp` is the one gfortran calls, and is held to
**zero ULP** there. Both runs are in the CI gate. Without the second, a
reassociation worth one or two ULP would pass silently.

Integer powers are a related trap and are not routed through `pow` at all.
gfortran expands `X**4` into multiplications, and the *shape* of the expansion
is load-bearing. Measured over 200,000 values with this project's pinned flags
(LOG entry B-0023): `(x*x)*(x*x)` matches gfortran on 200,000 of 200,000,
`((x*x)*x)*x` on 132,040, and `pow(x, 4.0)` on 99,523. So it is binary
exponentiation, squaring twice, and the two plausible alternatives are each
wrong a third to half of the time.
