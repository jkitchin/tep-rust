# Introduction

## What the Tennessee Eastman Process is

The Tennessee Eastman Process is the chemical engineering profession's standard
open benchmark for process control and for fault detection. Downs and Vogel
published it in 1993 as a challenge problem: a model of an industrial chemical
plant whose components are identified only by the letters A to H, distributed
as Fortran so that competing control schemes and competing detectors could be
compared on identical dynamics rather than on identical prose descriptions.
Thirty years of published results rest on that code.

The plant makes two liquid products, G and H, from four gaseous feeds, A, C, D
and E, with an inert B and a liquid byproduct F. Four irreversible, exothermic
gas-phase reactions run in the reactor's vapour space (`teprob.f:503-528`):

```text
1:  A + C + D -> G          2:  A + C + E -> H
3:  A + E     -> F          4:  3 D       -> 2 F
```

B takes part in none of them. It arrives with the mixed feed and leaves only
through the purge, which is the entire reason the plant has a purge. Around
those reactions sit five unit operations, described one page each later in this
book: a feed mixing zone, the reactor with its cooling coil and agitator, a
condenser and vapour-liquid separator, a compressor with a recycle loop, and a
steam stripper. The whole plant is fifty integrated states (`teprob.f:24-26`),
forty-one measurements, twelve manipulated variables, and twenty programmed
disturbances (`teprob.f:340`, which loops `DO 500 I=1,20`).

Two properties make it a good benchmark and a hard one. It is open-loop
unstable, so a run with the valves held still does not merely drift, it trips:
this port measures that trip at 3.060 simulated hours on reactor pressure (LOG
entry B-0041). And its composition analysers report their *previous* sample
rather than the current one (`teprob.f:711-761`), so the measurements a
controller sees carry real dead time.

## What this port is

`tepsim` is a pure-Rust reimplementation of that Fortran, ported from the
original source rather than from any later reimplementation. The reference
material is the Braatz group release, vendored unmodified under `reference/`:
`teprob.f`, 1594 lines holding the plant model and eight utility routines, and
`temain_mod.f`, 1413 lines holding the Euler driver and the decentralised PI
control suite.

Two design commitments shape everything. The first is that the port is checked
against the Fortran continuously rather than at the end: a development-only
crate, `tepsim-oracle`, compiles the unmodified Fortran with gfortran and links
it into the test binary, so a Rust test can force both implementations into the
same state and compare them in the same process. The second is determinism. The
core crate is `no_std`, forbids `unsafe`, uses no `f32`, no SIMD, no clock and
no randomness outside the model's own generator, and answers `exp`, `pow` and
`ln` from a vendored pure-Rust `libm`, so that the same scenario produces the
same bits on x86-64, on aarch64 and in a browser.

Every place the port knowingly departs from the original is written down, with
its class, its measured effect and the test that measured it, in the [quirk and
delta register](deltas.md).

## The headline result

A complete 48-hour closed-loop run, 172,800 integrator steps at a one-second
step, is **bit-identical to the Fortran in all 41 measurements and all 12
manipulated variables**, at every one of those steps, when both are given the
same `exp` and `pow`.

That measurement is from the `LOG.org` entry for B-0041, "Phase 4 acceptance:
the driver's full 48-hour run", with gfortran 15.2.0:

| `libm`   | worst `XMEAS` | worst `XMV` (fraction of range) | within `XNS` |
|----------|---------------|---------------------------------|--------------|
| platform | 0             | 0                               | 48.000 h     |
| vendored | 1.705e-10 at `XMEAS(1)` | 1.246e-11 at `XMV(10)`| 48.000 h     |

The two rows are the two halves of one claim. Under `libm-system`, where the
transcendental functions are the ones gfortran itself calls, the port and the
original agree exactly: not to a tolerance, to the bit. Under the vendored
`libm` that actually ships, the only difference is transcendental rounding, and
after two simulated days it has not reached a tenth of a billionth of any
instrument's noise standard deviation.

The size of that transcendental difference is measured too. The first
transcendental call in the model is the Antoine vapour pressure at
`teprob.f:485` and `teprob.f:488`; over the whole range of arguments this model
reaches, the vendored `libm` and gfortran's disagree on 9.945% of them, by
exactly one ULP each time (LOG entry B-0018). Everything downstream of a
condensible partial pressure therefore carries about 1.1e-16 of relative
difference that no care in the algebra removes. Rather than accept a tolerance
and lose the sharper claim, every differential test runs twice, once against
each `libm`, and the second run is held to zero ULP.

Since B-0052 that claim attaches to the public API and not only to the test
harness: `tepsim::Simulation` reproduces the validated loop bit for bit over the
nominal scenario and `IDV(1)`, `IDV(6)`, `IDV(13)` and `IDV(20)`, on all 60
samples and all 53 channels, in both `libm` configurations.

## Run it in the browser first

[**TEP Studio**](studio/) is this simulator compiled to WebAssembly, running
entirely in your browser with no server and nothing to install. Start the plant,
switch a disturbance on, and watch the controllers fight it. Every scenario is
shareable as a URL, because the whole scenario is in the link.

It is the same code the rest of this book documents, not a demonstration model:
the page prints a determinism digest that matches the one the native build
computes, and `apps/studio/node/deployed.test.mjs` asserts that on the module
that actually ships.

## How to read this book

[Getting started](getting-started.md) is the working API and the command line
tool, and is the fastest way to a CSV of plant data. [The plant](process/plant.md)
and [the right-hand side](process/rhs.md) explain how the fifty states, the
thirteen internal streams and the six hundred lines of `TEFUNC` are organised
here. The five unit-operation pages carry the governing equations, the variable
tables and the `teprob.f` line ranges they were derived from.
[Validation](validation.md) is the ten-tier ladder and the numbers actually
measured at each rung. [Status](status.md) says what is finished and what is
not, which for a port in progress is the page to read second.

Every claim about the model in this book cites the `teprob.f` or `temain_mod.f`
line range it came from, so a sceptical reader can check it against the
vendored source. Every validation number cites the `LOG.org` iteration that
measured it.
