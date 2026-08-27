# Quirk and delta register

Every deliberate deviation from the original Fortran gets an entry here: what
the original does, what this port does instead, the class, the **measured**
effect, and the test that measures it. An entry with a description but no
number is not finished.

The classes come from `PLAN.org`, "Quirk and delta register", and they decide
how much caution an entry warrants:

| Class | Meaning | Disposition |
|---|---|---|
| **A** | No numerical effect. Dead code, naming confusion, reentrancy. | Fixed without discussion. |
| **B** | Numerically observable, semantically clearly wrong. | Fixed, with a measured delta. Expected to measure as zero under normal operation, which is itself the thing to demonstrate. |
| **C** | Behaviour-defining and benchmark-relevant. | Implemented behind a flag, never the default, and released only on explicit sign-off after a full Tier 5 and Tier 6 delta report. |

Each entry is anchored in the source by a `@delta` comment immediately above
the item it applies to, in the same spirit as the `@port` claims that
`cargo xtask provenance` collects:

```rust
// @port  teprob.f:1415-1442
// @delta D-001 class=B teprob.f:1439-1440
pub fn temperature_from_enthalpy(/* ... */) { /* ... */ }
```

`PLAN.org` calls for this page to be generated from those annotations. It is
written by hand for now; the generator is a Phase 9 item, and the annotation
convention exists from the first entry so that there is something to generate
from when it arrives.

## Open questions

Read this section before re-litigating a quirk. It is the list of things that
have been *noticed* but not yet *decided*.

- Nothing open. D-001, D-002 and D-003 are decided and measured.

---

## D-001 — `TESUB2` reports success after failing to converge

**Class B.** `teprob.f:1439-1440`, in `TESUB2` (`teprob.f:1415-1442`).

### What the original does

The Newton loop is written so that the convergence test is the loop-terminal
statement:

```fortran
      TIN=T
      DO 250 J=1,100
      ...
      T=T+DT
 250  IF(DABS(DT).LT.1.D-12)GO TO 300
      T=TIN
 300  RETURN
```

On convergence, `GO TO 300` leaves the loop with the solved temperature. On
failure the loop simply runs out, control falls through to `T=TIN`, and the
routine **restores the caller's original guess and returns exactly as it does
on success**. There is no error code, no status flag and no output argument
that distinguishes the two. A caller receiving a temperature cannot tell
whether it solved the problem or is holding its own input back.

The consequence is not confined to one number. `TESUB2` is what converts the
reactor, separator and stripper energy states into temperatures
(`teprob.f:460-465`), so a silently abandoned solve seeds every downstream
pressure, flow and heat-transfer term for that step with a stale temperature,
and the run continues as though nothing happened.

### What this port does

[`temperature_from_enthalpy`] returns `Result<f64, TemperatureError>`.
Convergence returns `Ok`; exhausting the hundred iterations returns
`TemperatureError::DidNotConverge`, carrying the guess, the last iterate and
the final step, so a caller can report or recover. The iteration itself is
unchanged: same evaluation order, same step, same criterion, same cap.

### Measured effect

**Zero.** Across the full Tier 1 sweep the non-convergence path never fires,
so the two behaviours never differ on the physical domain.

| Basis | Start | Cases | Abandoned | Rust vs Fortran |
|---|---|---|---|---|
| `ITY=0` | warm | 9,987,490 | 0 | 0 ULP |
| `ITY=0` | cold | 9,987,490 | 0 | 0 ULP |
| `ITY=1` | warm | 9,987,490 | 0 | 0 ULP |
| `ITY=1` | cold | 9,987,490 | 0 | 0 ULP |
| `ITY=2` | warm | 9,987,490 | 0 | 0 ULP |
| `ITY=2` | cold | 9,987,490 | 0 | 0 ULP |

59,924,940 solves, none abandoned, every returned temperature bit-identical to
the Fortran's. "Warm" starts Newton from the temperature the target enthalpy
was built from, which is what every call site in the plant does; "cold" starts
it from the opposite end of the 0-175 °C range, which nothing in the plant
does and which is there to make the iteration work for its answer.

Round-trip accuracy, as a diagnostic rather than a gate: the warm start
recovers the original temperature exactly, and the cold start lands within
3.7e-13 °C of it. Newton is quadratic on this problem, so the 1e-12 criterion
on the *step* certifies an error far below itself.

### Disposition

Adopted as the default. The fix is free: it changes no number the model can
reach, and it converts a failure mode that is invisible by construction into
one the type system forces a caller to handle.

The measurement is not a one-off. It is the assertion
`total_abandoned == 0` in the Tier 1 test below, so if a future change to the
sweep, the constants or the iteration ever reaches the non-convergence path,
it fails there rather than going unnoticed.

### Measured by

`crates/tepsim-oracle/tests/tier1_temperature.rs`, run at full volume by
`cargo xtask validate --tiers 1`. Two tests: one sweeps and counts, the other
demonstrates the divergence directly on an unreachable target, asserting that
the Fortran hands back the guess verbatim while the port reports an error.

[`temperature_from_enthalpy`]: https://docs.rs/tepsim-core/latest/tepsim_core/thermo/fn.temperature_from_enthalpy.html

---

## D-002 — `R1F` and `R2F` name two unrelated quantities

**Class A.** `teprob.f:503-511`, in `TEFUNC`.

### What the original does

`R1F` and `R2F` are set at `teprob.f:415-416` from `TESUB8(7, TIME)` and
`TESUB8(8, TIME)`: the IDV(13) slow-drift multipliers on the reaction kinetics.
They are used in that meaning at `503-504`.

Five lines later they are reassigned in place:

```fortran
      RR(1)=DEXP(31.5859536-40000.0/1.987/TKR)*R1F     ! drift factor
      RR(2)=DEXP(3.00094014-20000.0/1.987/TKR)*R2F     ! drift factor
      ...
      R1F=PPR(1)**1.1544                               ! now a pressure power
      R2F=PPR(3)**0.3735                               ! now a pressure power
      RR(1)=RR(1)*R1F*R2F*PPR(4)
      RR(2)=RR(2)*R1F*R2F*PPR(5)
```

The two roles share nothing but the storage. `RR(1)` at line 510 carries the
drift factor once, from line 503, and the pressure powers once, from line 508.

### Why it matters

There is no numerical effect, but there is a large *reading* effect. Taking
`R1F` at line 510 to still be the drift factor gives

```text
r1 = f1 · e^(a1 - E1/T) · f1 · pC^0.3735 · pD · Vv
```

which is a coherent-looking rate law, second order in the disturbance and
missing the pressure order on A entirely. Nothing in the source contradicts it
locally, and a port written that way is self-consistent: it reproduces its own
answer on every run.

Measured: reading it that way moves `RR` by **100% relative** on the
adversarial pool, so the differential does catch it. But only because there is
a differential. This is the kind of misreading that would have propagated
through every later item in a port without one.

### What this port does

Gives the two roles separate names. The drift factors arrive as
`kinetics::ReactionDrift`, and the pressure powers are locals called
`order_a` and `order_c`.

### Measured effect

None. The arithmetic is identical; only the names differ.

The visible consequence is for harnesses, not for the model: after `TEFUNC`
returns, `COMMON`'s `R1F` holds a pressure power. Tier 2 therefore fetches the
drift from `TESUB8` directly, which is sound because nothing after
`teprob.f:406` writes the walk state.

### Measured by

`crates/tepsim-oracle/tests/tier2_kinetics.rs`. Under
`--features oracle,libm-system` the whole range is bit-identical to the
Fortran, which is what confirms the reading: a wrong-but-consistent reading
cannot be bit-identical to a right one.

---

## D-003 — `CRXR(2)` is read but never assigned

**Class A.** `teprob.f:521-527` and `teprob.f:763`.

### What the original does

Seven of the eight net-production slots are written:

```fortran
      CRXR(1)=-RR(1)-RR(2)-RR(3)
      CRXR(3)=-RR(1)-RR(2)
      CRXR(4)=-RR(1)-1.5D0*RR(4)
      CRXR(5)=-RR(2)-RR(3)
      CRXR(6)=RR(3)+RR(4)
      CRXR(7)=RR(1)
      CRXR(8)=RR(2)
```

`CRXR(2)` is not among them. It is read anyway, in the reactor component
balance at `teprob.f:763`:

```fortran
      YP(I)=FCM(I,7)-FCM(I,8)+CRXR(I)
```

for `I = 1..8`.

### Why it works

`CRXR` lives in `COMMON/TEPROC/`, which is static storage and therefore
zero-initialised, and nothing anywhere in the file ever writes slot 2. So B's
net production is zero for the life of the process.

That is the correct physics: B is inert and takes part in none of the four
reactions. But it is correct *by static initialisation*, not by statement. The
guarantee comes from the linker, not from the model.

### What this port does

States it. `Kinetics::production` is built from an explicit array of zeros and
B is simply never written, which reproduces the value while making the reason
local.

### Measured effect

None, and that is asserted rather than assumed: the oracle is required to
report exactly `0.0` for `CRXR(2)` on every sampled state. If it ever did not,
the benign reading here would be wrong and the reactor's B balance would be
picking up whatever was last left in that word.

### Measured by

`the_inert_has_no_net_production_in_either_implementation` in
`crates/tepsim-oracle/tests/tier2_kinetics.rs`. 200 nominal states, all exactly
zero.
