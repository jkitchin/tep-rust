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

- Nothing open. D-001 is decided and measured.

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
