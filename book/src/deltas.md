# Quirk and delta register

> **The two halves of this register are cross-checked.** `cargo xtask deltas`
> collects every `@delta` marker in the source, matches it against the `## D-0NN`
> headings below, and fails if an entry has no marker, if a marker has no entry,
> or if the two disagree about the class. The collected table is the generated
> [delta marker index](validation/delta-index.md). The prose below is still
> written by hand, which is what the closing paragraph of the next section is
> about; the cross-check and the index are not.

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

- **Nothing is open.** D-001 through D-011 are all decided and measured. The
  last three decisions were taken on 2026-08-28: D-007 (B-0065), D-011
  (B-0066), and what the Tier 5 battery does about a stuck valve (B-0067).
- D-007's *measure* question, which was open here, was settled on 2026-08-27:
  Tier 2 measures error against the scale of the terms entering a balance, not
  against the balance's result, because a balance is inflow minus outflow and
  cancels. See B-0026a and the module documentation of
  `crates/tepsim-oracle/tests/tier2_balances.rs`, which carries the numbers.

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

---

## D-004 — the mixed feed carries 1e-10 lbmol/h that nothing accounts for

**Class B.** `teprob.f:568-569`, in `TEFUNC`.

### What the original does

Every other valve-lagged flow is a clean proportionality:

```fortran
      FTM(1)=VPOS(1)*VRNG(1)/100.0
```

The mixed A and C feed is not:

```fortran
      FTM(4)=VPOS(4)*(1.D0-IDV(7)*0.2D0)
     .*VRNG(4)/100.0+1.D-10
```

The trailing `+1.D-10` has no physical meaning. It is there so that `FTM(4)`
is never exactly zero.

### Why it is there

`teprob.f:606` computes `FCM(I,4)=XST(I,4)*FTM(4)`, and the mixing zone
balance at `teprob.f:783-788` sums those component flows. None of that divides
by `FTM(4)`, so the guard is not protecting a division in this range.

It protects the *composition* measurement path. With the valve shut and no
epsilon, stream 4 contributes exactly nothing, and a controller reading a
composition ratio off it would see 0/0. The epsilon keeps the stream infinitely
dilute rather than absent.

### Why it is Class B and not Class A

Because it is numerically observable, permanently. The plant receives 1e-10
lbmol/h of A, B and C that no feed valve delivers and no mass balance
accounts for, for the entire run. Under normal operation `FTM(4)` is around
9.35 lbmol/h, so the addend is 1.1e-11 relative and utterly negligible; with
the valve shut it is the entire flow.

That is the shape of a Class B quirk exactly: wrong in principle, invisible in
practice, and the thing to demonstrate is that removing it measures as zero
rather than assuming so.

### What this port does

Reproduces it, as [`tepsim_core::flows::FEED_FLOW_EPSILON`], on the faithful
path. The eventual fix belongs in Phase 6 behind a flag, with a Tier 5 delta
report, like every other Class B item.

### Measured effect

Not yet. The delta is measured when the fix lands (Phase 6, Tier 10). What is
pinned now is the behaviour it produces: with valve 4 shut, `FTM(4)` is exactly
1e-10 and not zero.

### Measured by

`the_mixed_feed_never_reaches_exactly_zero` in
`crates/tepsim-core/src/flows.rs`, and the Tier 2 differential in
`crates/tepsim-oracle/tests/tier2_flows.rs`, which is bit-identical to the
Fortran under `libm-system` and would not be if the addend were dropped.

---

## D-005 — `SFR(4..8)`'s initial values are dead

**Class A.** `teprob.f:1129-1133`, set in `TEINIT`.

### What the original does

`TEINIT` sets all eight stripping factors:

```fortran
      SFR(1)=0.99500
      SFR(2)=0.99100
      SFR(3)=0.99000
      SFR(4)=0.91600
      SFR(5)=0.93600
      SFR(6)=0.93800
      SFR(7)=5.80000D-02
      SFR(8)=3.01000D-02
```

The first three are load-bearing: nothing ever writes them again, and
`teprob.f:643` reads all eight on every evaluation, so A, B and C strip at
those fixed fractions for the life of the run.

The last five are dead. `teprob.f:614-634` writes `SFR(4)` through `SFR(8)`
unconditionally, through whichever of its two branches it takes, before
`teprob.f:643` reads them. So the values on lines 1129-1133 are overwritten on
the first evaluation and never observed.

### Why it is worth an entry

Because the block looks uniform and is not. A reader checking the eight
`TEINIT` lines against the eight slots would reasonably conclude that all
eight are initial conditions, and a port that omitted lines 1129-1133 as dead
would be correct while a port that omitted 1126-1128 would be silently wrong
in the third decimal place of every product composition.

The split is also the only thing that explains why the loop at `teprob.f:643`
runs `I=1,8` when the branch above it writes only five slots.

Note the precision changes across the boundary too: `SFR(1..6)` are single
precision and `SFR(7..8)` are written `5.80000D-02` and `3.01000D-02`, in
double. Since the last five are dead, that inconsistency has no effect, which
is itself worth knowing before someone spends time on it.

### What this port does

Carries the three live values as
[`tepsim_core::stripper::NON_CONDENSIBLE_STRIPPING`] and does not carry the
five dead ones, with the reasoning stated where the constant is defined.

### Measured effect

None. The five values are unobservable.

The three live ones are asserted rather than assumed: the oracle must report
exactly the `TEINIT` constants for `SFR(1..3)` on every sampled state. If it
ever did not, they would be recomputed somewhere the port has not found.

### Measured by

`the_non_condensible_factors_never_move_in_the_fortran_either` in
`crates/tepsim-oracle/tests/tier2_stripper.rs`. 300 nominal states and all 21
adversarial boundaries, bit-identical.

---

## D-006 — `XMEAS(20)` is assigned twice, with different factors

**Class A.** `teprob.f:698-699`, in `TEFUNC`.

### What the original does

```fortran
      XMEAS(20)=CPDH*0.0003927D6
      XMEAS(20)=CPDH*0.29307D3
```

The first assignment is dead. The second overwrites it on the next line,
before anything reads the slot.

### Why it is not a harmless duplicate

The two factors are not the same number. `0.0003927D6` is 392.7 and
`0.29307D3` is 293.07: a ratio of 1.34. So this is a *superseded* conversion
rather than a repeated one, and a port that transcribed the first line would
report compressor work 34% high on every sample, forever, while every other
measurement stayed correct.

Both are plausible as unit conversions, which is what makes it a trap.
Compressor duty in the model's internal units times 293.07 gives kilowatts,
and `XMEAS(20)` is documented as "Compressor Work (kW)". The dead factor
appears to be an earlier attempt at the same conversion.

### What this port does

Takes the second, which is what the original computes, and states in
[`tepsim_core::measurements`] that the first is dead and why the difference
matters.

### Measured effect

None on the model: the value is overwritten before use, so the derivative and
every state are untouched. The effect would be entirely on a controller reading
measurement 20, and the delta against the correct value would be 34%.

### Measured by

`the_compressor_work_uses_the_second_conversion_factor` in
`crates/tepsim-core/src/measurements.rs`, which also asserts the two factors
are far enough apart that taking the dead one would be visible. The Tier 2
differential in `crates/tepsim-oracle/tests/tier2_measurements.rs` covers it
against the oracle, and is bit-identical under `libm-system`.

---

## D-007 — a shutdown freezes the plant instead of stopping it

**Class C.** `teprob.f:807-811`, in `TEFUNC`. **Fixed by default since the
sign-off of 2026-08-28 (B-0065); reproduced by `QuirkFixes::faithful`.**

### What the original does

```fortran
      IF(ISD.NE.0)THEN
      DO 9030 I=1,NN
      YP(I)=0.0
 9030 CONTINUE
      ENDIF
```

When any of the eight shutdown conditions holds, all fifty derivatives become
zero. The state stops moving, the clock keeps running, and the caller is told
nothing: `ISD` is a local, and the returned vector is indistinguishable from a
plant at perfect steady state.

### Why it is Class C rather than B

Because published results depend on it. Every `d00`-`d21` dataset was generated
by a driver that kept integrating through a trip, and a run that ended instead
would produce a shorter, different file. Changing this changes
benchmark comparability, which is the definition `PLAN.org` gives for the class.

### What this port does

Reproduces it, and says so. [`tepsim_core::balances::Balances`] carries the
trip and a `frozen` flag alongside the derivative, so a caller never has to
infer a freeze from a vector of zeros. B-0024a's typed
[`tepsim_core::measurements::ShutdownCause`] means it can also say *which*
limit fired, which the original cannot.

The fix is [`tepsim_core::balances::QuirkFixes::trip_ends_the_run`], **on by
default since 2026-08-28**. `PLAN.org` required "a full Tier 5 and Tier 6 delta
report and an explicit sign-off before it becomes the default"; both are below.

Reproducing the quirk is one call: `QuirkFixes::faithful()`, or
`Scenario::faithful()` for a whole run, or `tep run --freeze-on-trip`. Every
comparison against the Fortran or against published data uses it, and
`tier5::run_port` pins it so no differential can accidentally run the fix.

### The sign-off, and what decided it

Two facts, one from the delta measurement and one from the published files.

The fix is **pure truncation**. `d007_changes_nothing_before_the_trip` shows
every sample up to the trip is bit-identical either way, so the flag decides
how many numbers there are and never what they are. That removes the usual
worry about a Class C fix, which is that it quietly changes results.

And the frozen tail is not a small artefact. Four of the forty-four published
files carry one, 1,832 rows in total:

| file | rows | frozen tail | share | from row |
|---|---|---|---|---|
| `d06_te.dat` | 960 | 682 | 71.0% | 278 |
| `d06.dat` | 480 | 363 | 75.6% | 117 |
| `d18_te.dat` | 960 | 571 | 59.5% | 389 |
| `d18.dat` | 480 | 216 | 45.0% | 264 |

Across those rows twenty-one continuous channels repeat the same five-digit
values without any change at all, because the measurement noise at
`teprob.f:711-716` is inside the shutdown guard and stops with everything else.
Only the three dead-time analysers keep moving. So three quarters of `d06.dat`
is a stopped plant that reads as an unusually steady one, and a detector
trained on that file spends most of its evidence on it. Data that cannot be
told apart from good data, and is not, is worse than no data.

The counter-argument was restartability, and it does not survive contact: the
plant cannot be restarted after a freeze either. The freeze does not preserve
an option, it only withholds the fact that the run is over.

Note that B-0025's backlog entry originally said the freeze should *not* be the
default. That is the opposite of what `PLAN.org` says, and `PLAN.org` is the
design of record. It is also the only reading Tier 2 can live with: the
adversarial pool contains thirteen states that trip, and a port that did not
freeze would disagree with the oracle on all fifty components for each of them.

### Measured effect

- The freeze fires on exactly the states the Fortran freezes, over all three
  pools: 2,412 running and 13 frozen, with no disagreement.
- Turning the fix on changes the derivative on all 8 tripping adversarial
  boundaries and on none of the 17 that do not trip.
- Every sample before the trip is bit-identical with the fix on and off, so
  the delta is entirely truncation. `tier10_quirk_deltas.rs` tabulates which
  scenarios trip, at which step, and what fraction of the run is discarded.
- The Tier 5 battery **cannot** measure this delta, and that is a property of
  the delta rather than of the battery: every statistic compares two ensembles
  of the same shape, and the fix makes one of them shorter. A KS statistic
  between a 960-sample run and a 278-sample one measures the truncation.

### A second, open question this exposed

Comparing the assembled derivative revealed that **28 of the 50 components
exceed Tier 2's 1e-12 relative gate under the vendored libm**, worst `YP(2)` at
1.393e-4 — while the whole right-hand side is *bit-identical* to the Fortran
under `libm-system`.

That is cancellation, not error: a balance is inflow minus outflow, and near
steady state those nearly agree, so a one-ULP difference in each term is 1e-16
of the terms and 1e-4 of the result. The 22 components that do meet 1e-12 are
exactly the ones that do not cancel.

Choosing what Tier 2 should measure against for a cancelling quantity changes
what Tier 2 means, so it is recorded here and left open. **B-0026a.**

### Measured by

`crates/tepsim-oracle/tests/tier2_balances.rs`, and
`the_fix_is_off_by_default_and_changes_the_answer_when_on` in
`crates/tepsim-core/src/balances.rs`.

---

## D-008 — `CONTRL22` is defined, tuned, and never called

**Class A.** `temain_mod.f:1295-1332`, with its constants at `246-317`.

### What the original does

Twenty controller subroutines are defined. The main loop calls nineteen of
them. `CONTRL22` is not among the calls.

It is not a stub. The main program initialises its full tuning alongside every
other loop's:

```fortran
      SETPT(12)=2633.7
      GAIN22=-1.0	  * 5.
      TAUI22=1000./3600.
      ERROLD22=0.0
```

and the subroutine is a complete PI controller reading `XMEAS(13)`, the
separator pressure, and writing `XMV(6)`, the purge valve.

### What it looks like it was

`XMV(6)` is the valve `CONTRL6` already owns, and `CONTRL6` carries a latching
pressure override (`temain_mod.f:710-731`) that does the same job by a
different mechanism. The override's release threshold is **2633.7**, which is
exactly `CONTRL22`'s setpoint.

`CONTRL22` is also the only one of the twenty whose error is not normalised by
a span: every other loop computes `(SETPT - XMEAS) * 100 / span`, and this one
computes the raw difference.

Two loops on one valve, one of them differently shaped from all the others, and
one of them disconnected. The straightforward reading is that `CONTRL22` was
the separator-pressure controller, that the override in `CONTRL6` replaced it,
and that it was left in place with its tuning intact rather than deleted.

### Why it matters

Because "nineteen control loops" and "twenty control subroutines" are both true
and a port has to pick. Counting the subroutines gives a plant with two
controllers fighting over the purge valve. Counting the calls gives the
published behaviour.

It also shifts the reading of the control structure. Separator pressure is
*not* under continuous control in this scheme; it is under an on-off override
with a 650 kPa deadband. Anyone reasoning about why the plant behaves as it
does around the pressure limits needs that.

### What this port does

Ports it, and does not schedule it. It is in
[`tepsim_control`]'s tuning table and in the Tier 1 differential, because
covering it costs nothing and it is the only check that will ever exercise it,
but the scheduler does not call it.

### Measured effect

None: it is not called, so it computes nothing.

The differential covers it anyway, at 0 ULP over 500 tunings, which is the only
evidence anywhere that the subroutine is correct.

### Measured by

`every_controller_matches_the_fortran_over_a_sweep` in
`crates/tepsim-oracle/tests/tier1_control.rs` includes `CONTRL22`.
`crates/tepsim-oracle/tests/driver_binding.rs` shows it is callable and that
the driver's setpoint array has a slot for it.

## D-009 — the driver starts from rounded valve positions, not `TEINIT`'s

**Class A.** `temain_mod.f:322-332`.

### What the original does

`TEINIT` leaves the twelve valve commands at the values `YY(39..50)` carries,
which are written to eight significant figures:

```fortran
      DATA YY /
     .  ...
     .  6.3053638D+01, 5.3980356D+01, 2.4644630D+01, ...
```

The driver then overwrites eleven of the twelve, by hand, at five:

```fortran
 	XMV(1) = 63.053 + 0.
	XMV(2) = 53.980 + 0.
	XMV(3) = 24.644 + 0.
	...
	XMV(11)= 18.114 + 0.
```

`XMV(12)`, the agitator, is not in the list and keeps `TEINIT`'s exact 50.

Every literal is fixed-form Fortran with no exponent letter, so each is a
`REAL(4)` widened to double. `63.053` reaches the plant as
`63.053001403808594`, not as `63.053`.

### Why it matters

A closed-loop run and an open-loop run do not start from the same plant. The
difference is in the fourth decimal place of ten valve commands, which is far
too small to see in any plot and far too large to ignore in a differential: a
port that starts a closed-loop run from `TEINIT`'s values disagrees with the
Fortran from step 1, and the disagreement then grows through the controllers.

Ten of the eleven, not all eleven: `TEINIT` leaves `YY(43) = 22.21000000`,
which is *already* five significant figures, so `22.210` rounds to the same
`f32`. That coincidence is the clearest available evidence that these numbers
are `TEINIT`'s rounded rather than an independently chosen operating point.

The `+ 0.` on every line appears to be a placeholder for a perturbation. It
changes nothing, and it is transcribed rather than simplified away, because
`single(63.053 + 0.)` and `single(63.053) + single(0.)` are the same value only
because the addend is zero.

### What this port does

Reproduces it, as [`tepsim_control::DRIVER_INITIAL_VALVES`], with each literal
passed through `single()`. `Driver::new` starts there.

### Measured effect

Ten of the twelve valve commands differ from `TEINIT`'s, all by less than
one part in a thousand of range:

| | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| gap | 3.70e-4 | 2.94e-4 | 4.41e-4 | 7.63e-5 | 0 | 2.52e-4 | 3.43e-4 | 1.56e-4 | 2.63e-4 | 1.87e-4 | 5.09e-4 | 0 |

The largest is `XMV(11)` at 5.09e-4; the smallest non-zero is `XMV(4)` at
7.63e-5. `XMV(5)` and `XMV(12)` are bit-identical, for the two different
reasons above.

### Measured by

`the_driver_starts_from_rounded_valve_positions` in
`crates/tepsim-control/src/lib.rs`, which asserts the count of differing valves
and the two exceptions individually, and
`the_rounding_gaps_are_recorded` beside it, which prints the table above.

## D-010 — the controllers read the previous step's measurements

**Class A, and load-bearing.** `temain_mod.f:366-411`.

### What the original does

The main loop, in order:

```fortran
        DO 1000 I = 1, NPTS
	  TEST=MOD(I,3)
	  IF (TEST.EQ.0) THEN
		CALL CONTRL1
	  	...
	  ENDIF
          ...
	  CALL INTGTR(NN,TIME,DELTAT,YY,YP)
 	  CALL CONSHAND
 1000 CONTINUE
```

`XMEAS` is written by `TEFUNC`, which `INTGTR` calls. So on iteration `I` the
controllers read the measurements iteration `I - 1` produced. Every loop in the
scheme carries one plant step of dead time that is nowhere in any controller.

`CONSHAND`, the valve clamp, likewise runs *after* the integration rather than
after the controllers.

### Why it matters

This is not a subtlety that costs a few ULP. Feeding the controllers the
measurements of the step they are about to cause makes every loop one sample
tighter than the original. B-0039 measured it: `XMV(7)` lands on 35.62 instead
of 34.15 on the very first controller fire, a 1.5% of range error before the
plant has run three seconds, and `XMEAS(14)` is 23% out four hours later. That
is a different plant, not a rounding of the same one.

The clamp's placement is subtler. `TEFUNC` clamps its own copy of the valve
positions at `teprob.f:803-804`, so with no sticking fault active it makes no
observable difference where `CONSHAND` runs. It stops being unobservable under
`IDV(14)`, `IDV(15)` or `IDV(19)`: `teprob.f:801` only moves `VCV` toward `XMV`
when the two differ by more than the stick threshold, and an unclamped `XMV` of
105 crosses that threshold at a different moment than a clamped 100 does.

### What this port does

Reproduces both, in [`tepsim_control::Driver`] rather than in `Scheme`.
`Driver::control` takes the previous step's measurements and does not clamp;
`Driver::settle` is `CONSHAND` and is called with the plant already advanced.
The split exists so that the ordering is a thing the type makes explicit rather
than a convention a caller has to remember.

### Measured effect

With the two orderings otherwise identical, the first `CONTRL7` fire gives
34.147823842 (previous-step, matching the Fortran bit for bit) against
35.623679189 (current-step).

### Measured by

`the_controllers_read_the_previous_steps_measurements` in
`crates/tepsim-oracle/tests/tier4_closed_loop.rs`, which requires the two
orderings to give *different* answers before checking which one is right.

## D-011 — the driver switches `IDV(12)` on eight hours in, whatever you asked for

**Class C. Not reproduced by default since the sign-off of 2026-08-28 (B-0066);
reproduced by `Scenario::faithful`.** `temain_mod.f:366-368`, with `SSPTS` at
line 226.

### What the original does

The first statement in the simulation loop body, before any controller:

```fortran
        DO 1000 I = 1, NPTS
         IF (I.GE.SSPTS) THEN
                 IDV(12)=1
          ENDIF
```

with

```fortran
      SSPTS = 3600 * 8
```

At a one-second step that is eight simulated hours. `IDV(12)` is the condenser
cooling water inlet temperature random variation. The assignment is
unconditional: it does not consult the scenario, and there is no way to run
this driver past eight hours without it.

The header comment at line 98 explains `SSPTS` as "the number of data points to
simulate in steady state operation before implementing the disturbance", so a
disturbance was clearly meant. But which disturbance is not a parameter; it is
`IDV(12)`, written into the loop.

### Why it matters

Every published closed-loop dataset generated with this driver and longer than
eight hours carries `IDV(12)`, including the runs nominally labelled
fault-free. Anyone comparing against published data needs the quirk; anyone
running a controlled experiment on some *other* disturbance needs it gone,
because after hour eight they are running two disturbances and reporting one.

It also interacts with the fault the caller asked for. If the scenario is
`IDV(4)` (a step in reactor cooling water inlet temperature), then from hour
eight the run is `IDV(4)` **and** `IDV(12)`, and both act on cooling water.

### The effect does not begin at hour eight

`IDV(12)` reaches the plant through `IDVWLK(6)` (`teprob.f:351`), and the walk
that gates is only redrawn when `TIME` passes that channel's `TNEXT`
(`teprob.f:359-360`). So switching the flag on at step 28,800 changes nothing
until channel 6's next segment boundary, which for the nominal seed is step
**29,390**, about ten minutes later. Both the port and the Fortran part from
their fixed counterparts at exactly that step.

That matters for anyone trying to locate the quirk in a dataset by eye: the
visible onset is not at the eight-hour mark, and it moves with the seed.

### What this port does

Can do either, and does not by default. [`tepsim_control::Driver`] forces
`IDV(12)` on at [`tepsim_control::STEADY_STATE_STEPS`], and
[`tepsim_control::DriverQuirks::only_the_requested_disturbances`] turns that
off; [`tepsim::Scenario::driver_forces_idv12`] is the facade's control and is
`false` by default. `Scenario::faithful()` and `tep run --force-idv12` reproduce
the driver as shipped, and `tier5::run_port` pins it on so every differential
against the Fortran carries it.
[`tepsim_control::Driver::scenario_is_overridden`] reports when the driver has
gone beyond what was asked for, so a caller never has to infer it.

### Measured effect

Ten simulated hours, nominal scenario, one-second Euler, faithful against
fixed, everything else identical:

| | |
|---|---|
| first difference | step 29,390, `XMEAS(22)` (condenser cooling water outlet temperature) |
| identical before that | every bit of all 41 measurements, all 29,389 steps |
| worst over the run | 1.092e-1 relative, at `XMEAS(37)` |
| worst at hour ten | 3.262e-2 relative, at `XMEAS(38)` |

Ten percent on a product-composition measurement is not a rounding difference.
It is a different experiment.

### The sign-off, and what decided it

The question was whether the shipped default should be `temain_mod.f`'s
behaviour or the caller's scenario. Decided on 2026-08-28: **the caller's
scenario.**

**One of the two original arguments had already been shown to be false.** The
case for keeping the quirk was that reproducing published data requires it.
Tier 7 (B-0051) established that the published `d00` through `d21` files carry
`IDV(12)` only in `d12` and `d12_te`: they were generated with
`temain_mod.f:367` *replaced*, not kept.

Two independent predictions support that. `d12_te`'s spread jumps at row 160,
where `IDV(12)` would arrive, by at least 5.3 times what `d00_te`'s does on
every channel the disturbance drives. And keeping the line inflates the port's
spread against `d00` by up to 9.9 times, against 1.55 times for replacing it.
Removing it roughly triples the agreement across the whole comparison table,
for example `d17_te` from a Kolmogorov-Smirnov median of 0.658 to 0.044.

So the decision is now between two different goals rather than between fidelity
and convenience:

- **Reproduce the source.** `temain_mod.f` as shipped does force `IDV(12)`, and
  this project is a port of that source. Keeping the default means the port
  does what the code it was ported from does.
- **Reproduce the data.** Most people who use the Tennessee Eastman problem use
  the published datasets, not the driver. Matching them means turning it off.

The second was chosen, on three grounds. The evidence is one-sided: the only
argument for forcing was that the shipped line exists, and the published bytes
say the line was replaced when the data was made. The prose at
`temain_mod.f:101-102` agrees, telling the reader to "go to line 367" and put
their own disturbance there, which is an instruction to replace. And the cost
of being wrong is asymmetric: a caller who asks for `IDV(4)` and silently gets
`IDV(4)` and `IDV(12)` has no way to notice, whereas a caller who wants the
shipped driver's behaviour asks for it by name and gets exactly that.

Both remain reachable, which was the one firm constraint: `Scenario::faithful`,
`tep run --force-idv12`, `driver_forces_idv12=True` from Python, and the
`idv12` field of the scenario text, which still parses old links to the
scenario they always named.

The numbers above are Tier 4. The Tier 5 measurement `PLAN.org` asks for on a
Class C delta is `d011_the_forced_disturbance_moves_the_plant_by_a_tenth` in
`crates/tepsim-oracle/tests/tier10_quirk_deltas.rs`, which runs the paired
battery with the flag on and off and requires the shift to exceed a tenth of a
margin on some channel: forcing `IDV(12)` is not cosmetic, and the test fails if
it ever looks that way.

### Measured by

`the_driver_forces_idv12_at_the_eight_hour_mark` and
`the_forced_disturbance_changes_the_plant_measurably` in
`crates/tepsim-oracle/tests/tier4_closed_loop.rs`. The second cross-checks the
onset step against two Fortran runs that differ in the same way, so 29,390 is
ground truth rather than an artifact of the port.

`the_published_files_were_not_generated_with_the_forced_idv12` in
`crates/tepsim-oracle/tests/tier7_published.rs` is what established that the
published datasets do not contain it.
