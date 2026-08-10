# TEP-Rust: session operating manual

Pure-Rust port of the Tennessee Eastman Process simulator, with Python wheels and
a browser (wasm) app. The full design is in `PLAN.org`. **Do not read all of
`PLAN.org` every session.** Read the section named by the backlog item you are
working, plus this file.

This project is built by a series of independent sessions. Each session performs
**exactly one iteration**: one backlog item, taken all the way to done, with
tests and validation, committed green. This file is the protocol for that.

---

## The loop

### 1. Orient (do this first, always)

```bash
git -C . status --short && git log --oneline -8
sed -n '1,60p' BACKLOG.org          # current-state block + top of queue
tail -n 60 LOG.org                  # last 1-2 iterations, esp. their "Next" section
cargo xtask ci                      # confirm the baseline is green
```

Three branch points:

- **Baseline is red.** The iteration is "make it green", nothing else. Do not
  start new work on a red tree.
- **Working tree is dirty and an item is marked `NEXT`.** A previous session
  died mid-iteration. Read its partial work, then either finish that item or
  `git restore` and re-start it. Say which in the log.
- **Clean and green.** Proceed.

### 2. Fidelity preflight

Green says nothing about whether fidelity to the Fortran degraded. Each session
sees only its own item, but drift is a whole-project property, so it has to be
re-established at the top of every cycle. Skip this only while no model code
exists yet (through B-0004).

```bash
gfortran --version | head -1                    # record this in the log entry
cargo xtask fidelity                            # golden 100-step oracle diff, ~2 s
cargo xtask validate --tiers 1,2,3 --compare-to-log
cargo xtask provenance                          # unclaimed teprob.f line ranges
```

Then read the open questions in `book/src/deltas.md`, so you do not re-litigate a
quirk that already has a decision.

Four rules make this worth running:

- **A degradation inside tolerance is still this iteration's work.** If Tier 2
  max relative error moves 3e-14 to 8e-13, both pass the 1e-12 gate, but
  something broke. The only place that is visible is the logged history. Stop and
  find it; do not proceed to a new item.
- **`gfortran` changed since the last log entry?** Re-baselining Tiers 1 and 2 is
  its own logged event, not a regression hunt. The reference numbers depend on the
  compiler, so a Homebrew upgrade looks exactly like a bug and will otherwise eat
  a whole session.
- **No `gfortran`, no model work.** If the oracle cannot build on this machine,
  restrict the iteration to infrastructure, docs, Python, or wasm items.
- **Unclaimed `teprob.f` ranges are a finding.** `xtask provenance` is how a
  silently dropped term in a balance equation gets caught. No differential test
  will find a term that is never evaluated.

### 3. Pick exactly one item

Take the topmost item in `BACKLOG.org` that is `TODO` and whose `:DEPS:` are all
`DONE`. Flip it to `NEXT` and commit that flip immediately, so a crashed session
leaves a trace.

**Sizing rule.** If you cannot finish it *including tests and validation* in this
session, do not start it. Split it in place into child items, take the first
child, and record the split in the log. A half-finished item is the one failure
mode that makes the next session expensive.

### 4. Work it

Test first where it is natural: write the failing test, then make it pass. Where
it is not natural, still prove the test has teeth by stubbing out the change and
watching it fail before you commit.

Consult the `PLAN.org` section named in the item's `:PLAN:` property. If the plan
and the code disagree, the plan wins, or you record a decision entry and update
the plan. Never silently diverge.

**Read the Fortran, not the port.** For any item touching the model, read the
source line range in `reference/fortran/teprob.f` *in this session*, before
reading any existing Rust. Do not port from the upstream Python implementation,
and do not port by analogy to a Rust function that looks similar. This exists to
stop a misreading in one session from being copied forward through the next
thirty.

**Claim what you port.** Every ported function gets a claim on its own comment
line, immediately above it:

```rust
// @port teprob.f:505-522
fn reaction_rates(/* ... */) { /* ... */ }
```

`cargo xtask provenance` collects these and reports what nothing accounts for.
The marker is anchored at the start of the comment on purpose, so prose about
the convention and test fixtures cannot inflate coverage.

**Constants are asserted, never retyped.** Transcribe once, then prove equality
against the oracle's `COMMON` blocks. Digits read off a listing are a silent
failure mode with no test that catches them.

### 5. Close out

Every one of these, in order. No exceptions, no partial credit.

- [ ] `cargo xtask ci` green
- [ ] The item's `:TIER:` validation exists, runs, and passes
- [ ] **Numbers recorded**, not verdicts (see below), including the preflight
      numbers and the `gfortran` version they were produced with
- [ ] `cargo xtask provenance` reports no *new* unclaimed `teprob.f` ranges
- [ ] Every function ported this iteration carries its source line range
- [ ] Public API has rustdoc; unit-operation modules have equations, the
      variable table, and the `teprob.f` line range they came from
- [ ] Backlog item flipped `NEXT` → `DONE` with a `:RESULT:` property
- [ ] `LOG.org` entry appended, including the `** Next` section
- [ ] **Two commits**, in this order, working tree clean after:
      1. `[B-00NN] <summary>` — the work itself
      2. `[B-00NN] Record iteration` — `BACKLOG.org` and `LOG.org`, with
         `:COMMIT:` and `:RESULT:` naming the hash of commit 1

Two commits, not one, because a commit cannot contain its own hash. Write
`pending` in the hash fields, make the work commit, fill in the real hash, then
make the bookkeeping commit. Never `--amend` after recording a hash: it
invalidates the value you just wrote.

---

## Record numbers, not verdicts

"Tier 2 passing" is worthless to the next session. Write down what was actually
measured, so drift is visible across iterations:

> Tier 2: 42,000 sampled states, max relative error 3.1e-13, worst component
> `YP(27)` (stripper energy) at state `adversarial/tcc-near-170`, ULP histogram
> p99 = 2, p100 = 11.

This is the single highest-value habit in this project. Every validation number
goes in the log entry, and regressions get caught by comparing against the
previous entry rather than by a threshold that someone quietly relaxed.

---

## Rules that override convenience

- **Never loosen a tolerance, skip a test, add `#[ignore]`, or widen an
  equivalence margin to make something pass.** If a tolerance looks wrong, stop,
  log a `BLOCKED` item with the evidence, and move to the next item.
- **Never edit `reference/`.** The vendored Fortran and the `d00`–`d21` files are
  ground truth. Checksums are asserted in CI.
- **Class C quirk fixes need sign-off.** (Shutdown derivative-zeroing, binary-only
  IDV flags, fixed-step Euler, hard-coded `IDV(12)` in the driver.) Implement
  behind a flag, measure the delta with the full Tier 5 battery, log the numbers,
  mark the item `BLOCKED` on a decision, and move on. Do not make it the default.
- **`tepsim-core` is `#![forbid(unsafe_code)]`** and `no_std + alloc`. `unsafe` is
  allowed only in the PyO3 and wasm-bindgen glue.
- **Determinism is a hard invariant.** No `f32`, no SIMD or rayon inside the core,
  no reordered reductions, no `Date`/time/randomness outside `TepRng`. The core
  uses the vendored `libm` crate for `exp`/`pow`/`ln`.
- **`tepsim-oracle` never becomes a dependency** of `tepsim`, `tepsim-py`, or
  `tepsim-wasm`. It is dev-only, behind the `oracle` feature.
- **The oracle's compiler and flags are pinned.** `build.rs` fixes the gfortran
  flags explicitly and a test asserts them. Never add `-ffast-math` or anything
  that permits reassociation. Changing the flags invalidates every recorded Tier
  1 and Tier 2 number, so it is a logged re-baseline, never a casual edit.
- Do not pivot the architecture mid-iteration. If the design in `PLAN.org` is
  genuinely wrong, stop, write a decision entry, and ask.

---

## File map

| File | Role | Who writes it |
|---|---|---|
| `PLAN.org` | Design of record. Architecture, validation tiers, roadmap. | Changed only via a logged decision |
| `BACKLOG.org` | Ordered work queue + current-state block. One heading per iteration-sized item. | Every iteration |
| `LOG.org` | Append-only history. One heading per iteration. | Every iteration |
| `CLAUDE.md` | This protocol. | Rarely |
| `tep-rust.bib` | org-ref bibliography for `PLAN.org`. | Rarely |
| `reference/` | Vendored original Fortran and published datasets. | Never |
| `book/` | mdBook: theory, unit operations, generated validation report. | As features land |

---

## Commands

```bash
cargo xtask ci            # THE gate: fmt --check, clippy -D warnings, test, doc, deny
cargo xtask ci --fast     # skip the oracle differential job (local iteration only)
cargo xtask fidelity      # preflight: 100 steps vs a committed golden oracle trace
cargo xtask provenance    # teprob.f line ranges not claimed by any Rust function
cargo xtask validate      # full ladder; regenerates book/src/validation/
cargo xtask validate --tiers 1,2,3 --compare-to-log
cargo xtask bench         # criterion, with regression comparison
cargo test -p tepsim-oracle --features oracle    # needs gfortran
```

`xtask fidelity` deliberately does **not** need gfortran: it diffs against a
golden trace committed to the repo, so it runs in seconds on every machine and
catches gross breakage. The live oracle diff in `xtask validate` is the periodic
deep check.

Availability: `xtask ci` and `xtask provenance` arrive in **B-0002**, `xtask
fidelity` in **B-0004**, `xtask validate` incrementally from **B-0008**. Before
B-0002 the gate is `cargo fmt --check && cargo clippy --all-targets -- -D
warnings && cargo test --workspace`, and the fidelity preflight is skipped
because there is no model code to protect.

---

## Environment

Verified present on this machine: `git` 2.50, `cargo`/`rustc` 1.89 (Homebrew),
`gfortran` 15.2 (Homebrew GCC), `python3` 3.13, `maturin` 1.12.4, `cargo-deny`.

Missing, and needed later:

- **`rustup`** is not installed; Rust comes from Homebrew. This blocks
  `rust-toolchain.toml` pinning and the `wasm32-unknown-unknown` target.
  Handled by **B-0001**. Installing it changes the user's system, so **ask
  before doing it.**
- **`trunk`** is not installed. Needed only for Phase 8 (the browser app).

The oracle needs `gfortran`, which is present here. In CI it runs on Linux and
macOS runners only, never Windows.

---

## Validation tiers (quick reference; details in `PLAN.org`)

| Tier | What it proves | Gate |
|---|---|---|
| 1 | `TESUB1`–`TESUB8` match the oracle | rel err < 1e-13; `TESUB7` bit-exact |
| 2 | Single-step derivatives match | rel err < 1e-12 over all three sampling pools |
| 3 | RNG call **order** matches | trace diff is empty |
| 4 | Trajectories (diagnostic, not a gate) | error < `XNS(i)` for the first hours; divergence attributed to libm |
| 5 | Statistical equivalence | TOST on means **and variances**, KS, ACF, Welch spectra, correlation matrix |
| 6 | Downstream detectors cannot tell the sources apart | cross-source ≈ within-source |
| 7 | Reproduces published `d00`–`d21` | per-file Tier 5 report |
| 8 | Differential fuzzing | nightly, no counterexamples |
| 9 | Cross-platform determinism | identical BLAKE3 digests incl. wasm in a browser |
| 10 | Every quirk fix has a measured delta | no statistic outside its equivalence margin |

Tier 4 is **diagnostic**. Do not treat long-horizon trajectory divergence as a
bug; it is expected from `exp`/`pow` ULP differences. Tier 5 and 6 are the gates.

---

## Log entry format

Append to the end of `LOG.org`:

```org
* 2026-08-10 B-0009 Tier 1 harness: sweep generators and ULP reporting
:PROPERTIES:
:ITEM:     B-0009
:PHASE:    1
:TIER:     1
:COMMIT:   a1b2c3d
:GFORTRAN: 15.2.0
:STATUS:   DONE
:END:
** Preflight
Tier 1/2/3 headline numbers as found at the start, and whether they moved
against the previous entry. Unclaimed teprob.f ranges, if any.
** What changed
Two or three sentences. What exists now that did not before.
** Validation numbers
The actual measurements. Counts, max errors, percentiles, p-values.
** Surprises and decisions
Anything the next session would be annoyed not to know. Empty is fine.
** Next
The single most useful thing to do next, and why. This is what the next
session reads first.
```

Use `STATUS: BLOCKED` and a `** Blocked on** section when an item cannot finish.
State the question precisely enough that a one-line answer unblocks it.

---

## When blocked

Do not guess, and do not work around it. Append a `BLOCKED` log entry with the
question and the evidence, flip the backlog item to `BLOCKED` with a `:QUESTION:`
property, take the next unblocked item, and surface the question in your closing
message to the user.
