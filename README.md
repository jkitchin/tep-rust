# TEP-Rust

A pure-Rust implementation of the Tennessee Eastman Process (TEP) simulator,
delivered three ways: a Rust crate, a Python package shipped as binary wheels,
and a self-contained WebAssembly application that runs entirely in the browser
with no server.

TEP is the standard benchmark for plant-wide process control and for
fault-detection research: a reactor, condenser, vapor-liquid separator,
recycle compressor, and stripper, with eight chemical species, four reactions,
50 states, 41 measurements, 12 manipulated variables, and 20 canonical
disturbances.

## Status

**Phase 0, under construction. There is no working simulator yet.**

This repository is being built one complete, validated increment at a time. See
[`BACKLOG.org`](BACKLOG.org) for what is next and [`LOG.org`](LOG.org) for what
has been done and what was measured.

## Why another TEP

Existing ports are faithful in spirit but not verifiable in practice. The
canonical Fortran is fully deterministic given a seed, which means a port can be
held to a far sharper standard than "the statistics look similar". This one aims
for:

- **Provable equivalence.** A ten-tier validation ladder running against the
  original Fortran as a live oracle, from bit-exact utility routines through
  derivative-level agreement to statistical equivalence testing and, finally, to
  showing that fault detectors cannot tell which simulator produced their data.
- **Reproducibility as a hard invariant.** Identical results across x86-64,
  aarch64, and wasm32, enforced in CI by golden digests.
- **Faults as data, not code.** The 20 canonical disturbances are expressed as
  (injection point, profile) pairs, so users can define arbitrary custom faults,
  set continuous magnitudes, compose them, and schedule them.
- **Machine-readable ground truth,** so detection-delay metrics stop being
  guesswork.

## Planned interfaces

Simulation mode sets conditions and runs forward. Control mode lets users write
controllers in Rust, in Python, or, in the browser app, in a small scripting
language with no toolchain required. Disturbances can be scheduled, or defined
from scratch against any injection point in the plant.

## Documents

| File | What it is |
|---|---|
| [`PLAN.org`](PLAN.org) | Design of record: architecture, validation strategy, roadmap |
| [`BACKLOG.org`](BACKLOG.org) | Ordered work queue and current state |
| [`LOG.org`](LOG.org) | Iteration history with measured validation numbers |
| [`CLAUDE.md`](CLAUDE.md) | Development protocol for this repository |
| [`NOTICE.md`](NOTICE.md) | Attribution, upstream license, and citation requirements |

## License

MIT ([`LICENSE`](LICENSE)). Portions are derived from Fortran licensed under the
University of Illinois/NCSA Open Source License ([`LICENSE-NCSA`](LICENSE-NCSA)),
whose attribution conditions apply to source and binary redistribution alike and
cannot be dropped. The combined work is `MIT AND NCSA` in SPDX terms. See
[`NOTICE.md`](NOTICE.md) before redistributing.

The original process and control problem are due to J. J. Downs and E. F. Vogel
(Tennessee Eastman Company); the widely used modified code is due to the Large
Scale Systems Research Laboratory at the University of Illinois under Prof.
Richard D. Braatz.
