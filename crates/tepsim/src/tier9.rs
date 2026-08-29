//! Tier 9: cross-platform determinism.
//!
//! The claim this project rests on is that a [`Scenario`] is a *description* of
//! a dataset rather than a recipe that happens to produce one on the machine it
//! was run on. `no_std`, the vendored `libm`, no `f32`, no SIMD, no reordered
//! reduction and no clock all exist to make that true. This module is where the
//! claim is stated as a number and checked.
//!
//! # How it is checked
//!
//! [`CASES`] is a table of fixed scenarios, each with a **committed constant**:
//! the digest that scenario produced when it was first measured. Every platform
//! runs the same table and compares against the same constants. That is the
//! whole design, and the alternative (compute a digest twice in one process and
//! compare) would prove only that the machine agrees with itself.
//!
//! Because the constants are committed, a platform completes its half of the
//! claim on its own. `cargo test -p tepsim` fails on a machine whose native
//! arithmetic disagrees. `cargo xtask tier9` additionally builds the module for
//! `wasm32-unknown-unknown` and has a WebAssembly runtime evaluate the same
//! table, so one command covers the host architecture and wasm32 together.
//! Running it on a second architecture is what extends the claim to that
//! architecture; nothing here needs changing to do it.
//!
//! # Why FNV-1a and not BLAKE3
//!
//! `PLAN.org` says BLAKE3. This is FNV-1a, deliberately, for four reasons.
//!
//! The digest has to be computed *inside* the artifact on every target,
//! including a browser, because a digest computed outside would be a digest of
//! whatever crossed the boundary rather than of what the simulator produced.
//! That puts the hash in a `no_std` crate that compiles to wasm32 under an
//! explicit size budget (`PLAN.org`, "Budgets"), and BLAKE3 is a large thing to
//! spend that budget on.
//!
//! Moving BLAKE3 into `xtask` instead would not help. What `xtask` compares is
//! two 64-bit numbers produced on two different machines. Hashing a `u64` with
//! a stronger hash does not make comparing it stronger.
//!
//! FNV-1a over [`f64::to_bits`] is integer arithmetic end to end, so the digest
//! cannot itself become a source of the cross-platform disagreement it exists
//! to detect. A hash that did any floating-point arithmetic could.
//!
//! And the threat model is drift, not an adversary. Nobody is constructing a
//! second trajectory that collides with the first; the failure being watched
//! for is a compiler, a libm or an instruction selection quietly changing one
//! bit. A 64-bit hash over every bit of every emitted value catches that with a
//! probability that rounds to one.
//!
//! It is not cryptographic and not collision-resistant, and saying so is better
//! than implying a guarantee that is not being provided.
//!
//! # What is hashed
//!
//! Each recorded sample contributes `[hours, XMEAS(1..41), XMV(1..12)]`, in
//! that order, as IEEE 754 bit patterns, little-endian. That is exactly the row
//! `tepsim-wasm` hands a browser, so the digest covers the numbers that
//! actually leave the library rather than an internal form that a packing bug
//! could diverge from.
//!
//! Bit patterns rather than values, and with no normalisation: `-0.0` and `0.0`
//! are different here on purpose. [`tepsim_scenario::Digest`] normalises them
//! because it identifies an *experiment*, and two descriptions of the same
//! experiment should hash alike. This identifies an *output*, where a signed
//! zero appearing on one architecture and not another is precisely the finding
//! Tier 9 exists to surface.

use crate::Integrator;
use crate::recorder::Recorder;
use crate::run::Sample;
use crate::scenario::Scenario;
use crate::sim::Simulation;

/// The 64-bit FNV-1a hash of a stream of `f64` bit patterns.
///
/// Used for two jobs that must agree: the canonical Tier 9 digest here, and
/// `tepsim-wasm`'s running checksum over the chunks it hands a browser. One
/// implementation, so the two cannot drift apart.
///
/// See the module documentation for what this is and is not.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Fnv1a64(u64);

impl Fnv1a64 {
    /// The FNV-1a 64-bit offset basis.
    const OFFSET_BASIS: u64 = 0xcbf2_9ce4_8422_2325;
    /// The FNV-1a 64-bit prime.
    const PRIME: u64 = 0x0000_0100_0000_01b3;

    /// An empty digest.
    #[must_use]
    pub const fn new() -> Self {
        Self(Self::OFFSET_BASIS)
    }

    /// Absorb one byte.
    #[inline]
    pub const fn write_u8(&mut self, byte: u8) {
        self.0 ^= byte as u64;
        self.0 = self.0.wrapping_mul(Self::PRIME);
    }

    /// Absorb the eight little-endian bytes of a value's bit pattern.
    ///
    /// Bits, not the value, so byte order is fixed by IEEE 754 rather than by
    /// the host, and so `-0.0` stays distinguishable from `0.0`.
    #[inline]
    pub const fn write_f64(&mut self, value: f64) {
        let bytes = value.to_bits().to_le_bytes();
        let mut i = 0;
        while i < bytes.len() {
            self.write_u8(bytes[i]);
            i += 1;
        }
    }

    /// Absorb every value of a slice, in order.
    #[inline]
    pub fn write_slice(&mut self, values: &[f64]) {
        for value in values {
            self.write_f64(*value);
        }
    }

    /// Absorb a boolean as one byte.
    #[inline]
    pub const fn write_bool(&mut self, value: bool) {
        self.write_u8(value as u8);
    }

    /// Absorb the eight little-endian bytes of an integer.
    #[inline]
    pub const fn write_u64(&mut self, value: u64) {
        let bytes = value.to_le_bytes();
        let mut i = 0;
        while i < bytes.len() {
            self.write_u8(bytes[i]);
            i += 1;
        }
    }

    /// Absorb a string's bytes, then its length.
    ///
    /// The length is in there so that `"ab"` then `"c"` cannot collide with
    /// `"a"` then `"bc"`, which matters for [`suite_digest`], where case names
    /// are hashed one after another.
    #[inline]
    pub const fn write_str(&mut self, value: &str) {
        let bytes = value.as_bytes();
        let mut i = 0;
        while i < bytes.len() {
            self.write_u8(bytes[i]);
            i += 1;
        }
        self.write_u64(bytes.len() as u64);
    }

    /// The digest so far.
    #[must_use]
    pub const fn finish(self) -> u64 {
        self.0
    }
}

impl Default for Fnv1a64 {
    fn default() -> Self {
        Self::new()
    }
}

/// Absorbs `[hours, XMEAS(1..41), XMV(1..12)]` from every sample, in order.
///
/// This is what makes [`digest`] a streaming computation: a nine-hour run is
/// 32,400 steps and its samples never have to be in memory at once. It is also
/// the layout `tepsim_wasm::Runner` hashes, which is what lets one constant
/// serve both.
impl Recorder for Fnv1a64 {
    fn record(&mut self, sample: &Sample) {
        self.write_f64(sample.hours);
        self.write_slice(&sample.row());
    }
}

/// One fixed scenario and the digest it is committed to produce.
///
/// The digest is a measurement that was written down, not a target that was
/// chosen. If one moves, the platform, the compiler or the model changed the
/// numbers, and the response is to find out which rather than to update the
/// constant. Moving one is a logged re-baseline, the same as a `gfortran` or
/// toolchain change.
#[derive(Clone, Copy, Debug)]
pub struct Case {
    /// A stable identifier, used in reports and hashed into [`suite_digest`].
    pub name: &'static str,
    /// What numerical path this case is in the table to exercise.
    pub covers: &'static str,
    /// Builds the scenario. A function rather than a value because
    /// [`Scenario`]'s builders are not all `const`.
    pub build: fn() -> Scenario,
    /// The committed digest. See the type documentation.
    pub digest: u64,
}

impl Case {
    /// The scenario this case runs.
    #[must_use]
    pub fn scenario(&self) -> Scenario {
        (self.build)()
    }

    /// Run the scenario and digest its output, here and now.
    #[must_use]
    pub fn compute(&self) -> u64 {
        digest(self.scenario())
    }

    /// Whether this platform reproduces the committed digest.
    #[must_use]
    pub fn agrees(&self) -> bool {
        self.compute() == self.digest
    }
}

fn baseline_1h() -> Scenario {
    Scenario::baseline().with_hours(1.0)
}

fn fault_1_2h() -> Scenario {
    Scenario::fault(1).with_hours(2.0)
}

fn open_loop_4h() -> Scenario {
    Scenario::baseline().with_hours(4.0).open_loop()
}

fn rk4_1h() -> Scenario {
    Scenario::baseline()
        .with_hours(1.0)
        .with_integrator(Integrator::Rk4)
}

fn dormand_prince_1h() -> Scenario {
    Scenario::baseline()
        .with_hours(1.0)
        .with_integrator(Integrator::DormandPrince)
}

fn baseline_9h() -> Scenario {
    Scenario::baseline().with_hours(9.0)
}

/// The Tier 9 table.
///
/// Six cases rather than one, because a single short fault-free run leaves most
/// of the arithmetic untested and a determinism check that only covers the easy
/// path is worth very little. Between them these reach the disturbance table,
/// the uncontrolled plant, the two multi-stage integrators, and a horizon long
/// enough to cross the driver's forced `IDV(12)` at eight hours.
///
/// Every duration is a whole number of hours, and the cadence divides an hour
/// exactly, so the sample count is unambiguous and a chunked consumer and a
/// batch one record the same rows.
pub const CASES: &[Case] = &[
    Case {
        name: "baseline-1h",
        covers: "closed loop, fault free, Euler: the controllers, the \
                 measurement noise generator and the analysers' dead time",
        build: baseline_1h,
        digest: 0xc8a2_6889_992f_1719,
    },
    Case {
        name: "fault-1-2h",
        covers: "IDV(1), a step in the A/C feed ratio, with the controllers \
                 pushing back against it",
        build: fault_1_2h,
        digest: 0x29ff_575a_59cf_fd30,
    },
    Case {
        name: "open-loop-4h",
        covers: "the valves held at their initial positions, so nothing damps \
                 the divergence: reactor pressure trips at 3.06 h and the \
                 plant then freezes and keeps reporting (teprob.f:807-811, \
                 delta D-007), which is the shutdown path no other case \
                 reaches",
        build: open_loop_4h,
        digest: 0x50e7_14b2_b49a_ce16,
    },
    Case {
        name: "rk4-1h",
        covers: "four derivative evaluations and a stage-weighted sum per \
                 step: four times the arithmetic, in a different order",
        build: rk4_1h,
        digest: 0xa6a6_77f8_c0a7_cebe,
    },
    Case {
        name: "dormand-prince-1h",
        covers: "seven stages and the densest floating-point path in the \
                 project, including the embedded error estimate",
        build: dormand_prince_1h,
        digest: 0x92a0_55e5_c4c2_5017,
    },
    Case {
        name: "baseline-9h",
        covers: "32,400 steps, past the driver's forced IDV(12) at eight \
                 hours (delta D-011), and the longest reach of the generator \
                 in this table",
        build: baseline_9h,
        digest: 0xce52_3c23_6499_1886,
    },
];

/// Run a scenario and return the digest of everything it emits.
///
/// Streaming: the samples are hashed as they are produced and never collected,
/// so the cost is independent of the run's length.
#[must_use]
pub fn digest(scenario: Scenario) -> u64 {
    let mut hash = Fnv1a64::new();
    // The outcome is deliberately not hashed. A trip freezes the plant and the
    // frozen samples keep coming (`teprob.f:807-811`, delta D-007), so a run
    // that ended differently already digests differently through its rows.
    let _ = Simulation::new(scenario).run_into(&mut hash);
    hash.finish()
}

/// The whole table as one number, for a log entry or a status line.
///
/// Each case contributes its name and its freshly computed digest, so a
/// rename, a reordering, an addition or a changed number all move it. Compare
/// against [`SUITE_DIGEST`].
#[must_use]
pub fn suite_digest() -> u64 {
    let mut hash = Fnv1a64::new();
    for case in CASES {
        hash.write_str(case.name);
        hash.write_u64(case.compute());
    }
    hash.finish()
}

/// The committed value of [`suite_digest`].
pub const SUITE_DIGEST: u64 = 0x9538_92e8_61fd_5e68;

/// A case whose digest is not the one committed.
///
/// This is the most important value this project can produce. It means two
/// platforms, or two builds, do not agree about what the simulator computes.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Mismatch {
    /// The case that disagreed, by [`Case::name`].
    pub case: &'static str,
    /// The committed digest.
    pub expected: u64,
    /// What this platform produced.
    pub computed: u64,
}

impl core::fmt::Display for Mismatch {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(
            f,
            "tier 9 case `{}`: expected {:016x}, this platform computed {:016x}",
            self.case, self.expected, self.computed
        )
    }
}

/// Check every case against its committed digest.
///
/// # Errors
///
/// The first [`Mismatch`], in table order. The first one is the informative
/// one: the cases share almost all of their code, so a difference in the
/// baseline explains the rest.
pub fn check() -> Result<(), Mismatch> {
    for case in CASES {
        let computed = case.compute();
        if computed != case.digest {
            return Err(Mismatch {
                case: case.name,
                expected: case.digest,
                computed,
            });
        }
    }
    Ok(())
}
