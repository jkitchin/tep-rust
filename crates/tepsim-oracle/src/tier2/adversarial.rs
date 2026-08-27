//! The third Tier 2 sampling pool: states placed deliberately at a
//! discontinuity or a clamp.
//!
//! Tier 1 made the case for this pool empirically. Its boundary generator was
//! an eighth of the sweep and produced most of the worst cases, because the
//! interesting numerics live where the model changes its mind, not in the
//! interior. Every branch in `TEFUNC` is enumerable from the listing, so this
//! pool is built by hand rather than sampled.
//!
//! # Constructed by bisection, not by algebra
//!
//! Placing a state at `VLR/7.8 = 50` means inverting the model: the reactor
//! liquid volume is the total liquid holdup over a density that itself depends
//! on the composition and the temperature that the holdups determine. Doing
//! that in closed form would be a second implementation of the very code under
//! test, and a wrong one would put the state somewhere plausible but not on the
//! boundary, which is worse than having no state at all because it *looks* like
//! coverage.
//!
//! So the oracle is used as the forward map. A [`Knob`] perturbs one scalar
//! degree of freedom of a scenario, a [`Target`] reads one number back out of
//! the resulting snapshot, and [`seek`] bisects until the number lands.
//! Whatever the model does, the state ends up where it was asked to be, and
//! [`Boundary::verify`] then confirms it rather than assuming it.
//!
//! # One branch in the catalogue is unreachable, and that is a finding
//!
//! `teprob.f:585-586` clamps the purge flow when `PTS` falls below 760 mmHg.
//! No state built from a nominal separator composition reaches it. Cooling the
//! separator to 0 C and removing its entire vapour inventory leaves `PTS` at
//! *811.49 mmHg*, because the separator's own liquid exerts that much vapour
//! pressure at 0 C and the ideal-gas term floors at `TKS = 273.15` rather than
//! at zero. The clamp is 7% away and cannot be crossed from that side.
//!
//! It is therefore not in [`catalogue`]. The port must still implement it, and
//! B-0021 covers it with a unit test at a composition chosen to make it fire,
//! rather than with a physically-reachable state. The measurement is pinned by
//! a test so that this stays a known fact rather than a forgotten omission.
//!
//! # A boundary state does not exercise the branch it bounds
//!
//! Every comparison in `TEFUNC` is strict: `.GT.`, `.LT.`. A state placed
//! exactly on a threshold therefore takes the *other* side. That is valuable
//! on its own, because misreading `.GT.` as `.GE.` is invisible anywhere else,
//! but it means the branch behind the threshold stays unexercised unless a
//! second state is placed past it.
//!
//! Two pairs in this catalogue are built that way, and both were found the
//! same way: by a coverage test in the item that implemented the branch,
//! failing because the branch had never been entered. `PR = CPPRMX` pairs with
//! "PR above CPPRMX", and `TCC = 5.292` pairs with "TCC below". Any future
//! threshold should get both from the start.
//!
//! # The compressor clamps arrived late, on purpose
//!
//! `teprob.f:590-591` clamps the compressor pressure ratio at 1 and at
//! `CPPRMX`. Both were left out of the original catalogue because `CPPRMX` is
//! a compressor constant and the compressor did not exist yet; B-0021 added
//! them. The `CPPRMX` target aims at the *single-precision* value
//! 1.2999999523162842, not at 1.3: the comparison at line 591 is against what
//! gfortran stored, and a state placed at the double literal sits on the wrong
//! side of it.
//!
//! # Knobs that hold the temperature fixed
//!
//! Scaling a vessel's liquid holdups alone would change its specific energy,
//! `ESR = ETR / UTLR` (`teprob.f:456`), and so move its temperature as a side
//! effect. The volume knobs therefore scale the holdups *and* the energy state
//! together, which leaves the specific energy, and so the temperature,
//! untouched. The temperature knobs move the energy alone, for the same reason
//! in reverse.

use tepsim_core::state::N_STATES;

use super::{Pool, Scenario, Snapshot};
use crate::Oracle;

/// One scalar degree of freedom of a scenario.
#[derive(Clone, Copy, Debug)]
pub struct Knob {
    /// What it moves, for the report.
    pub name: &'static str,
    /// Applies the knob at setting `k` to a fresh copy of the base scenario.
    pub apply: fn(&mut Scenario, f64),
}

impl Knob {
    /// Scale a vessel's liquid holdups and its energy together, leaving the
    /// specific energy and therefore the temperature unchanged.
    fn scale_holdup(
        scenario: &mut Scenario,
        slots: core::ops::Range<usize>,
        energy: usize,
        k: f64,
    ) {
        for slot in slots {
            scenario.state[slot] *= k;
        }
        scenario.state[energy] *= k;
    }

    /// Reactor liquid inventory: states 4-8 and the energy at 9.
    pub const REACTOR_VOLUME: Self = Self {
        name: "reactor liquid inventory",
        apply: |s, k| Self::scale_holdup(s, 3..8, 8, k),
    };

    /// Separator liquid inventory: states 13-17 and the energy at 18.
    pub const SEPARATOR_VOLUME: Self = Self {
        name: "separator liquid inventory",
        apply: |s, k| Self::scale_holdup(s, 12..17, 17, k),
    };

    /// Stripper liquid inventory: states 19-26 and the energy at 27.
    pub const STRIPPER_VOLUME: Self = Self {
        name: "stripper liquid inventory",
        apply: |s, k| Self::scale_holdup(s, 18..26, 26, k),
    };

    /// Stripper energy alone, which moves `TCC` without moving the inventory.
    pub const STRIPPER_ENERGY: Self = Self {
        name: "stripper specific energy",
        apply: |s, k| s.state[26] *= k,
    };

    /// Reactor energy alone, which moves `TCR`.
    pub const REACTOR_ENERGY: Self = Self {
        name: "reactor specific energy",
        apply: |s, k| s.state[8] *= k,
    };

    /// Mixing-zone inventory, which moves `PTV` (`teprob.f:491`).
    pub const MIXING_INVENTORY: Self = Self {
        name: "mixing zone inventory",
        apply: |s, k| {
            for slot in 27..35 {
                s.state[slot] *= k;
            }
            s.state[35] *= k;
        },
    };

    /// Reactor vapour holdups, states 1-3, which set the ideal-gas part of
    /// `PTR` (`teprob.f:479`).
    pub const REACTOR_VAPOUR: Self = Self {
        name: "reactor vapour inventory",
        apply: |s, k| {
            for slot in 0..3 {
                s.state[slot] *= k;
            }
        },
    };

    /// Separator vapour holdups, states 10-12, which set the ideal-gas part of
    /// `PTS` (`teprob.f:481`).
    pub const SEPARATOR_VAPOUR: Self = Self {
        name: "separator vapour inventory",
        apply: |s, k| {
            for slot in 9..12 {
                s.state[slot] *= k;
            }
        },
    };

    /// Separator energy alone, which moves `TCS` and so the Antoine part of
    /// `PTS` (`teprob.f:488`).
    pub const SEPARATOR_ENERGY: Self = Self {
        name: "separator specific energy",
        apply: |s, k| s.state[17] *= k,
    };

    /// Separator vapour holdup and energy together, which is the most `PTS` can
    /// be driven down by.
    ///
    /// Used to establish that the purge clamp at `teprob.f:585-586` is *not*
    /// reachable rather than to reach it; see
    /// `tests/tier2_adversarial.rs`.
    pub const SEPARATOR_DEPRESSURISE: Self = Self {
        name: "separator vapour and energy",
        apply: |s, k| {
            for slot in 9..12 {
                s.state[slot] *= k;
            }
            s.state[17] *= k;
        },
    };

    /// The separator underflow valve position, `VPOS(7)`, which sets `FTM(11)`
    /// (`teprob.f:571`).
    pub const UNDERFLOW_VALVE: Self = Self {
        name: "separator underflow valve",
        apply: |s, k| s.state[44] = k,
    };
}

/// One number read back out of an evaluation, and the boundary it must reach.
#[derive(Clone, Copy, Debug)]
pub struct Target {
    /// What is being placed, for the report.
    pub name: &'static str,
    /// The branch this boundary belongs to, and where to read it.
    pub why: &'static str,
    /// Reads the observable out of a snapshot.
    pub observe: fn(&Snapshot) -> f64,
    /// The value the observable must take.
    pub value: f64,
    /// The magnitude to judge the miss against.
    ///
    /// Usually the target value itself, but several boundaries are *pressure
    /// differences* whose boundary value is zero, and a relative error against
    /// zero is not a number. For those the scale is the pressure itself.
    pub scale: f64,
}

/// A constructed adversarial state, with the evidence that it landed.
#[derive(Clone, Debug)]
pub struct Boundary {
    /// What it was built for.
    pub target: Target,
    /// The knob setting that got there.
    pub setting: f64,
    /// The scenario itself.
    pub scenario: Scenario,
    /// What the observable actually came out at.
    pub reached: f64,
    /// Whether evaluating it trips the shutdown, in which case its derivative
    /// is all zeros (`teprob.f:807-811`) and it is coverage, not evidence.
    pub tripped: bool,
}

impl Boundary {
    /// How far the constructed state sits from the boundary it was aimed at,
    /// relative to the boundary's own magnitude.
    #[must_use]
    pub fn miss(&self) -> f64 {
        (self.reached - self.target.value).abs() / self.target.scale.abs()
    }

    /// Fail unless the state actually sits on its boundary.
    ///
    /// A state built for `TCC = 5.292` that came out a degree away is worse
    /// than no state at all: it exercises the ordinary branch while appearing
    /// in the report as boundary coverage.
    ///
    /// # Panics
    ///
    /// If the observable missed its target by more than `tolerance`, relative.
    pub fn verify(&self, tolerance: f64) {
        assert!(
            self.miss() <= tolerance,
            "{}: aimed at {} = {}, reached {} (off by {:e} relative). {}",
            self.target.name,
            self.target.name,
            self.target.value,
            self.reached,
            self.miss(),
            self.target.why
        );
    }
}

/// Bisect `knob` until `target`'s observable reaches its value.
///
/// Requires the observable to be monotone in the knob over `bracket`, which is
/// checked: if the two ends do not straddle the target, this returns `None`
/// rather than converging to something arbitrary. A silently mis-placed
/// adversarial state is the one failure this pool cannot afford.
///
/// Sixty iterations, which is enough to exhaust an `f64`'s worth of bracket
/// regardless of its width.
pub fn seek(
    oracle: &mut Oracle,
    base: &Scenario,
    knob: Knob,
    target: Target,
    bracket: (f64, f64),
) -> Option<Boundary> {
    let evaluate = |oracle: &mut Oracle, k: f64| -> (f64, Snapshot) {
        let mut scenario = base.clone();
        (knob.apply)(&mut scenario, k);
        let snapshot = scenario.force(oracle);
        ((target.observe)(&snapshot), snapshot)
    };

    let (mut low, mut high) = bracket;
    let (at_low, _) = evaluate(oracle, low);
    let (at_high, _) = evaluate(oracle, high);
    if !at_low.is_finite() || !at_high.is_finite() {
        return None;
    }
    // The target must lie between the ends, in whichever direction they run.
    let straddles = (at_low - target.value) * (at_high - target.value) <= 0.0;
    if !straddles {
        return None;
    }
    let increasing = at_high > at_low;

    for _ in 0..60 {
        let mid = 0.5 * (low + high);
        let (value, _) = evaluate(oracle, mid);
        if !value.is_finite() {
            return None;
        }
        if (value < target.value) == increasing {
            low = mid;
        } else {
            high = mid;
        }
    }

    let setting = 0.5 * (low + high);
    let mut scenario = base.clone();
    (knob.apply)(&mut scenario, setting);
    let snapshot = scenario.force(oracle);
    Some(Boundary {
        target,
        setting,
        reached: (target.observe)(&snapshot),
        tripped: snapshot.tripped,
        scenario,
    })
}

/// Every boundary in the catalogue, with the knob and bracket that reaches it.
///
/// Each entry names the `teprob.f` line whose branch it sits on. The list is
/// the one enumerated in `PLAN.org`, "Tier 2", plus the shutdown conditions.
#[must_use]
pub fn catalogue() -> Vec<(Knob, Target, (f64, f64))> {
    // The heat-transfer ramp: `teprob.f:663-668` switches at VLR/7.8 = 10 and
    // = 50, so the boundary values of VLR are 78 and 390.
    let vlr = |s: &Snapshot| s.common.vlr;
    let vls = |s: &Snapshot| s.common.vls;
    let vlc = |s: &Snapshot| s.common.vlc;
    let tcc = |s: &Snapshot| s.common.tcc;
    let tcr = |s: &Snapshot| s.common.tcr;
    let ftm11 = |s: &Snapshot| s.common.ftm[10];
    let ptr = |s: &Snapshot| s.common.ptr;

    vec![
        (
            Knob::REACTOR_VOLUME,
            Target {
                name: "VLR at the 10% heat-transfer breakpoint",
                why: "teprob.f:665, UARLEV pins to 0 below VLR/7.8 = 10",
                observe: vlr,
                value: 78.0,
                scale: 78.0,
            },
            (0.1, 3.0),
        ),
        (
            Knob::REACTOR_VOLUME,
            Target {
                name: "VLR at the 50% heat-transfer breakpoint",
                why: "teprob.f:663, UARLEV pins to 1 above VLR/7.8 = 50",
                observe: vlr,
                value: 390.0,
                scale: 390.0,
            },
            (0.1, 8.0),
        ),
        (
            Knob::STRIPPER_ENERGY,
            Target {
                name: "TCC at the lower stripping-factor branch",
                why: "teprob.f:617, TMPFAC pins to 0.1 below 5.292 C",
                observe: tcc,
                value: 5.292,
                scale: 5.292,
            },
            (0.02, 2.0),
        ),
        // And one *below* it. `teprob.f:617` tests `.LT.`, so the boundary
        // state above sits on 5.292 and takes the hyperbolic branch; it proves
        // the comparison is strict and leaves the pinned branch unexercised.
        // The same pairing as the `CPPRMX` clamp: for a strict comparison, the
        // boundary and the taken side are two different states.
        (
            Knob::STRIPPER_ENERGY,
            Target {
                name: "TCC below the lower stripping-factor branch",
                why: "teprob.f:618, the pinned branch the boundary state does not take",
                observe: tcc,
                value: 4.0,
                scale: 4.0,
            },
            (0.01, 2.0),
        ),
        (
            Knob::STRIPPER_ENERGY,
            Target {
                name: "TCC at the upper stripping-factor branch",
                why: "teprob.f:615, TMPFAC goes linear above 170 C",
                observe: tcc,
                value: 170.0,
                scale: 170.0,
            },
            (0.5, 6.0),
        ),
        (
            Knob::STRIPPER_ENERGY,
            Target {
                name: "TCC approaching the 177 C pole",
                why: "teprob.f:620, TMPFAC = 363.744/(177 - TCC) is singular",
                observe: tcc,
                value: 176.0,
                scale: 176.0,
            },
            (0.5, 6.0),
        ),
        (
            Knob::REACTOR_ENERGY,
            Target {
                name: "TCR at the shutdown limit",
                why: "teprob.f:706, ISD fires above 175 C",
                observe: tcr,
                value: 175.0,
                scale: 175.0,
            },
            (0.5, 3.0),
        ),
        (
            Knob::UNDERFLOW_VALVE,
            Target {
                name: "FTM(11) at the stripping-factor threshold",
                why: "teprob.f:614, the whole SFR block switches at 0.1",
                observe: ftm11,
                value: 0.1,
                scale: 0.1,
            },
            (0.0, 5.0),
        ),
        (
            Knob::REACTOR_VOLUME,
            Target {
                name: "VLR at the upper shutdown limit",
                why: "teprob.f:704, ISD fires above VLR/35.3145 = 24",
                observe: vlr,
                value: 24.0 * 35.3145,
                scale: 24.0 * 35.3145,
            },
            (0.1, 12.0),
        ),
        (
            Knob::REACTOR_VOLUME,
            Target {
                name: "VLR at the lower shutdown limit",
                why: "teprob.f:705, ISD fires below VLR/35.3145 = 2",
                observe: vlr,
                value: 2.0 * 35.3145,
                scale: 2.0 * 35.3145,
            },
            (0.05, 3.0),
        ),
        (
            Knob::SEPARATOR_VOLUME,
            Target {
                name: "VLS at the upper shutdown limit",
                why: "teprob.f:707, ISD fires above VLS/35.3145 = 12",
                observe: vls,
                value: 12.0 * 35.3145,
                scale: 12.0 * 35.3145,
            },
            (0.1, 20.0),
        ),
        (
            Knob::SEPARATOR_VOLUME,
            Target {
                name: "VLS at the lower shutdown limit",
                why: "teprob.f:708, ISD fires below VLS/35.3145 = 1",
                observe: vls,
                value: 35.3145,
                scale: 35.3145,
            },
            (0.02, 3.0),
        ),
        (
            Knob::STRIPPER_VOLUME,
            Target {
                name: "VLC at the upper shutdown limit",
                why: "teprob.f:709, ISD fires above VLC/35.3145 = 8",
                observe: vlc,
                value: 8.0 * 35.3145,
                scale: 8.0 * 35.3145,
            },
            (0.1, 20.0),
        ),
        (
            Knob::STRIPPER_VOLUME,
            Target {
                name: "VLC at the lower shutdown limit",
                why: "teprob.f:710, ISD fires below VLC/35.3145 = 1",
                observe: vlc,
                value: 35.3145,
                scale: 35.3145,
            },
            (0.02, 3.0),
        ),
        (
            Knob::REACTOR_VAPOUR,
            Target {
                name: "PTR at the reactor pressure shutdown limit",
                why: "teprob.f:703, ISD fires above XMEAS(7) = 3000 kPa gauge",
                observe: ptr,
                value: 3000.0 / 101.325 * 760.0 + 760.0,
                scale: 3000.0 / 101.325 * 760.0 + 760.0,
            },
            (0.5, 8.0),
        ),
        (
            Knob::MIXING_INVENTORY,
            Target {
                name: "PTV = PTR, the mixing-to-reactor flow clamp",
                why: "teprob.f:576-577, DLP is clamped at zero below this",
                observe: |s| s.common.ptv - s.common.ptr,
                value: 0.0,
                scale: 21_000.0,
            },
            (0.1, 2.0),
        ),
        (
            Knob::SEPARATOR_VAPOUR,
            Target {
                name: "PTR = PTS, the reactor-to-separator flow clamp",
                why: "teprob.f:580-581, DLP is clamped at zero below this",
                observe: |s| s.common.ptr - s.common.pts,
                value: 0.0,
                scale: 21_000.0,
            },
            (0.5, 4.0),
        ),
        (
            Knob::MIXING_INVENTORY,
            Target {
                name: "PTV = PTS, the recycle flow clamp",
                why: "teprob.f:596-597, DLP is clamped at zero below this",
                observe: |s| s.common.ptv - s.common.pts,
                value: 0.0,
                scale: 21_000.0,
            },
            (0.1, 2.0),
        ),
        // The two compressor pressure-ratio clamps, deferred here from B-0016
        // because they are branches against `CPPRMX`, which does not exist
        // until the compressor curve does (B-0021).
        //
        // The first lands on the same state as "PTV = PTS" above: a ratio of
        // one *is* a zero pressure difference. It is listed separately anyway,
        // because the catalogue is how a reader finds the state that covers a
        // given line, and `teprob.f:590` is not `teprob.f:597`.
        (
            Knob::MIXING_INVENTORY,
            Target {
                name: "PR = 1, the compressor reverse-flow clamp",
                why: "teprob.f:590, PR is clamped up to 1 below this",
                observe: |s| s.common.ptv / s.common.pts,
                value: 1.0,
                scale: 1.0,
            },
            (0.1, 2.0),
        ),
        (
            Knob::MIXING_INVENTORY,
            Target {
                name: "PR = CPPRMX, the compressor maximum-ratio clamp",
                why: "teprob.f:591, PR is clamped down to CPPRMX above this",
                observe: |s| s.common.ptv / s.common.pts,
                // Single precision, so 1.2999999523162842 rather than 1.3.
                // Aiming at the double literal would place the state on the
                // wrong side of the comparison.
                value: tepsim_core::flows::MAX_PRESSURE_RATIO,
                scale: tepsim_core::flows::MAX_PRESSURE_RATIO,
            },
            (0.5, 6.0),
        ),
        // And one state *inside* the clamped region. The boundary state above
        // sits exactly on `CPPRMX`, and `teprob.f:591` tests `.GT.`, so it
        // does not clamp: it proves the comparison is strict, which is the
        // subtle half. It leaves the taken branch unexercised, and a branch no
        // sampled state enters is not evidence, so this entry covers it.
        (
            Knob::MIXING_INVENTORY,
            Target {
                name: "PR above CPPRMX, inside the clamped region",
                why: "teprob.f:591, the branch the boundary state does not take",
                observe: |s| s.common.ptv / s.common.pts,
                value: 1.4,
                scale: 1.4,
            },
            (0.5, 6.0),
        ),
    ]
}

/// Build every boundary the catalogue can reach from `base`.
///
/// Returns the ones that landed and the names of the ones that did not, so a
/// gap in the coverage is reported rather than silently absent.
pub fn build(oracle: &mut Oracle, base: &Scenario) -> (Vec<Boundary>, Vec<&'static str>) {
    let mut built = Vec::new();
    let mut missed = Vec::new();
    for (knob, target, bracket) in catalogue() {
        match seek(oracle, base, knob, target, bracket) {
            Some(boundary) => built.push(boundary),
            None => missed.push(target.name),
        }
    }
    (built, missed)
}

/// A boundary's scenario, tagged as coming from the adversarial pool.
#[must_use]
pub fn pool() -> Pool {
    Pool::Adversarial
}

/// Every state slot the knobs in this module touch, for a coverage check.
#[must_use]
pub fn touched_slots() -> Vec<usize> {
    let mut slots: Vec<usize> = (0..9).chain(9..18).chain(18..27).chain(27..36).collect();
    slots.push(44);
    slots.retain(|s| *s < N_STATES);
    slots.sort_unstable();
    slots.dedup();
    slots
}
