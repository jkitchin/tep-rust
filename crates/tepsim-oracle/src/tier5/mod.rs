//! Tier 5: the run harness.
//!
//! A *run* is one scenario at one seed, closed loop for a fixed horizon,
//! sampled at the cadence `temain_mod.f` itself samples at. B-0047a builds the
//! runs; B-0047b judges them.
//!
//! # What a sample is
//!
//! `temain_mod.f:401` writes output on `MOD(I, 180) == 0`, that is every 180
//! simulated seconds. Over 48 hours that is 960 rows, and each row is the 41
//! measurements and the 12 manipulated variables: the 53 variables `PLAN.org`
//! names for the correlation matrix.
//!
//! Sampling rather than keeping every step is not only a cost decision, though
//! it is that: it is what makes the data *comparable to the published files*,
//! which are at exactly this cadence. Tier 7 will compare against `d00`
//! through `d21` and would otherwise be comparing different things.
//!
//! # How many scenarios
//!
//! Twenty-one, not the twenty-two `PLAN.org` says. The vendored `teprob.f`
//! carries twenty disturbances: `teprob.f:340` loops `DO 500 I=1,20` and
//! `crate::FAULTS` has twenty entries. Later versions of the model add
//! `IDV(21)`, and the count in `PLAN.org` was written from that literature
//! rather than from this source. So: nominal plus `IDV(1)` through `IDV(20)`.
//!
//! # Cost
//!
//! Measured in release: 766 ms for one 48-hour port run of 172,800 steps. The
//! full battery, 21 scenarios by 100 seeds, is about 27 minutes per source.
//! That is right for `cargo xtask validate --tiers 5` and wrong for `ci`,
//! which is why the size is selected by an environment variable in the same
//! way Tier 1's sweep and Tier 4's horizon are.

pub mod battery;

use tepsim_control::{DRIVER_INITIAL_VALVES, PRESET, STEADY_STATE_STEPS};
use tepsim_core::TemperatureSeeds;

use crate::Oracle;

/// The integrator step, one second in hours. `temain_mod.f:231`.
pub const DT: f64 = 1.0 / 3600.0;

/// Output cadence: every 180 steps, as `temain_mod.f:401` does.
pub const SAMPLE_EVERY: usize = 180;

/// Variables recorded per sample: `XMEAS(1..41)` then `XMV(1..12)`.
pub const VARIABLES: usize = 53;

/// How many disturbances the vendored `teprob.f` has. See the module docs for
/// why this is twenty rather than twenty-one.
pub const FAULTS: usize = 20;

/// Scenarios: nominal plus one per disturbance.
pub const SCENARIOS: usize = FAULTS + 1;

/// `temain_mod.f:369-394`, transcribed. The order within each group is the
/// source's.
///
/// `temain_mod.f`'s main loop cannot be called: `instrument.rs` turns the
/// program into a subroutine nothing calls, precisely so it does not run. So
/// the loop is reproduced here by calling the same subroutines on the same
/// schedule. That is not circular, because the port's schedule is written
/// independently in [`tepsim_control::Scheme`] and the two are compared.
pub fn run_fortran_controllers(oracle: &mut Oracle, step: usize) {
    if step % 3 == 0 {
        oracle.contrl1();
        oracle.contrl2();
        oracle.contrl3();
        oracle.contrl4();
        oracle.contrl5();
        oracle.contrl6();
        oracle.contrl7();
        oracle.contrl8();
        oracle.contrl9();
        oracle.contrl10();
        oracle.contrl11();
        oracle.contrl16();
        oracle.contrl17();
        oracle.contrl18();
    }
    if step % 360 == 0 {
        oracle.contrl13();
        oracle.contrl14();
        oracle.contrl15();
        oracle.contrl19();
    }
    if step % 900 == 0 {
        oracle.contrl20();
    }
}

/// `CONSHAND`, `temain_mod.f:1401-1404`.
pub fn conshand(oracle: &mut Oracle) {
    let mut xmv = oracle.manipulated();
    // Not `f64::clamp`: `CONSHAND` tests `.LE. 0.0`, so it normalises negative
    // zero where `clamp`'s strict `<` leaves it alone. See
    // `tepsim_control::Scheme::clamp`.
    #[allow(clippy::manual_clamp, reason = "CONSHAND normalises negative zero")]
    for valve in xmv.iter_mut().take(11) {
        if *valve <= 0.0 {
            *valve = 0.0;
        }
        if *valve >= 100.0 {
            *valve = 100.0;
        }
    }
    oracle.set_manipulated(&xmv);
}

/// Load the Braatz preset into the oracle's `COMMON`.
pub fn load_preset(oracle: &mut Oracle) {
    let mut all = oracle.ctrlall();
    all.deltat = DT;
    for entry in &PRESET {
        all.setpt[entry.setpoint_index - 1] = entry.setpoint;
    }
    oracle.set_ctrlall(&all);
    oracle.set_flag6(0);

    // Written out rather than generated: twenty `COMMON` blocks are twenty
    // types, eight of them without a `taui` field, and that difference is
    // exactly what B-0037 made structural.
    macro_rules! p {
        ($n:literal) => {
            tepsim_control::preset($n).expect("a preset").tuning
        };
    }
    oracle.set_ctrl1(&crate::Ctrl1 {
        gain: p!(1).gain,
        errold: 0.0,
    });
    oracle.set_ctrl2(&crate::Ctrl2 {
        gain: p!(2).gain,
        errold: 0.0,
    });
    oracle.set_ctrl3(&crate::Ctrl3 {
        gain: p!(3).gain,
        errold: 0.0,
    });
    oracle.set_ctrl4(&crate::Ctrl4 {
        gain: p!(4).gain,
        errold: 0.0,
    });
    oracle.set_ctrl5(&crate::Ctrl5 {
        gain: p!(5).gain,
        taui: p!(5).reset.expect("PI"),
        errold: 0.0,
    });
    oracle.set_ctrl6(&crate::Ctrl6 {
        gain: p!(6).gain,
        errold: 0.0,
    });
    oracle.set_ctrl7(&crate::Ctrl7 {
        gain: p!(7).gain,
        errold: 0.0,
    });
    oracle.set_ctrl8(&crate::Ctrl8 {
        gain: p!(8).gain,
        errold: 0.0,
    });
    oracle.set_ctrl9(&crate::Ctrl9 {
        gain: p!(9).gain,
        errold: 0.0,
    });
    oracle.set_ctrl10(&crate::Ctrl10 {
        gain: p!(10).gain,
        taui: p!(10).reset.expect("PI"),
        errold: 0.0,
    });
    oracle.set_ctrl11(&crate::Ctrl11 {
        gain: p!(11).gain,
        taui: p!(11).reset.expect("PI"),
        errold: 0.0,
    });
    oracle.set_ctrl13(&crate::Ctrl13 {
        gain: p!(13).gain,
        taui: p!(13).reset.expect("PI"),
        errold: 0.0,
    });
    oracle.set_ctrl14(&crate::Ctrl14 {
        gain: p!(14).gain,
        taui: p!(14).reset.expect("PI"),
        errold: 0.0,
    });
    oracle.set_ctrl15(&crate::Ctrl15 {
        gain: p!(15).gain,
        taui: p!(15).reset.expect("PI"),
        errold: 0.0,
    });
    oracle.set_ctrl16(&crate::Ctrl16 {
        gain: p!(16).gain,
        taui: p!(16).reset.expect("PI"),
        errold: 0.0,
    });
    oracle.set_ctrl17(&crate::Ctrl17 {
        gain: p!(17).gain,
        taui: p!(17).reset.expect("PI"),
        errold: 0.0,
    });
    oracle.set_ctrl18(&crate::Ctrl18 {
        gain: p!(18).gain,
        taui: p!(18).reset.expect("PI"),
        errold: 0.0,
    });
    oracle.set_ctrl19(&crate::Ctrl19 {
        gain: p!(19).gain,
        taui: p!(19).reset.expect("PI"),
        errold: 0.0,
    });
    oracle.set_ctrl20(&crate::Ctrl20 {
        gain: p!(20).gain,
        taui: p!(20).reset.expect("PI"),
        errold: 0.0,
    });
    oracle.set_ctrl22(&crate::Ctrl22 {
        gain: p!(22).gain,
        taui: p!(22).reset.expect("PI"),
        errold: 0.0,
    });
}

/// One scenario: the nominal plant, or a single disturbance.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Scenario {
    /// One-based `IDV` index, or zero for the fault-free plant.
    pub fault: usize,
}

impl Scenario {
    /// The fault-free plant.
    pub const NOMINAL: Self = Self { fault: 0 };

    /// `IDV(n)`, one-based.
    ///
    /// # Panics
    ///
    /// If `n` is zero or above [`FAULTS`].
    #[must_use]
    pub const fn fault(n: usize) -> Self {
        assert!(n >= 1 && n <= FAULTS, "IDV index out of range");
        Self { fault: n }
    }

    /// Every scenario, nominal first.
    pub fn all() -> impl Iterator<Item = Self> {
        (0..=FAULTS).map(|fault| Self { fault })
    }

    /// The `IDV` vector this scenario asks for.
    #[must_use]
    pub fn disturbances(self) -> [f64; 20] {
        core::array::from_fn(|i| f64::from(u8::from(i + 1 == self.fault)))
    }

    /// The same, as the oracle's integer flags.
    #[must_use]
    pub fn disturbance_flags(self) -> [i32; 20] {
        core::array::from_fn(|i| i32::from(i + 1 == self.fault))
    }

    /// A short label for a report row.
    #[must_use]
    pub fn label(self) -> String {
        if self.fault == 0 {
            "nominal".into()
        } else {
            format!("IDV({})", self.fault)
        }
    }
}

/// The generator word for run `index`.
///
/// `TESUB7` is a multiplicative congruential generator, `G <- G * 9228907 mod
/// 2^32` (`teprob.f:1551`). Its period is maximal only for an odd seed, and
/// zero is a fixed point, so a seed has to be odd and non-zero. These are
/// produced by a SplitMix64 mix of the index, masked to 32 bits and forced
/// odd: unrelated to the TEP generator's own recurrence, so consecutive
/// indices do not give consecutive positions in one stream.
///
/// Index 0 is [`crate::golden::SEED`] itself, so the first run of every
/// battery is the one every other test in this repository uses.
#[must_use]
pub fn seed(index: usize) -> f64 {
    if index == 0 {
        return crate::golden::SEED;
    }
    let mut z = (index as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^= z >> 31;
    // Odd and inside 2^32.
    let word = (z & 0xFFFF_FFFF) | 1;
    word as f64
}

/// One run: a scenario at a seed, sampled.
#[derive(Clone, Debug, PartialEq)]
pub struct Run {
    /// Which scenario.
    pub scenario: Scenario,
    /// Which generator word it started from.
    pub seed: f64,
    /// `samples[k]` is the 53 variables at step `(k + 1) * SAMPLE_EVERY`.
    pub samples: Vec<[f64; VARIABLES]>,
    /// The step at which the plant tripped, if it did.
    ///
    /// A tripped run is *kept*, not discarded. `teprob.f:807-811` freezes the
    /// plant rather than stopping it, so the samples after a trip are constant
    /// and that is itself the behaviour under test. Discarding them would hide
    /// a port that tripped where the original did not.
    pub tripped: Option<usize>,
}

impl Run {
    /// One variable's series across the run.
    #[must_use]
    pub fn series(&self, variable: usize) -> Vec<f64> {
        self.samples.iter().map(|row| row[variable]).collect()
    }

    /// All 53 series.
    #[must_use]
    pub fn all_series(&self) -> Vec<Vec<f64>> {
        (0..VARIABLES).map(|v| self.series(v)).collect()
    }
}

/// The warm start and priming measurements a run begins from.
///
/// Taken from the oracle once, because `TEINIT` calls `TEFUNC` internally and
/// the four Newton warm starts it leaves are not the nominal literals; see
/// B-0034. [`Oracle::init_cold`] makes this independent of what ran before.
#[derive(Clone, Copy, Debug)]
pub struct Start {
    /// The four vessel temperatures after `TEINIT`.
    pub seeds: TemperatureSeeds,
    /// `XMEAS` after `TEINIT`, which the driver's first controller fire reads.
    pub measurements: [f64; 41],
    /// `YY` after `TEINIT`.
    pub state: [f64; 50],
}

/// Capture the common starting point for both sources.
pub fn start(oracle: &mut Oracle) -> Start {
    let (_, state) = oracle.init_cold();
    oracle.set_disturbances(&[0; 20]);
    let common = oracle.teproc();
    Start {
        seeds: TemperatureSeeds {
            reactor: common.tcr,
            separator: common.tcs,
            stripper: common.tcc,
            mixing: common.tcv,
        },
        measurements: oracle.measurements(),
        state,
    }
}

/// Run the port, closed loop, and sample it.
///
/// Delegates to [`tepsim::Simulation`]. There was a second copy of this loop
/// here until B-0052, written before the facade existed; two transcriptions of
/// `temain_mod.f`'s main loop is one more than the project can keep honest, and
/// `facade_equivalence.rs` is what proved they were the same before this one
/// was deleted.
///
/// `start` is accepted for symmetry with [`run_fortran`] and for the
/// generator word; the warm start comes from
/// [`tepsim_core::TemperatureSeeds::after_initialisation`], which
/// `facade_equivalence.rs` asserts equals `start.seeds` bit for bit.
#[must_use]
pub fn run_port(start: &Start, scenario: Scenario, seed: f64, hours: usize) -> Run {
    let _ = start;
    // `faithful`, never `baseline`. Every differential in Tiers 5, 6 and 7 runs
    // through here, and each one asks whether this port *is* `teprob.f`. It
    // cannot answer that while carrying a deliberate difference from it, so the
    // two Class C quirks the sign-off of 2026-08-28 turned off by default are
    // turned back on: the plant freezes on a trip rather than the run ending
    // (D-007), and the driver forces `IDV(12)` at eight hours (D-011). The
    // Fortran on the other side of the comparison does both unconditionally.
    let mut described = tepsim::Scenario::faithful()
        .with_hours(hours as f64)
        .with_seed(seed);
    if scenario.fault != 0 {
        described = described.with_fault(scenario.fault);
    }
    described.sample_every = SAMPLE_EVERY;

    let finished = tepsim::Simulation::new(described).run();
    Run {
        scenario,
        seed,
        samples: finished.samples.iter().map(tepsim::Sample::row).collect(),
        tripped: match finished.outcome {
            tepsim::Outcome::Tripped { step, .. } | tepsim::Outcome::SolveFailed { step } => {
                Some(step)
            }
            tepsim::Outcome::Completed => None,
        },
    }
}

/// Run the Fortran, closed loop, and sample it the same way.
#[must_use]
pub fn run_fortran(
    oracle: &mut Oracle,
    start: &Start,
    scenario: Scenario,
    seed: f64,
    hours: usize,
) -> Run {
    let (_, _) = oracle.init_cold();
    oracle.set_disturbances(&scenario.disturbance_flags());
    oracle.set_rng(seed);
    load_preset(oracle);
    oracle.set_manipulated(&DRIVER_INITIAL_VALVES);

    let mut yy = start.state;
    let mut samples = Vec::with_capacity(hours * 3_600 / SAMPLE_EVERY);
    let mut tripped = None;
    let mut t = 0.0;

    for step in 1..=hours * 3_600 {
        // `temain_mod.f:366-368`, the driver's forced IDV(12), on top of
        // whatever the scenario asked for.
        if step >= STEADY_STATE_STEPS {
            let mut idv = scenario.disturbance_flags();
            idv[11] = 1;
            oracle.set_disturbances(&idv);
        }
        run_fortran_controllers(oracle, step);
        let yp = oracle.derivatives(t, &yy);
        if oracle.shutdown_flag() != 0 && tripped.is_none() {
            tripped = Some(step);
        }
        // Read out only on a sample step. `XMEAS` and `XMV` are read from
        // `COMMON` across the FFI boundary, and at 172,800 steps per run that
        // is 1.7 million array copies of which 959 in 960 are discarded.
        // `CONSHAND` still runs every step, because it *writes*.
        let sampling = step % SAMPLE_EVERY == 0;
        let measurements = if sampling {
            oracle.measurements()
        } else {
            [0.0; 41]
        };
        conshand(oracle);

        for (slot, rate) in yy.iter_mut().zip(yp) {
            // Two roundings, not a fused multiply-add; see `tier2`.
            #[allow(clippy::suboptimal_flops, reason = "matches INTGTR's rounding")]
            {
                *slot += DT * rate;
            }
        }
        t += DT;

        if sampling {
            samples.push(row(&measurements, &oracle.manipulated()));
        }
    }

    Run {
        scenario,
        seed,
        samples,
        tripped,
    }
}

/// `XMEAS(1..41)` then `XMV(1..12)`, in that order.
fn row(measurements: &[f64; 41], valves: &[f64; 12]) -> [f64; VARIABLES] {
    let mut out = [0.0; VARIABLES];
    out[..41].copy_from_slice(measurements);
    out[41..].copy_from_slice(valves);
    out
}

/// How big a battery to run, from `TEP_TIER5`.
///
/// `full` is the gate: 21 scenarios, 100 seeds, 48 hours, about 27 minutes per
/// source. Anything else, including the variable being absent, is the smoke
/// battery that `ci` runs.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Battery {
    /// How many scenarios, nominal first.
    pub scenarios: usize,
    /// How many seeds per scenario.
    pub seeds: usize,
    /// How many simulated hours per run.
    pub hours: usize,
}

impl Battery {
    /// The environment variable that selects between [`Battery::SMOKE`] and
    /// [`Battery::FULL`].
    pub const ENV: &'static str = "TEP_TIER5";

    /// What `ci` runs: enough to exercise every code path, not enough to be a
    /// statistical claim.
    pub const SMOKE: Self = Self {
        scenarios: 3,
        seeds: 4,
        hours: 2,
    };

    /// What `cargo xtask validate --tiers 5` runs.
    pub const FULL: Self = Self {
        scenarios: SCENARIOS,
        seeds: 100,
        hours: 48,
    };

    /// Which battery to run.
    #[must_use]
    pub fn selected() -> Self {
        match std::env::var(Self::ENV).as_deref() {
            Ok("full") => Self::FULL,
            _ => Self::SMOKE,
        }
    }

    /// How many runs this battery is, per source.
    #[must_use]
    pub const fn runs(&self) -> usize {
        self.scenarios * self.seeds
    }

    /// How many samples each run produces.
    #[must_use]
    pub const fn samples(&self) -> usize {
        self.hours * 3_600 / SAMPLE_EVERY
    }
}

/// The generator words `teprob.f` carries in comments, one per published
/// dataset.
///
/// `teprob.f:1187-1256`. The compiled-in seed is `4651207995`; everything else
/// in that block is commented out and names the file it was used to generate.
/// They are recorded here because B-0051 (Tier 7) has to reproduce `d00`
/// through `d21` and these are the only statement anywhere of how they were
/// made.
///
/// Three things to know before using them.
///
/// **Many exceed 2^32.** `TESUB7` computes `G <- DMOD(G * 9228907, 2^32)`, so
/// the first draw reduces the seed regardless; a value above the modulus is a
/// perfectly ordinary starting point and not a transcription error.
///
/// **Some are even.** `4346024432`, `6678322168` and others. A multiplicative
/// generator modulo a power of two keeps the factors of two it starts with, so
/// an even seed has a shorter period and its low bits never move. That is what
/// the original did, and reproducing the published files means doing the same.
///
/// **One does not compile.** `d19_te` is written `G=9090909232.DO`, with the
/// letter O rather than a zero. It is inside a comment, so nothing ever caught
/// it. Recorded here as `9090909232` on the assumption that the digit was
/// meant, and flagged so that a Tier 7 failure on `d19` has an obvious first
/// suspect.
pub mod published_seeds {
    /// The seed actually compiled into `teprob.f:1187`.
    pub const COMPILED_IN: f64 = 4_651_207_995.0;

    /// Labelled "original" at `teprob.f:1189`.
    pub const ORIGINAL: f64 = 1_431_655_765.0;

    /// Labelled "d00_tr_new" at `teprob.f:1188`.
    pub const D00_TRAINING_NEW: f64 = 5_687_912_315.0;

    /// `d00_tr` through `d26_tr`, in order.
    pub const TRAINING: [f64; 27] = [
        4_243_534_565.0,
        7_854_912_354.0,
        3_456_432_354.0,
        1_731_738_903.0,
        4_346_024_432.0,
        5_784_921_734.0,
        6_678_322_168.0,
        7_984_782_901.0,
        8_934_302_332.0,
        9_873_223_412.0,
        1_089_278_833.0,
        1_940_284_333.0,
        2_589_274_931.0,
        3_485_834_345.0,
        4_593_493_842.0,
        5_683_213_434.0,
        6_788_343_442.0,
        1_723_234_455.0,
        8_943_243_993.0,
        9_445_382_439.0,
        9_902_234_324.0,
        2_144_342_545.0,
        3_433_249_064.0,
        4_356_565_463.0,
        8_998_485_332.0,
        7_654_534_567.0,
        5_457_789_234.0,
    ];

    /// An alternative `d18_tr`, labelled `dd18_tr` at `teprob.f:1209`. Why
    /// there are two is not stated anywhere in the source.
    pub const D18_TRAINING_ALTERNATIVE: f64 = 1_234_567_890.0;

    /// `d00_te` through `d26_te`, in order.
    ///
    /// Entry 19 is the one written `.DO` in the source; see the module docs.
    pub const TESTING: [f64; 27] = [
        1_254_545_354.0,
        2_994_833_239.0,
        2_891_123_453.0,
        3_420_494_299.0,
        4_598_956_239.0,
        5_658_678_765.0,
        6_598_593_453.0,
        7_327_843_434.0,
        8_943_242_344.0,
        9_343_430_004.0,
        1_039_839_281.0,
        1_134_345_551.0,
        2_232_323_236.0,
        3_454_354_353.0,
        4_545_445_883.0,
        5_849_489_384.0,
        6_284_545_932.0,
        4_342_232_344.0,
        5_635_346_588.0,
        9_090_909_232.0,
        8_322_308_324.0,
        2_132_432_423.0,
        5_454_589_923.0,
        6_923_255_678.0,
        8_493_323_434.0,
        9_338_398_429.0,
        1_997_072_199.0,
    ];

    /// The index in [`TESTING`] whose source line is malformed.
    pub const MALFORMED_TESTING_INDEX: usize = 19;
}
