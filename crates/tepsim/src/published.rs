//! The published `d00`-`d21` datasets: their shape, their seeds, and how to
//! generate more data in the same form.
//!
//! Thirty years of fault-detection papers were written against forty-four
//! files. This module is what a caller needs to write a forty-fifth: the run
//! geometry, the seed each file was made with, and the column order the files
//! are actually stored in.
//!
//! # This module does not claim to reproduce the published bytes
//!
//! It cannot, and `crates/tepsim-oracle/src/tier7.rs` measures how far off it
//! is. Four things stand in the way, none of them fixable by trying harder:
//! the generating toolchain is unrecorded and its `exp` is not this one's; the
//! output was rounded to five significant figures before anyone saw it; `d21`
//! names a disturbance this revision of `teprob.f` does not contain; and the
//! protocol behind the twenty-two *training* files is written down nowhere.
//! See [`Unavailable`] and [`TRAINING_IS_A_HYPOTHESIS`].
//!
//! What this module does claim is that a run built here has the same
//! *geometry* as the file it is named after: the same number of rows at the
//! same spacing, with the disturbance arriving at the same row.
//!
//! # What the driver states outright
//!
//! | Fact | `temain_mod.f` | Value |
//! |---|---|---|
//! | Run length | line 220 | `NPTS = 172800` steps, 48 h at one second |
//! | Settling time | line 226 | `SSPTS = 3600 * 8`, eight hours |
//! | Output cadence | line 401 | `MOD(I,180)`, so 960 rows |
//! | Output precision | line 1358 | `FORMAT(1X,E13.5)`, five significant digits |
//!
//! The seeds are the commented block at `teprob.f:1187-1256`.
//! `tests/published.rs` reads both files and asserts every constant below
//! against them, so none of this is transcribed on trust.

use alloc::format;
use alloc::string::String;

use tepsim_scenario::Event;

use crate::run::{MANIPULATED, MEASUREMENTS, Sample};
use crate::scenario::Scenario;

/// Steps between recorded rows, `MOD(I,180)` at `temain_mod.f:401`.
///
/// At the one-second step that is three minutes.
pub const SAMPLE_EVERY: usize = 180;

/// `NPTS`, `temain_mod.f:220`: the run length the driver ships with.
pub const DRIVER_STEPS: usize = 172_800;

/// `SSPTS`, `temain_mod.f:226`: how long the driver runs before the fault.
pub const SETTLING_STEPS: usize = 3600 * 8;

/// Rows in a `dNN_te.dat` file.
pub const TESTING_ROWS: usize = DRIVER_STEPS / SAMPLE_EVERY;

/// Rows in a `dNN.dat` training file, for `NN >= 1`.
pub const TRAINING_ROWS: usize = 480;

/// Rows in `d00.dat`, the only training file that is not [`TRAINING_ROWS`].
///
/// Upstream's own `reference/data/README.md` says 480 for this file too. It is
/// wrong, and those twenty extra rows are the only real evidence about how the
/// training files were made. See [`TRAINING_IS_A_HYPOTHESIS`].
pub const D00_TRAINING_ROWS: usize = 500;

/// Simulated hours in a testing run.
pub const TESTING_HOURS: f64 = 48.0;

/// When the disturbance arrives in a testing run.
pub const TESTING_FAULT_HOURS: f64 = 8.0;

/// Simulated hours a training run is generated for, before discarding.
///
/// Part of the hypothesis, not of the record. See [`TRAINING_IS_A_HYPOTHESIS`].
pub const TRAINING_HOURS: f64 = 25.0;

/// How much of a training run is discarded from the front, in rows.
///
/// Zero for `d00`, which is why `d00.dat` is 500 rows and the rest are 480.
pub const TRAINING_DISCARDED_ROWS: usize = D00_TRAINING_ROWS - TRAINING_ROWS;

/// Columns in a published file: `XMEAS(1..41)` then `XMV(1..11)`.
///
/// Fifty-two, not the fifty-three a [`crate::Run`] carries. `XMV(12)`, the
/// agitator speed, is in no published file. No controller in `temain_mod.f`
/// writes it, so it holds `TEINIT`'s `YY(50) = 50` for every run.
pub const COLUMNS: usize = MEASUREMENTS + MANIPULATED - 1;

/// Significant digits surviving `FORMAT(1X,E13.5)` at `temain_mod.f:1358`.
pub const PUBLISHED_DIGITS: usize = 5;

/// **The training protocol is a hypothesis, not a record.**
///
/// `temain_mod.f` generates 960 rows with the fault at row 160. Every training
/// file is 480 rows (500 for `d00`) with the fault present from row 0. Whoever
/// made them edited constants that were never committed.
///
/// One hypothesis explains both row counts at once, and it is the one
/// [`File::scenario`] builds: run 25 hours with the disturbance on from the
/// first step, record 500 rows, and drop the first hour from every file except
/// `d00`. Twenty-five hours at three minutes is 500 rows exactly, and 500 minus
/// the twenty rows of one hour is 480 exactly.
///
/// A competing hypothesis, 24 hours with the fault from step zero and nothing
/// discarded, fits the 480 just as well and leaves `d00`'s 500 unexplained. It
/// is rejected on that, and on nothing stronger.
pub const TRAINING_IS_A_HYPOTHESIS: bool = true;

/// Which half of the dataset a file belongs to.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Split {
    /// `dNN.dat`: fault from the first step, 480 rows, undocumented protocol.
    Training,
    /// `dNN_te.dat`: fault at hour eight, 960 rows, straight from the driver.
    Testing,
}

impl Split {
    /// The filename suffix before `.dat`.
    #[must_use]
    pub const fn suffix(self) -> &'static str {
        match self {
            Self::Training => "",
            Self::Testing => "_te",
        }
    }
}

/// Why a published file cannot be regenerated.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Unavailable {
    /// `d21` names a disturbance this revision of `teprob.f` does not have.
    ///
    /// `teprob.f:340` loops `DO 500 I=1,20` and the header at
    /// `teprob.f:172-191` lists twenty. Yet `d21.dat` and `d21_te.dat` ship and
    /// `teprob.f:1207,1233` carry seeds for both. Those files were made with a
    /// revision that is not the vendored one, so no scenario here reproduces
    /// them.
    FaultNotInThisRevision {
        /// The disturbance the file names.
        fault: usize,
    },
}

/// One published file.
#[derive(Clone, Copy, Debug)]
pub struct File {
    /// The disturbance number, `0` for fault-free.
    pub fault: usize,
    /// Which half.
    pub split: Split,
    /// Rows in the shipped file.
    pub rows: usize,
    /// The generator word `teprob.f:1187-1256` records for it.
    pub seed: f64,
}

impl File {
    /// The filename without its extension, such as `d01_te`.
    #[must_use]
    pub fn stem(&self) -> String {
        format!("d{:02}{}", self.fault, self.split.suffix())
    }

    /// Rows generated and then discarded from the front.
    ///
    /// Always zero for a testing file, and for `d00` training. See
    /// [`TRAINING_IS_A_HYPOTHESIS`].
    #[must_use]
    pub const fn discarded_rows(&self) -> usize {
        match self.split {
            Split::Testing => 0,
            Split::Training => {
                if self.fault == 0 {
                    0
                } else {
                    TRAINING_DISCARDED_ROWS
                }
            }
        }
    }

    /// The scenario that generates this file's geometry.
    ///
    /// `driver_forces_idv12` is **off**. The published bytes settle that
    /// question: every `_te` file other than `d12_te` sits at the nominal
    /// operating point across hour eight, which it could not do if the driver's
    /// `IDV(12)=1` at `temain_mod.f:367` had been left in. Delta D-011.
    ///
    /// # Errors
    ///
    /// [`Unavailable::FaultNotInThisRevision`] for `d21`.
    pub fn scenario(&self) -> Result<Scenario, Unavailable> {
        if self.fault > crate::DISTURBANCES {
            return Err(Unavailable::FaultNotInThisRevision { fault: self.fault });
        }
        let mut scenario = Scenario::baseline()
            .with_seed(self.seed)
            .sampling_every(SAMPLE_EVERY);
        scenario.driver_forces_idv12 = false;
        // Pinned rather than left to the default. `teprob.f:807-811` freezes a
        // tripped plant and keeps reporting, so every published file has its
        // full row count whether or not the plant went down. A build that ends
        // a run on a trip would produce a short file, which would not be the
        // geometry this module promises. Delta D-007.
        scenario.quirks.trip_ends_the_run = false;

        scenario = match self.split {
            Split::Testing => {
                let base = scenario.with_hours(TESTING_HOURS);
                if self.fault == 0 {
                    base
                } else {
                    base.with_event(Event::start(TESTING_FAULT_HOURS, self.fault))
                }
            }
            Split::Training => {
                let base = scenario.with_hours(TRAINING_HOURS);
                if self.fault == 0 {
                    base
                } else {
                    base.with_fault(self.fault)
                }
            }
        };
        Ok(scenario)
    }
}

/// Every published file, `d00` through `d21`, training then testing.
///
/// Forty-four entries. The seeds are `teprob.f:1187-1256` in order;
/// `tests/published.rs` reads them back out of the vendored source.
pub const FILES: [File; 44] = {
    // `d00_tr` through `d21_tr`, the commented block at `teprob.f:1190-1213`.
    //
    // `teprob.f:1188` carries a second candidate for `d00`, labelled
    // `d00_tr_new`, and `teprob.f:1209` a second for `d18`, labelled
    // `dd18_tr`. Which of each pair made the shipped file is stated nowhere,
    // so the plainly-labelled one is used and the alternative is
    // [`D00_TRAINING_ALTERNATIVE`] / [`D18_TRAINING_ALTERNATIVE`].
    const TR: [f64; 22] = [
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
    ];
    // `d00_te` through `d21_te`, `teprob.f:1220-1241`.
    const TE: [f64; 22] = [
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
    ];

    let mut files = [File {
        fault: 0,
        split: Split::Training,
        rows: D00_TRAINING_ROWS,
        seed: TR[0],
    }; 44];
    let mut fault = 0;
    while fault < 22 {
        files[fault] = File {
            fault,
            split: Split::Training,
            rows: if fault == 0 {
                D00_TRAINING_ROWS
            } else {
                TRAINING_ROWS
            },
            seed: TR[fault],
        };
        files[22 + fault] = File {
            fault,
            split: Split::Testing,
            rows: TESTING_ROWS,
            seed: TE[fault],
        };
        fault += 1;
    }
    files
};

/// The other seed offered for `d00_tr`, labelled `d00_tr_new` at
/// `teprob.f:1188`.
pub const D00_TRAINING_ALTERNATIVE: f64 = 5_687_912_315.0;

/// The other seed offered for `d18_tr`, labelled `dd18_tr` at `teprob.f:1209`.
pub const D18_TRAINING_ALTERNATIVE: f64 = 1_234_567_890.0;

/// The seed compiled into `teprob.f:1187`, which is what an unmodified build
/// runs with.
pub const COMPILED_IN_SEED: f64 = 4_651_207_995.0;

/// The seed labelled "original" at `teprob.f:1189`.
pub const ORIGINAL_SEED: f64 = 1_431_655_765.0;

/// `XMV(12)`, the agitator, as it stands in every published run.
///
/// Not a guess: no controller in `temain_mod.f` writes `XMV(12)`, so it holds
/// `TEINIT`'s `YY(50)` for the whole run.
pub const UNRECORDED_AGITATOR: f64 = 50.0;

/// One sample in published column order: `XMEAS(1..41)` then `XMV(1..11)`.
///
/// `XMV(12)` is dropped, because no published file has it.
#[must_use]
pub fn row(sample: &Sample) -> [f64; COLUMNS] {
    let mut out = [0.0; COLUMNS];
    out[..MEASUREMENTS].copy_from_slice(&sample.measurements);
    out[MEASUREMENTS..].copy_from_slice(&sample.manipulated[..MANIPULATED - 1]);
    out
}

/// Column headings in published order, for the CSV form.
#[must_use]
pub fn column_names() -> [&'static str; COLUMNS] {
    let all = crate::run::channel_names();
    let mut names = [""; COLUMNS];
    names[..MEASUREMENTS].copy_from_slice(&all[..MEASUREMENTS]);
    names[MEASUREMENTS..].copy_from_slice(&all[MEASUREMENTS..MEASUREMENTS + MANIPULATED - 1]);
    names
}

/// Round as `FORMAT(1X,E13.5)` did, to [`PUBLISHED_DIGITS`] significant digits.
///
/// # Why generated data is degraded on purpose
///
/// The shipped files went through this and cannot be un-rounded. A generated
/// file that carried seventeen digits would not be the same kind of object, and
/// anything fitted to it would see structure in digits the published data never
/// had. `tep dataset --format csv` is the way to ask for full precision, and it
/// says in its header that it is not published-shaped.
#[must_use]
pub fn round_as_published(value: f64) -> f64 {
    if !value.is_finite() {
        return value;
    }
    // Through the decimal text rather than by scaling, because `v * 10f64.powi(k)`
    // rounds twice and disagrees with what a formatter does in the last digit.
    let text = format!("{:.*e}", PUBLISHED_DIGITS - 1, value);
    text.parse().unwrap_or(value)
}

/// The Rieth ensemble: 500 runs per fault, and a different training protocol.
///
/// `rieth-2017-addit`, the Harvard Dataverse set most current machine-learning
/// work on this process is trained against. It is not the original `d00`-`d21`
/// distribution and it is not generated the same way, which is the reason this
/// module exists separately from the constants above.
///
/// # Where it differs from the original, and why that matters
///
/// The testing protocol is the same: 48 hours, 960 rows, the fault at hour
/// eight. The *training* protocol is not, and the difference is easy to miss
/// because both are 25 hours and both split twenty rows against four hundred
/// and eighty.
///
/// [`TRAINING_HOURS`] here runs 25 hours with the fault arriving at hour one
/// and keeps all 500 rows, so the first 20 are fault-free and the remaining 480
/// are faulted. [`File::scenario`] instead runs 25 hours with the fault
/// live from the first step and *discards* the first hour, which is the only
/// hypothesis that explains `d00.dat`'s 500 rows against every other training
/// file's 480. Those are different runs, and one file cannot stand in for the
/// other.
///
/// # This generates their shape, not their data
///
/// The same caveat the rest of this module carries. Their seeds are not
/// recorded anywhere this repository can read, the toolchain that produced
/// their files is not recorded either, and `IDV(21)` is not in this revision of
/// the model.
pub mod rieth {
    use super::{SAMPLE_EVERY, Split, Unavailable};
    use crate::scenario::Scenario;
    use tepsim_scenario::Event;

    /// Simulation runs per fault, per split, in the published distribution.
    pub const RUNS: usize = 500;

    /// Hours in a training run.
    pub const TRAINING_HOURS: f64 = 25.0;

    /// Rows in a training run: the whole 25 hours at three minutes.
    pub const TRAINING_ROWS: usize = 500;

    /// When the fault arrives in a faulted training run.
    pub const TRAINING_ONSET_HOURS: f64 = 1.0;

    /// Fault-free rows at the head of a faulted training run.
    pub const TRAINING_NORMAL_ROWS: usize = 20;

    /// Hours in a testing run.
    pub const TESTING_HOURS: f64 = 48.0;

    /// Rows in a testing run.
    pub const TESTING_ROWS: usize = 960;

    /// When the fault arrives in a faulted testing run.
    pub const TESTING_ONSET_HOURS: f64 = 8.0;

    /// Fault-free rows at the head of a faulted testing run.
    pub const TESTING_NORMAL_ROWS: usize = 160;

    /// The generator word for one run.
    ///
    /// # Why this is not called independent
    ///
    /// The dataset's own description says its runs use independent,
    /// non-overlapping seeds. This function does not make that claim, because
    /// with this generator it is not one that can be checked cheaply.
    /// `teprob.f`'s `TESUB7` is `G = mod(G * 9228907, 2^32)` evaluated in double
    /// precision, and the product exceeds 2^53 for most states, so the modulus
    /// is inexact and the sequence is not the clean multiplicative congruential
    /// generator it looks like. A cycle search from the compiled-in seed found
    /// no return to its starting state within 2e9 states, which rules out the
    /// group-theoretic period of such a generator and settles nothing about
    /// whether two streams overlap.
    ///
    /// So what is promised here is only what can be delivered: the seeds are
    /// distinct, spread across the generator's range, and a deterministic
    /// function of `(fault, split, run)`, so a run is reproducible from its
    /// coordinates alone.
    #[must_use]
    pub fn seed(fault: usize, split: Split, run: usize) -> f64 {
        // Odd, distinct per coordinate, and comfortably inside the range the
        // original's own recorded seeds occupy (about 1e9 to 9.9e9). The strides
        // are coprime with each other so no two coordinates collide.
        let split_offset = match split {
            Split::Training => 0_u64,
            Split::Testing => 1,
        };
        let index = (fault as u64) * 2 * RUNS as u64 + split_offset * RUNS as u64 + run as u64;
        let word = 1_000_000_007_u64 + index * 7_919_u64 * 2;
        word as f64
    }

    /// The scenario for one run of the ensemble.
    ///
    /// # Errors
    ///
    /// [`Unavailable::FaultNotInThisRevision`] for `IDV(21)`, which this
    /// revision of `teprob.f` does not contain.
    pub fn scenario(fault: usize, split: Split, run: usize) -> Result<Scenario, Unavailable> {
        if fault > crate::DISTURBANCES {
            return Err(Unavailable::FaultNotInThisRevision { fault });
        }
        let (hours, onset) = match split {
            Split::Training => (TRAINING_HOURS, TRAINING_ONSET_HOURS),
            Split::Testing => (TESTING_HOURS, TESTING_ONSET_HOURS),
        };
        let mut scenario = Scenario::baseline()
            .with_seed(seed(fault, split, run))
            .with_hours(hours)
            .sampling_every(SAMPLE_EVERY);
        // The same two pins `File::scenario` carries and for the same reasons:
        // the published bytes say IDV(12) was not forced, and a published file
        // has its full row count whether or not the plant tripped.
        scenario.driver_forces_idv12 = false;
        scenario.quirks.trip_ends_the_run = false;
        if fault != 0 {
            scenario = scenario.with_event(Event::start(onset, fault));
        }
        Ok(scenario)
    }

    /// Rows in a run of this split.
    #[must_use]
    pub const fn rows(split: Split) -> usize {
        match split {
            Split::Training => TRAINING_ROWS,
            Split::Testing => TESTING_ROWS,
        }
    }

    /// Fault-free rows at the head of a faulted run of this split.
    #[must_use]
    pub const fn normal_rows(split: Split) -> usize {
        match split {
            Split::Training => TRAINING_NORMAL_ROWS,
            Split::Testing => TESTING_NORMAL_ROWS,
        }
    }
}
