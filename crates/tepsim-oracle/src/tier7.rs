//! Tier 7: the published `d00`-`d21` datasets.
//!
//! B-0051. Tier 5 asks whether the port and the Fortran, run under *our*
//! protocol, are statistically the same simulator. Tier 7 asks a harder and
//! less forgiving question: can the port regenerate the forty-four files that
//! thirty years of fault-detection papers were written against?
//!
//! # This module is a measurement, not a gate
//!
//! `PLAN.org` is explicit about the discipline, and it is the whole point of
//! the item:
//!
//! > Where reproduction is imperfect, document exactly which protocol detail
//! > is unknown rather than tuning until the numbers match.
//!
//! So nothing here is fitted. Every number in [`Protocol`] is either read off
//! `temain_mod.f` or derived from a row count, and which of the two it is is
//! recorded in [`Protocol::source`]. A protocol detail that cannot be pinned
//! down is an [`Unknown`], not a free parameter.
//!
//! # What the source actually states
//!
//! `temain_mod.f` is the driver that made these files, and it states four
//! things outright:
//!
//! | Fact | Line | Value |
//! |---|---|---|
//! | Run length | `temain_mod.f:220` | `NPTS = 172800` steps, 48 h at one second |
//! | Settling time | `temain_mod.f:226` | `SSPTS = 3600 * 8`, eight hours |
//! | Output cadence | `temain_mod.f:401` | every 180 steps, so 960 rows |
//! | Output format | `temain_mod.f:1358` | `FORMAT(1X,E13.5)`, five significant digits |
//!
//! and one thing in prose, at `temain_mod.f:101-102`: *"Go to line 367,
//! implement any of the 21 programmed disturbances."* Line 367 is inside
//! `IF (I.GE.SSPTS)`, so the fault arrives at eight hours rather than at the
//! start. The data agrees: every `dNN_te.dat` sits at the nominal operating
//! point for its first 160 rows and departs from it in row 160.
//!
//! That accounts for all twenty-two `_te` files exactly. It accounts for none
//! of the twenty-two training files, which are 480 rows and 500 rows against
//! the driver's 960, and the source says nothing whatever about how they were
//! made. See [`Protocol::source`] and [`Unknown::TrainingProtocolIsUndocumented`].
//!
//! # Structure
//!
//! [`FILES`] is the inventory. [`Published::load`] reads one, transposing
//! `d00.dat` because it alone is stored the other way up. [`generate`] runs
//! the port under a [`Protocol`]. Both come back as a [`tier5::Run`] so that
//! [`tier5::battery`] can judge them without a second implementation of
//! anything.

use std::path::{Path, PathBuf};

use crate::tier5::{self, Run, Scenario, VARIABLES};

/// Columns in a published file: `XMEAS(1..41)` then `XMV(1..11)`.
///
/// Fifty-two, not the fifty-three a [`Run`] carries. `XMV(12)`, the agitator
/// speed, is not recorded in any published file. See [`UNRECORDED_AGITATOR`].
pub const PUBLISHED_COLUMNS: usize = 52;

/// `XMV(12)` is absent from every published file, so the published side of a
/// comparison has to supply it from somewhere.
///
/// It is supplied as a constant, and the constant is not a guess: no
/// controller in `temain_mod.f` writes `XMV(12)`, so it holds `TEINIT`'s
/// `YY(50) = 50` for the whole run. `the_agitator_is_constant_in_the_port`
/// asserts the port does the same rather than taking it on faith, and
/// [`tier5::battery::VariableReport::constant`] then skips the moment gates
/// for it on both sides.
pub const UNRECORDED_AGITATOR: f64 = 50.0;

/// Rows in a `dNN_te.dat` file. 48 h at one row per 180 s.
pub const TESTING_ROWS: usize = 960;

/// Rows in a `dNN.dat` training file, for `NN >= 1`.
pub const TRAINING_ROWS: usize = 480;

/// Rows in `d00.dat`, which is the only training file that is not 480.
///
/// Upstream's own `reference/data/README.md` says 480. It is wrong, and the
/// twenty extra rows are the single most useful piece of evidence about the
/// undocumented training protocol; see [`Protocol::training`].
pub const D00_TRAINING_ROWS: usize = 500;

/// `NPTS`, `temain_mod.f:220`: the run length the driver ships with, in
/// one-second steps.
pub const DRIVER_STEPS: usize = 172_800;

/// `SSPTS`, `temain_mod.f:226`: how long the driver runs before the fault.
pub const DRIVER_SETTLING_STEPS: usize = 3600 * 8;

/// Significant digits in a published value. `FORMAT(1X,E13.5)`,
/// `temain_mod.f:1358-1360`.
///
/// This is a real limit on what Tier 7 can conclude. A published row is the
/// simulator's output rounded to five figures, so two runs agreeing to five
/// figures are indistinguishable in these files however far apart they are in
/// the sixth.
pub const SIGNIFICANT_DIGITS: u32 = 5;

/// Which half of the published pair a file belongs to.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Split {
    /// `dNN.dat`, the training half.
    Training,
    /// `dNN_te.dat`, the testing half.
    Testing,
}

impl Split {
    /// The suffix this split's filenames carry.
    #[must_use]
    pub const fn suffix(self) -> &'static str {
        match self {
            Self::Training => "",
            Self::Testing => "_te",
        }
    }
}

/// Where a number in a [`Protocol`] came from.
///
/// The distinction is the deliverable. A protocol quantity the driver states
/// is a fact about how the files were made; one inferred from a row count is a
/// hypothesis that happens to fit, and a reader has to be able to tell them
/// apart without reading this module's history.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Provenance {
    /// Read off `temain_mod.f`.
    Stated,
    /// Derived from the shape of the published files, because the source is
    /// silent.
    Inferred,
}

/// How one published file is claimed to have been generated.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Protocol {
    /// Simulated hours.
    pub hours: f64,
    /// When the file's own fault is switched on, in hours. `None` for `d00`,
    /// which has none.
    pub onset_hours: Option<f64>,
    /// Leading samples produced but not published.
    pub discarded: usize,
    /// Rows the file should then hold.
    pub rows: usize,
    /// Whether the numbers above are stated by the source or inferred from the
    /// data.
    pub source: Provenance,
}

impl Protocol {
    /// The testing protocol, entirely from `temain_mod.f`.
    ///
    /// 48 hours, fault at eight, nothing discarded, 960 rows. Every one of
    /// those four numbers is in the driver, and the row count of all
    /// twenty-two `_te` files matches.
    #[must_use]
    pub const fn testing(has_fault: bool) -> Self {
        Self {
            hours: 48.0,
            onset_hours: if has_fault { Some(8.0) } else { None },
            discarded: 0,
            rows: TESTING_ROWS,
            source: Provenance::Stated,
        }
    }

    /// The training protocol, which the source does not state.
    ///
    /// Reconstructed from two observations and one piece of outside
    /// knowledge:
    ///
    /// - `d00.dat` holds 500 rows, which is 25 hours at 180 seconds. Nothing
    ///   in the driver produces 500 of anything; `NPTS` would have to have
    ///   been edited, which `temain_mod.f:95-96` tells the user to do.
    /// - Every other training file holds 480, which is 500 minus 20, and 20
    ///   rows is exactly one hour.
    /// - Russell, Chiang and Braatz describe the training sets as 25-hour runs
    ///   with the fault introduced after one hour, the first hour discarded.
    ///
    /// So: 25 hours, fault at one hour, first 20 rows dropped, 480 published.
    /// `d00` has no fault, so nothing is dropped and all 500 rows ship.
    ///
    /// **This is a hypothesis that fits two row counts.** It is not stated
    /// anywhere in the vendored source, and it is not the only arithmetic that
    /// produces 480 rows: a 24-hour run with the fault on from the first step
    /// produces 480 rows too, and puts the same three minutes between fault
    /// onset and the first published row. See
    /// [`Unknown::TrainingProtocolIsUndocumented`].
    #[must_use]
    pub const fn training(has_fault: bool) -> Self {
        Self {
            hours: 25.0,
            onset_hours: if has_fault { Some(1.0) } else { None },
            discarded: if has_fault { 20 } else { 0 },
            rows: if has_fault {
                TRAINING_ROWS
            } else {
                D00_TRAINING_ROWS
            },
            source: Provenance::Inferred,
        }
    }

    /// Samples the run produces before anything is discarded.
    #[must_use]
    pub fn samples(&self) -> usize {
        self.rows + self.discarded
    }
}

/// One published file.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Published {
    /// `IDV` index, or zero for `d00`.
    pub fault: usize,
    /// Which half.
    pub split: Split,
}

impl Published {
    /// The file's stem, as papers cite it: `d00`, `d13_te`.
    #[must_use]
    pub fn name(&self) -> String {
        format!("d{:02}{}", self.fault, self.split.suffix())
    }

    /// Its path under `reference/data/`.
    #[must_use]
    pub fn path(&self) -> PathBuf {
        data_dir().join(format!("{}.dat", self.name()))
    }

    /// The scenario this file claims to be.
    ///
    /// `IDV(21)` has no entry in this model, so [`Published::is_representable`]
    /// gates the call; see [`Unknown::Idv21DoesNotExistInThisRevision`].
    ///
    /// # Panics
    ///
    /// If the fault is outside `1..=tier5::FAULTS`.
    #[must_use]
    pub fn scenario(&self) -> Scenario {
        if self.fault == 0 {
            Scenario::NOMINAL
        } else {
            Scenario::fault(self.fault)
        }
    }

    /// Whether this model can express the file's fault at all.
    ///
    /// False for `d21` alone. The vendored `teprob.f:340` loops `DO 500 I=1,20`
    /// and `IDV(21)` is not in it, so there is no seed, no run and no
    /// comparison to make. The file still ships and is still inventoried.
    #[must_use]
    pub const fn is_representable(&self) -> bool {
        self.fault <= tier5::FAULTS
    }

    /// The generator word `teprob.f` records for this file.
    ///
    /// `teprob.f:1187-1256`, transcribed in [`tier5::published_seeds`]. There
    /// is one for every file including `d21`, so this is defined even where
    /// [`Published::is_representable`] is false.
    #[must_use]
    pub fn seed(&self) -> f64 {
        match self.split {
            Split::Training => tier5::published_seeds::TRAINING[self.fault],
            Split::Testing => tier5::published_seeds::TESTING[self.fault],
        }
    }

    /// The protocol this file is claimed to have been generated under.
    #[must_use]
    pub fn protocol(&self) -> Protocol {
        match self.split {
            Split::Training => Protocol::training(self.fault != 0),
            Split::Testing => Protocol::testing(self.fault != 0),
        }
    }

    /// Read it, as rows of 52.
    ///
    /// `d00.dat` is stored transposed, 52 rows of 500 columns where every
    /// other file is samples-by-52, and is transposed back here. That is not a
    /// guess: 52 is the column count of all forty-three other files and the
    /// width of the observation vector, and `reference/README.org` records the
    /// same finding independently.
    ///
    /// # Panics
    ///
    /// If the file is missing, unparsable, or not rectangular.
    #[must_use]
    pub fn load(&self) -> Vec<[f64; PUBLISHED_COLUMNS]> {
        let path = self.path();
        let text = std::fs::read_to_string(&path)
            .unwrap_or_else(|e| panic!("reading {}: {e}", path.display()));
        let grid: Vec<Vec<f64>> = text
            .lines()
            .filter(|line| !line.trim().is_empty())
            .map(|line| {
                line.split_whitespace()
                    .map(|token| {
                        token
                            .parse::<f64>()
                            .unwrap_or_else(|e| panic!("{}: {token:?}: {e}", path.display()))
                    })
                    .collect()
            })
            .collect();

        let width = grid.first().map_or(0, Vec::len);
        assert!(
            grid.iter().all(|row| row.len() == width),
            "{} is ragged",
            path.display()
        );

        // The transpose is decided by the shape, not by the filename, so a
        // second transposed file would be handled and a `d00.dat` that was
        // ever regenerated the usual way up would not silently break.
        let rows: Vec<Vec<f64>> = if width == PUBLISHED_COLUMNS {
            grid
        } else {
            assert_eq!(
                grid.len(),
                PUBLISHED_COLUMNS,
                "{}: neither {PUBLISHED_COLUMNS} columns nor {PUBLISHED_COLUMNS} rows",
                path.display()
            );
            (0..width)
                .map(|c| grid.iter().map(|row| row[c]).collect())
                .collect()
        };

        rows.into_iter()
            .map(|row| {
                let mut out = [0.0; PUBLISHED_COLUMNS];
                out.copy_from_slice(&row);
                out
            })
            .collect()
    }

    /// Read it as a [`Run`], so the Tier 5 battery can consume it.
    ///
    /// The fifty-third channel is [`UNRECORDED_AGITATOR`]; see that constant
    /// for why supplying it is not fabricating data.
    #[must_use]
    pub fn run(&self) -> Run {
        let samples = self
            .load()
            .into_iter()
            .map(|row| {
                let mut out = [0.0; VARIABLES];
                out[..PUBLISHED_COLUMNS].copy_from_slice(&row);
                out[PUBLISHED_COLUMNS] = UNRECORDED_AGITATOR;
                out
            })
            .collect();
        Run {
            scenario: if self.is_representable() {
                self.scenario()
            } else {
                Scenario::NOMINAL
            },
            seed: self.seed(),
            samples,
            tripped: None,
        }
    }
}

/// Every published file, `d00` through `d21`, training then testing.
#[must_use]
pub fn files() -> Vec<Published> {
    let mut out = Vec::with_capacity(44);
    for split in [Split::Training, Split::Testing] {
        for fault in 0..=21 {
            out.push(Published { fault, split });
        }
    }
    out
}

/// How many published files there are.
pub const FILES: usize = 44;

/// Run the port under a protocol, and hand back what the file would have held.
///
/// The fault is switched on at [`Protocol::onset_hours`] through
/// [`tepsim::Simulation::request_disturbance`], not at step zero, and the
/// leading [`Protocol::discarded`] samples are dropped afterwards, so the
/// returned run is row-for-row what a published file claims to be.
///
/// `Scenario::driver_forces_idv12` is left at its default, which is on. That
/// is not a choice made here: it is `temain_mod.f:366-368`, the line the
/// driver ships with, and delta D-011. The other reading of that line is
/// [`generate_without_forced_idv12`], which is a diagnostic and not a knob.
#[must_use]
pub fn generate(file: &Published, seed: f64, protocol: &Protocol) -> Run {
    run_under(file, seed, protocol, true)
}

/// The same, with `temain_mod.f:367` read as an *example to replace* rather
/// than a line to keep.
///
/// This is the second half of [`Unknown::Idv12MayOrMayNotBeForced`], and it
/// exists so that the two readings can be measured against the published data
/// instead of argued about. It is a diagnostic and never the default:
/// [`generate`] keeps the line, because the line is what the driver ships
/// with.
#[must_use]
pub fn generate_without_forced_idv12(file: &Published, seed: f64, protocol: &Protocol) -> Run {
    run_under(file, seed, protocol, false)
}

fn run_under(file: &Published, seed: f64, protocol: &Protocol, forced_idv12: bool) -> Run {
    // `faithful`: four of the forty-four published files carry a frozen tail
    // from `teprob.f:807-811`, 1,832 rows in total and 75.6% of `d06.dat`, so a
    // build that ended the run at the trip would have nothing to compare them
    // against. D-011 is then set from the argument, because which way it went is
    // the question this function exists to answer.
    let mut described = tepsim::Scenario::faithful()
        .with_hours(protocol.hours)
        .with_seed(seed);
    described.driver_forces_idv12 = forced_idv12;
    // Floored at one so that an onset of zero hours means "live for the first
    // step", which is what the competing training hypothesis needs; step zero
    // does not exist, `temain_mod.f`'s `I` being one-based.
    let onset_step = protocol
        .onset_hours
        .map(|h| ((h * 3600.0).round() as usize).max(1))
        .unwrap_or(usize::MAX);

    let mut sim = tepsim::Simulation::new(described);
    let mut samples = Vec::with_capacity(protocol.samples());
    let mut tripped = None;
    while !sim.is_halted() {
        // One-based, as `temain_mod.f`'s `I` is, and switched *before* the
        // step so that the comparison `I .GE. SSPTS` matches the driver's.
        if file.fault != 0 && sim.steps_taken() + 1 == onset_step {
            sim.request_disturbance(file.fault, true);
        }
        // `step` returns `None` on the 179 steps in 180 that are not sample
        // steps as well as at the end of the run, so the loop is driven by
        // `is_halted` and not by the option.
        if let Some(sample) = sim.step() {
            samples.push(sample.row());
        }
    }
    if let Some(tepsim::Outcome::Tripped { step, .. } | tepsim::Outcome::SolveFailed { step }) =
        sim.outcome()
    {
        tripped = Some(step);
    }

    assert_eq!(
        samples.len(),
        protocol.samples(),
        "{} produced the wrong number of samples",
        file.name()
    );
    samples.drain(..protocol.discarded);

    Run {
        scenario: if file.is_representable() {
            file.scenario()
        } else {
            Scenario::NOMINAL
        },
        seed,
        samples,
        tripped,
    }
}

/// A protocol detail Tier 7 could not pin down.
///
/// This enumeration *is* the deliverable of B-0051. Each variant names one
/// thing that is not known about how the published files were generated, and
/// says what evidence would settle it. None of them is a knob: nothing in this
/// module was adjusted to make a number improve.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub enum Unknown {
    /// Nothing in the vendored source describes the training runs.
    ///
    /// `temain_mod.f` ships with `NPTS = 172800` and `SSPTS = 28800`, which
    /// produce 960 rows with the fault at row 160. The training files are 480
    /// rows (500 for `d00`) with the fault at row 0. Whoever made them edited
    /// two constants and the edit is not recorded. 25 hours with the first
    /// hour dropped and 24 hours with the fault from the first step both fit
    /// the 480, and only the first also explains `d00`'s 500.
    TrainingProtocolIsUndocumented,

    /// Whether `IDV(12)` was left forced on when another fault was requested.
    ///
    /// `temain_mod.f:367` is literally `IDV(12)=1`, and `temain_mod.f:101-102`
    /// tells the user to "go to line 367" and type the disturbance they want.
    /// Replacing that line gives a run with one fault; adding a line beside it
    /// gives a run with two, the requested one and `IDV(12)`. The instruction
    /// reads as replacement and the shipped line reads as an example, and the
    /// source decides between them nowhere. It matters for all forty-three
    /// files other than `d12`.
    ///
    /// **The published bytes decide it**, which is why this variant carries a
    /// [`Unknown::settled_by_evidence`] and the others mostly do not. See
    /// `the_published_files_were_not_generated_with_the_forced_idv12`.
    Idv12MayOrMayNotBeForced,

    /// `d21` names a fault this revision of `teprob.f` does not have.
    ///
    /// `teprob.f:340` loops `DO 500 I=1,20`, the header at
    /// `teprob.f:172-191` lists twenty, and `IDV(21)` appears nowhere in the
    /// model. Yet `d21.dat` and `d21_te.dat` ship, `reference/data/README.md`
    /// calls the fault "valve position constant (Stream 4)", and
    /// `teprob.f:1207,1233` carry seeds for both. The files were made with a
    /// revision that is not the vendored one.
    Idv21DoesNotExistInThisRevision,

    /// The five-figure output format destroys the evidence bit-exactness
    /// would need.
    ///
    /// `FORMAT(1X,E13.5)` at `temain_mod.f:1358` keeps five significant
    /// digits. Tier 4 shows trajectories separating on `exp` and `pow`
    /// rounding within hours, so even the correct protocol on the correct
    /// revision would not reproduce a 48-hour trajectory to five figures. Tier
    /// 7 is therefore a distributional question by construction, not an
    /// exactness one.
    OutputIsRoundedToFiveFigures,

    /// Neither the compiler nor the machine that generated the files is
    /// recorded.
    ///
    /// The files predate the vendored repository by two decades. Their `exp`
    /// and `pow` are some 1990s Fortran runtime's, and the port's are the
    /// vendored `libm`'s. Tier 2 measures that difference at one ULP on about
    /// ten percent of arguments, and Tier 4 measures what one ULP does to a
    /// 48-hour trajectory.
    GeneratingToolchainIsUnrecorded,

    /// Two seeds are recorded for `d00_tr` and two for `d18_tr`.
    ///
    /// `teprob.f:1188` labels one `d00_tr_new` and `teprob.f:1209` labels one
    /// `dd18_tr`. Which of each pair produced the shipped file is not stated,
    /// so `d00` and `d18` each have two candidate seeds.
    TwoSeedsAreRecordedForSomeFiles,

    /// `d19_te`'s seed does not compile.
    ///
    /// `teprob.f:1239` reads `G=9090909232.DO`, ending in the letter O. It is
    /// inside a comment so nothing ever caught it, and the intended digit is
    /// unrecoverable from the source: `9090909232` is the reading assumed
    /// here, and nothing confirms it.
    D19TestingSeedIsMalformed,

    /// The setpoints the runs were made at are not recorded either.
    ///
    /// `temain_mod.f:245-330` sets nineteen setpoints and twenty controller
    /// gains, and `temain_mod.f:334` shows a commented-out setpoint move. A
    /// dataset generated after any such edit is not reproducible from the
    /// shipped constants, and nothing in a published file says which
    /// constants it was made with.
    ControlSchemeConstantsAreNotRecordedPerFile,
}

impl Unknown {
    /// Every unknown, for a report to enumerate.
    #[must_use]
    pub const fn all() -> [Self; 8] {
        [
            Self::TrainingProtocolIsUndocumented,
            Self::Idv12MayOrMayNotBeForced,
            Self::Idv21DoesNotExistInThisRevision,
            Self::OutputIsRoundedToFiveFigures,
            Self::GeneratingToolchainIsUnrecorded,
            Self::TwoSeedsAreRecordedForSomeFiles,
            Self::D19TestingSeedIsMalformed,
            Self::ControlSchemeConstantsAreNotRecordedPerFile,
        ]
    }

    /// What the published data settles, where the source cannot.
    ///
    /// `None` means the question is still open: the source is silent and the
    /// bytes do not decide it either. Those are the ones that stop
    /// reproduction, and they are the deliverable of B-0051.
    ///
    /// A `Some` does not delete the unknown. The source still fails to state
    /// the answer, so a reader of `temain_mod.f` alone would still not know
    /// it, and the next dataset generated from that driver would still get it
    /// wrong. What it records is that Tier 7 measured the answer rather than
    /// guessed it.
    #[must_use]
    pub const fn settled_by_evidence(self) -> Option<&'static str> {
        match self {
            Self::Idv12MayOrMayNotBeForced => Some(
                "REPLACED. d12_te's spread jumps at row 160 by at least 5.3x what \
                 d00_te's does on every channel IDV(12) drives, and keeping the line \
                 inflates the port's spread against d00 by up to 9.9x against 1.55x \
                 for replacing it.",
            ),
            Self::TrainingProtocolIsUndocumented
            | Self::Idv21DoesNotExistInThisRevision
            | Self::OutputIsRoundedToFiveFigures
            | Self::GeneratingToolchainIsUnrecorded
            | Self::TwoSeedsAreRecordedForSomeFiles
            | Self::D19TestingSeedIsMalformed
            | Self::ControlSchemeConstantsAreNotRecordedPerFile => None,
        }
    }

    /// Which files it affects.
    #[must_use]
    pub const fn scope(self) -> &'static str {
        match self {
            Self::TrainingProtocolIsUndocumented => "all 22 training files",
            Self::Idv12MayOrMayNotBeForced => "all but d12 and d12_te",
            Self::Idv21DoesNotExistInThisRevision => "d21, d21_te",
            Self::OutputIsRoundedToFiveFigures
            | Self::GeneratingToolchainIsUnrecorded
            | Self::ControlSchemeConstantsAreNotRecordedPerFile => "all 44",
            Self::TwoSeedsAreRecordedForSomeFiles => "d00, d18",
            Self::D19TestingSeedIsMalformed => "d19_te",
        }
    }
}

/// `reference/data/`, from this crate's manifest.
fn data_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../reference/data")
        .canonicalize()
        .expect("reference/data exists")
}
