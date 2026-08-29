//! What a run produces: samples, labels, and how it ended.

use alloc::vec::Vec;

use tepsim_core::ShutdownCause;

use crate::scenario::{DISTURBANCES, Scenario};

/// Measurements per sample: `XMEAS(1..41)`.
pub const MEASUREMENTS: usize = 41;
/// Manipulated variables per sample: `XMV(1..12)`.
pub const MANIPULATED: usize = 12;
/// Channels recorded per sample.
pub const CHANNELS: usize = MEASUREMENTS + MANIPULATED;

/// One recorded instant.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Sample {
    /// Integrator step, one-based as `temain_mod.f`'s `I` is.
    pub step: usize,
    /// Simulated time, in hours.
    pub hours: f64,
    /// `XMEAS(1..41)`, with noise and with the analysers' dead time.
    pub measurements: [f64; MEASUREMENTS],
    /// `XMV(1..12)`.
    pub manipulated: [f64; MANIPULATED],
    /// Ground truth: what was actually wrong with the plant at this instant.
    pub labels: Labels,
}

impl Sample {
    /// The 53 channels in one row, measurements first.
    ///
    /// The layout every downstream consumer uses: the correlation matrix, the
    /// detectors in Tier 6, and the published files.
    #[must_use]
    pub fn row(&self) -> [f64; CHANNELS] {
        let mut out = [0.0; CHANNELS];
        out[..MEASUREMENTS].copy_from_slice(&self.measurements);
        out[MEASUREMENTS..].copy_from_slice(&self.manipulated);
        out
    }
}

/// Ground truth for one sample.
///
/// The original records nothing of the sort: a published dataset is a matrix of
/// numbers and a filename. Detection-delay figures in the literature are
/// therefore computed against whatever onset time the author assumed. Recording
/// it removes the guesswork.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Labels {
    /// Which disturbances were active, one-based: index 0 is `IDV(1)`.
    pub active: [bool; DISTURBANCES],
    /// Hours since each disturbance came on, or `None` if it never did.
    ///
    /// Not simply "time since the run began": the driver switches `IDV(12)` on
    /// at hour eight whatever the scenario asked for (delta D-011), so one of
    /// these onsets can be later than the others and is not the caller's doing.
    pub since_onset: [Option<f64>; DISTURBANCES],
}

impl Labels {
    /// Nothing wrong yet.
    #[must_use]
    pub const fn none() -> Self {
        Self {
            active: [false; DISTURBANCES],
            since_onset: [None; DISTURBANCES],
        }
    }

    /// Whether any disturbance was active.
    #[must_use]
    pub fn faulted(&self) -> bool {
        self.active.iter().any(|on| *on)
    }

    /// The active faults, one-based.
    pub fn faults(&self) -> impl Iterator<Item = usize> + '_ {
        self.active
            .iter()
            .enumerate()
            .filter(|(_, on)| **on)
            .map(|(index, _)| index + 1)
    }
}

/// How a run ended.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Outcome {
    /// It ran to the end of its scenario.
    Completed,
    /// The plant tripped.
    ///
    /// The run *stops* here by default. `teprob.f:807-811` instead freezes the
    /// plant and keeps reporting, which is what [`crate::Scenario::faithful`]
    /// and a cleared [`tepsim_core::QuirkFixes::trip_ends_the_run`] give, and
    /// what any comparison against published data needs. Delta D-007, signed
    /// off 2026-08-28.
    Tripped {
        /// The step it tripped at.
        step: usize,
        /// Simulated hours at that step.
        hours: f64,
        /// The first condition that fired.
        cause: Option<ShutdownCause>,
    },
    /// A temperature solve failed to converge.
    ///
    /// The original cannot report this: `TESUB2` returns its guess and claims
    /// success (delta D-001). The port declines to invent an answer.
    SolveFailed {
        /// The step it failed at.
        step: usize,
    },
}

impl Outcome {
    /// Whether the plant stayed up for the whole scenario.
    #[must_use]
    pub const fn is_completed(&self) -> bool {
        matches!(self, Self::Completed)
    }
}

/// A finished run: its scenario, its samples, and how it ended.
#[derive(Clone, Debug, PartialEq)]
pub struct Run {
    /// What was asked for.
    pub scenario: Scenario,
    /// The samples, in order.
    pub samples: Vec<Sample>,
    /// How it ended.
    pub outcome: Outcome,
}

impl Run {
    /// One channel across the whole run, as a column.
    ///
    /// Zero-based over the 53 channels: 0 to 40 are `XMEAS(1..41)` and 41 to
    /// 52 are `XMV(1..12)`.
    ///
    /// # Panics
    ///
    /// If `channel` is not below [`CHANNELS`].
    #[must_use]
    pub fn column(&self, channel: usize) -> Vec<f64> {
        assert!(channel < CHANNELS, "channel {channel} is out of range");
        self.samples.iter().map(|s| s.row()[channel]).collect()
    }

    /// All 53 columns.
    #[must_use]
    pub fn columns(&self) -> Vec<Vec<f64>> {
        (0..CHANNELS).map(|c| self.column(c)).collect()
    }

    /// One measurement across the run, one-based as `XMEAS(n)` is.
    ///
    /// # Panics
    ///
    /// If `n` is not in `1..=41`.
    #[must_use]
    pub fn measurement(&self, n: usize) -> Vec<f64> {
        assert!((1..=MEASUREMENTS).contains(&n), "XMEAS index out of range");
        self.samples.iter().map(|s| s.measurements[n - 1]).collect()
    }

    /// One manipulated variable across the run, one-based as `XMV(n)` is.
    ///
    /// # Panics
    ///
    /// If `n` is not in `1..=12`.
    #[must_use]
    pub fn manipulated(&self, n: usize) -> Vec<f64> {
        assert!((1..=MANIPULATED).contains(&n), "XMV index out of range");
        self.samples.iter().map(|s| s.manipulated[n - 1]).collect()
    }

    /// The step at which the plant tripped, if it did.
    #[must_use]
    pub const fn tripped_at(&self) -> Option<usize> {
        match self.outcome {
            Outcome::Tripped { step, .. } => Some(step),
            _ => None,
        }
    }
}

/// Names for the 53 channels, in row order.
///
/// Short and stable, for a CSV header or a column label. The measurement names
/// follow Downs and Vogel's Table 4 and the manipulated ones their Table 3.
#[must_use]
pub fn channel_names() -> [&'static str; CHANNELS] {
    [
        "XMEAS_1_A_feed",
        "XMEAS_2_D_feed",
        "XMEAS_3_E_feed",
        "XMEAS_4_total_feed",
        "XMEAS_5_recycle_flow",
        "XMEAS_6_reactor_feed_rate",
        "XMEAS_7_reactor_pressure",
        "XMEAS_8_reactor_level",
        "XMEAS_9_reactor_temperature",
        "XMEAS_10_purge_rate",
        "XMEAS_11_separator_temperature",
        "XMEAS_12_separator_level",
        "XMEAS_13_separator_pressure",
        "XMEAS_14_separator_underflow",
        "XMEAS_15_stripper_level",
        "XMEAS_16_stripper_pressure",
        "XMEAS_17_stripper_underflow",
        "XMEAS_18_stripper_temperature",
        "XMEAS_19_stripper_steam_flow",
        "XMEAS_20_compressor_work",
        "XMEAS_21_reactor_cw_outlet",
        "XMEAS_22_condenser_cw_outlet",
        "XMEAS_23_feed_A",
        "XMEAS_24_feed_B",
        "XMEAS_25_feed_C",
        "XMEAS_26_feed_D",
        "XMEAS_27_feed_E",
        "XMEAS_28_feed_F",
        "XMEAS_29_purge_A",
        "XMEAS_30_purge_B",
        "XMEAS_31_purge_C",
        "XMEAS_32_purge_D",
        "XMEAS_33_purge_E",
        "XMEAS_34_purge_F",
        "XMEAS_35_purge_G",
        "XMEAS_36_purge_H",
        "XMEAS_37_product_D",
        "XMEAS_38_product_E",
        "XMEAS_39_product_F",
        "XMEAS_40_product_G",
        "XMEAS_41_product_H",
        "XMV_1_D_feed_flow",
        "XMV_2_E_feed_flow",
        "XMV_3_A_feed_flow",
        "XMV_4_total_feed_flow",
        "XMV_5_compressor_recycle",
        "XMV_6_purge_valve",
        "XMV_7_separator_underflow",
        "XMV_8_stripper_underflow",
        "XMV_9_stripper_steam",
        "XMV_10_reactor_cw_flow",
        "XMV_11_condenser_cw_flow",
        "XMV_12_agitator_speed",
    ]
}
