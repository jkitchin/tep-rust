//! Column metadata: what a browser needs to label an axis or a CSV header.
//!
//! Two vocabularies, and they are not interchangeable.
//!
//! [`column_ids`] are `tepsim::channel_names`, the short stable identifiers the
//! rest of the project uses for the 53 channels: `XMEAS_7_reactor_pressure`.
//! They belong in a file header and in code, because they never change and they
//! sort.
//!
//! [`column_labels`] are the descriptions from the `teprob.f` header table,
//! qualified by their Fortran index, which is what belongs on a plot legend.
//! They come from `tepsim_core::variables`, which an integration test in the
//! core compares against the original's comments, so a mistyped unit fails a
//! build rather than propagating into every figure for years.
//!
//! Both lists are one longer than the 53 channels, because a packed row leads
//! with the sample's time. That column is added here rather than in the
//! facade's tables, which describe the plant and know nothing about how these
//! bindings pack a row.

use tepsim::channel_names;
use tepsim::run::{CHANNELS, MANIPULATED as MANIPULATED_COUNT, MEASUREMENTS as MEASUREMENT_COUNT};
use tepsim::tepsim_core::variables::{MANIPULATED, MEASUREMENTS, MeasIndex};

use crate::runner::ROW_WIDTH;

/// The identifier of the leading column of every packed row.
pub const TIME_COLUMN: &str = "time_hours";

/// The unit of the leading column of every packed row.
pub const TIME_UNIT: &str = "hr";

/// Short stable identifiers for the columns of a packed row.
///
/// `time_hours`, then `tepsim::channel_names` in order. For a CSV header, a
/// data frame, or anything a program will read back.
#[must_use]
pub fn column_ids() -> Vec<String> {
    let mut ids = Vec::with_capacity(ROW_WIDTH);
    ids.push(TIME_COLUMN.to_string());
    ids.extend(channel_names().iter().map(|name| (*name).to_string()));
    ids
}

/// Human-readable labels for the columns of a packed row.
///
/// The `teprob.f` descriptions, each qualified by its Fortran index. The
/// qualification is not decoration: "Component A" is the description of both
/// `XMEAS(23)` and `XMEAS(29)`, and a legend that showed the description alone
/// would be ambiguous in nineteen places.
#[must_use]
pub fn column_labels() -> Vec<String> {
    let mut labels = Vec::with_capacity(ROW_WIDTH);
    labels.push("Time".to_string());
    for info in &MEASUREMENTS {
        labels.push(format!("XMEAS({}) {}", info.index, info.description));
    }
    for info in &MANIPULATED {
        labels.push(format!("XMV({}) {}", info.index, info.description));
    }
    labels
}

/// Units for the columns of a packed row, aligned with [`column_ids`].
///
/// Every manipulated variable is a valve position in percent of span,
/// `teprob.f:94-107`, the agitator included: the original reports its speed on
/// the same scale as the rest.
#[must_use]
pub fn column_units() -> Vec<String> {
    let mut units = Vec::with_capacity(ROW_WIDTH);
    units.push(TIME_UNIT.to_string());
    for info in &MEASUREMENTS {
        units.push(info.unit.fortran_spelling().to_string());
    }
    for _ in &MANIPULATED {
        units.push("%".to_string());
    }
    units
}

/// Zero-based offsets into a packed row of the analyser-sampled measurements,
/// `XMEAS(23..=41)`.
///
/// These hold their value between analyser reports, so a chart should draw them
/// as steps rather than interpolating between samples. The browser needs to
/// know which ones before it has any data.
#[must_use]
pub fn sampled_columns() -> Vec<u32> {
    (1..=MeasIndex::COUNT)
        .filter_map(MeasIndex::new)
        .filter(|m| m.is_sampled())
        // Column 0 is time, so `XMEAS(n)` sits at offset `n`.
        .filter_map(|m| u32::try_from(m.get()).ok())
        .collect()
}

/// How many channels a sample carries: 41 measurements and 12 manipulated.
#[must_use]
pub const fn channel_count() -> usize {
    CHANNELS
}

/// How many measurements the plant reports. `XMEAS(41)`.
#[must_use]
pub const fn measurement_count() -> usize {
    MEASUREMENT_COUNT
}

/// How many manipulated variables the plant accepts. `XMV(12)`.
#[must_use]
pub const fn manipulated_count() -> usize {
    MANIPULATED_COUNT
}
