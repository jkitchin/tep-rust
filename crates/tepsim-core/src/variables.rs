//! The 41 measurements, 12 manipulated variables and 12 valves, with names and
//! units.
//!
//! Names and units are transcribed verbatim from the header comment tables in
//! `teprob.f:109-168` and `teprob.f:94-107`. Verbatim on purpose: an
//! integration test parses those comments out of the vendored source and
//! compares them against the tables here, so a mistyped unit fails the build
//! rather than propagating into every plot and dataset column for years.

/// The engineering unit a measurement is reported in.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Unit {
    /// Thousand standard cubic metres per hour.
    Kscmh,
    /// Kilograms per hour.
    KgPerHour,
    /// Kilopascals, gauge.
    KPaGauge,
    /// Percent of span.
    Percent,
    /// Degrees Celsius.
    DegC,
    /// Cubic metres per hour.
    M3PerHour,
    /// Kilowatts.
    Kw,
    /// Mole percent, for the composition analysers.
    MolePercent,
}

impl Unit {
    /// Exactly how the unit is spelled in the `teprob.f` header table.
    ///
    /// The parse test compares against this, so it must not be prettified.
    #[must_use]
    pub const fn fortran_spelling(self) -> &'static str {
        match self {
            Unit::Kscmh => "kscmh",
            Unit::KgPerHour => "kg/hr",
            Unit::KPaGauge => "kPa gauge",
            Unit::Percent => "%",
            Unit::DegC => "Deg C",
            Unit::M3PerHour => "m3/hr",
            Unit::Kw => "kW",
            Unit::MolePercent => "Mole %",
        }
    }
}

/// Which analyser produced a sampled composition measurement.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Analyzer {
    /// Reactor feed, paper stream 6. Samples every 0.1 h with 0.1 h dead time.
    ReactorFeed,
    /// Purge gas, paper stream 9. Samples every 0.1 h with 0.1 h dead time.
    PurgeGas,
    /// Product, paper stream 11. Samples every 0.25 h with 0.25 h dead time.
    Product,
}

impl Analyzer {
    /// Sampling interval in hours, `teprob.f:136-137, 148-149, 161-162`.
    #[must_use]
    pub const fn sampling_interval_hours(self) -> f64 {
        match self {
            Analyzer::ReactorFeed | Analyzer::PurgeGas => 0.1,
            Analyzer::Product => 0.25,
        }
    }

    /// Dead time in hours. Equal to the sampling interval for all three, which
    /// is why the original implements it by reporting the previous sample.
    #[must_use]
    pub const fn dead_time_hours(self) -> f64 {
        self.sampling_interval_hours()
    }
}

/// One entry of the measurement table.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MeasurementInfo {
    /// One-based index, as in `XMEAS(i)`.
    pub index: usize,
    /// The description, verbatim from the `teprob.f` header table.
    pub description: &'static str,
    /// The reported unit.
    pub unit: Unit,
    /// The analyser, for the 19 sampled composition measurements.
    pub analyzer: Option<Analyzer>,
}

/// One entry of the manipulated variable table.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ManipulatedInfo {
    /// One-based index, as in `XMV(i)`.
    pub index: usize,
    /// The description, verbatim from the `teprob.f` header table.
    pub description: &'static str,
}

/// A one-based measurement index, `XMEAS(1..=41)`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct MeasIndex(usize);

/// A one-based manipulated variable index, `XMV(1..=12)`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct MvIndex(usize);

/// A one-based valve index. There is one valve per manipulated variable, and
/// its position is a lagged state, `YY(38 + i)`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct ValveId(usize);

impl MeasIndex {
    /// How many measurements there are.
    pub const COUNT: usize = 41;
    /// The first index at which measurements are sampled rather than continuous.
    pub const FIRST_SAMPLED: usize = 23;

    /// Build from a one-based index, or `None` if out of range.
    #[must_use]
    pub const fn new(index: usize) -> Option<Self> {
        if index >= 1 && index <= Self::COUNT {
            Some(Self(index))
        } else {
            None
        }
    }

    /// The one-based index.
    #[must_use]
    pub const fn get(self) -> usize {
        self.0
    }

    /// The zero-based index into a Rust array.
    #[must_use]
    pub const fn zero_based(self) -> usize {
        self.0 - 1
    }

    /// The table entry for this measurement.
    #[must_use]
    pub const fn info(self) -> &'static MeasurementInfo {
        &MEASUREMENTS[self.0 - 1]
    }

    /// Whether this measurement comes from an analyser rather than continuously.
    #[must_use]
    pub const fn is_sampled(self) -> bool {
        self.0 >= Self::FIRST_SAMPLED
    }
}

impl MvIndex {
    /// How many manipulated variables there are.
    pub const COUNT: usize = 12;

    /// Build from a one-based index, or `None` if out of range.
    #[must_use]
    pub const fn new(index: usize) -> Option<Self> {
        if index >= 1 && index <= Self::COUNT {
            Some(Self(index))
        } else {
            None
        }
    }

    /// The one-based index.
    #[must_use]
    pub const fn get(self) -> usize {
        self.0
    }

    /// The zero-based index into a Rust array.
    #[must_use]
    pub const fn zero_based(self) -> usize {
        self.0 - 1
    }

    /// The table entry for this variable.
    #[must_use]
    pub const fn info(self) -> &'static ManipulatedInfo {
        &MANIPULATED[self.0 - 1]
    }

    /// The valve this variable drives. One each, same numbering.
    #[must_use]
    pub const fn valve(self) -> ValveId {
        ValveId(self.0)
    }
}

impl ValveId {
    /// How many valves there are.
    pub const COUNT: usize = 12;

    /// Build from a one-based index, or `None` if out of range.
    #[must_use]
    pub const fn new(index: usize) -> Option<Self> {
        if index >= 1 && index <= Self::COUNT {
            Some(Self(index))
        } else {
            None
        }
    }

    /// The one-based index.
    #[must_use]
    pub const fn get(self) -> usize {
        self.0
    }

    /// The state index whose value is this valve's position.
    ///
    /// `teprob.f:437` reads `VPOS(I) = YY(I + 38)`, so valve 1 is state 39.
    #[must_use]
    pub const fn state_index(self) -> usize {
        self.0 + 38
    }
}

/// The 41 measurements, verbatim from `teprob.f:109-168`.
//
// Claims the whole documentation block. These tables are a transcription of it,
// and `tests/header_tables.rs` proves the transcription is exact, so the block
// is genuinely accounted for rather than merely referenced.
//
// @port teprob.f:109-168
pub const MEASUREMENTS: [MeasurementInfo; MeasIndex::COUNT] = {
    const fn cont(index: usize, description: &'static str, unit: Unit) -> MeasurementInfo {
        MeasurementInfo {
            index,
            description,
            unit,
            analyzer: None,
        }
    }
    const fn samp(index: usize, description: &'static str, analyzer: Analyzer) -> MeasurementInfo {
        MeasurementInfo {
            index,
            description,
            unit: Unit::MolePercent,
            analyzer: Some(analyzer),
        }
    }
    [
        cont(1, "A Feed (stream 1)", Unit::Kscmh),
        cont(2, "D Feed (stream 2)", Unit::KgPerHour),
        cont(3, "E Feed (stream 3)", Unit::KgPerHour),
        cont(4, "A and C Feed (stream 4)", Unit::Kscmh),
        cont(5, "Recycle Flow (stream 8)", Unit::Kscmh),
        cont(6, "Reactor Feed Rate (stream 6)", Unit::Kscmh),
        cont(7, "Reactor Pressure", Unit::KPaGauge),
        cont(8, "Reactor Level", Unit::Percent),
        cont(9, "Reactor Temperature", Unit::DegC),
        cont(10, "Purge Rate (stream 9)", Unit::Kscmh),
        cont(11, "Product Sep Temp", Unit::DegC),
        cont(12, "Product Sep Level", Unit::Percent),
        cont(13, "Prod Sep Pressure", Unit::KPaGauge),
        cont(14, "Prod Sep Underflow (stream 10)", Unit::M3PerHour),
        cont(15, "Stripper Level", Unit::Percent),
        cont(16, "Stripper Pressure", Unit::KPaGauge),
        cont(17, "Stripper Underflow (stream 11)", Unit::M3PerHour),
        cont(18, "Stripper Temperature", Unit::DegC),
        cont(19, "Stripper Steam Flow", Unit::KgPerHour),
        cont(20, "Compressor Work", Unit::Kw),
        cont(21, "Reactor Cooling Water Outlet Temp", Unit::DegC),
        cont(22, "Separator Cooling Water Outlet Temp", Unit::DegC),
        samp(23, "Component A", Analyzer::ReactorFeed),
        samp(24, "Component B", Analyzer::ReactorFeed),
        samp(25, "Component C", Analyzer::ReactorFeed),
        samp(26, "Component D", Analyzer::ReactorFeed),
        samp(27, "Component E", Analyzer::ReactorFeed),
        samp(28, "Component F", Analyzer::ReactorFeed),
        samp(29, "Component A", Analyzer::PurgeGas),
        samp(30, "Component B", Analyzer::PurgeGas),
        samp(31, "Component C", Analyzer::PurgeGas),
        samp(32, "Component D", Analyzer::PurgeGas),
        samp(33, "Component E", Analyzer::PurgeGas),
        samp(34, "Component F", Analyzer::PurgeGas),
        samp(35, "Component G", Analyzer::PurgeGas),
        samp(36, "Component H", Analyzer::PurgeGas),
        samp(37, "Component D", Analyzer::Product),
        samp(38, "Component E", Analyzer::Product),
        samp(39, "Component F", Analyzer::Product),
        samp(40, "Component G", Analyzer::Product),
        samp(41, "Component H", Analyzer::Product),
    ]
};

/// The 12 manipulated variables, verbatim from `teprob.f:94-107`.
///
/// The first three carry a "(Corrected Order)" annotation upstream, recording
/// that the 1991 revision fixed their documented order. The annotation is not
/// part of the name and is dropped here; the parse test accounts for it.
//
// @port teprob.f:94-107
pub const MANIPULATED: [ManipulatedInfo; MvIndex::COUNT] = {
    const fn mv(index: usize, description: &'static str) -> ManipulatedInfo {
        ManipulatedInfo { index, description }
    }
    [
        mv(1, "D Feed Flow (stream 2)"),
        mv(2, "E Feed Flow (stream 3)"),
        mv(3, "A Feed Flow (stream 1)"),
        mv(4, "A and C Feed Flow (stream 4)"),
        mv(5, "Compressor Recycle Valve"),
        mv(6, "Purge Valve (stream 9)"),
        mv(7, "Separator Pot Liquid Flow (stream 10)"),
        mv(8, "Stripper Liquid Product Flow (stream 11)"),
        mv(9, "Stripper Steam Valve"),
        mv(10, "Reactor Cooling Water Flow"),
        mv(11, "Condenser Cooling Water Flow"),
        mv(12, "Agitator Speed"),
    ]
};

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testing::assert_exact;

    #[test]
    fn measurement_indices_are_dense_and_one_based() {
        for (i, info) in MEASUREMENTS.iter().enumerate() {
            assert_eq!(info.index, i + 1);
        }
        assert!(MeasIndex::new(0).is_none());
        assert!(MeasIndex::new(42).is_none());
        assert_eq!(MeasIndex::new(41).map(MeasIndex::get), Some(41));
    }

    #[test]
    fn manipulated_indices_are_dense_and_one_based() {
        for (i, info) in MANIPULATED.iter().enumerate() {
            assert_eq!(info.index, i + 1);
        }
        assert!(MvIndex::new(0).is_none());
        assert!(MvIndex::new(13).is_none());
    }

    #[test]
    fn the_sampled_measurements_are_exactly_23_through_41() {
        for i in 1..=MeasIndex::COUNT {
            let m = MeasIndex::new(i).expect("in range");
            assert_eq!(
                m.is_sampled(),
                m.info().analyzer.is_some(),
                "XMEAS({i}) disagrees about whether it is sampled"
            );
        }
        let sampled = (1..=MeasIndex::COUNT)
            .filter(|i| MeasIndex::new(*i).expect("in range").is_sampled())
            .count();
        assert_eq!(sampled, 19, "there are 19 composition measurements");
    }

    #[test]
    fn analyser_timing_matches_the_header_table() {
        assert_exact(
            Analyzer::ReactorFeed.sampling_interval_hours(),
            0.1,
            "reactor feed",
        );
        assert_exact(
            Analyzer::PurgeGas.sampling_interval_hours(),
            0.1,
            "purge gas",
        );
        assert_exact(Analyzer::Product.sampling_interval_hours(), 0.25, "product");
        for a in [Analyzer::ReactorFeed, Analyzer::PurgeGas, Analyzer::Product] {
            assert_exact(
                a.dead_time_hours(),
                a.sampling_interval_hours(),
                "dead time",
            );
        }
    }

    #[test]
    fn each_analyser_covers_the_species_the_header_lists() {
        let count = |a: Analyzer| {
            MEASUREMENTS
                .iter()
                .filter(|m| m.analyzer == Some(a))
                .count()
        };
        assert_eq!(count(Analyzer::ReactorFeed), 6, "A through F");
        assert_eq!(count(Analyzer::PurgeGas), 8, "A through H");
        assert_eq!(count(Analyzer::Product), 5, "D through H");
    }

    /// `teprob.f:437` reads `VPOS(I) = YY(I + 38)`.
    #[test]
    fn valve_positions_are_states_39_through_50() {
        assert_eq!(ValveId::new(1).expect("valve 1").state_index(), 39);
        assert_eq!(ValveId::new(12).expect("valve 12").state_index(), 50);
        assert!(ValveId::new(13).is_none());
    }

    #[test]
    fn every_manipulated_variable_drives_the_valve_of_the_same_number() {
        for i in 1..=MvIndex::COUNT {
            let mv = MvIndex::new(i).expect("in range");
            assert_eq!(mv.valve().get(), i);
        }
    }
}
