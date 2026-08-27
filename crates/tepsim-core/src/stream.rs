//! The thirteen internal process streams.
//!
//! # The trap
//!
//! **The Fortran's stream indices are not the stream numbers in the paper.**
//! `FTM(1)` is the D feed, which Downs and Vogel call stream 2. `FTM(3)` is the
//! A feed, which they call stream 1. Nothing in the source says so; the mapping
//! has to be recovered by matching each `FTM` index against the manipulated
//! variable that drives it and the measurement that reports it.
//!
//! Every reimplementation of TEP has to rediscover this, and getting it wrong
//! produces a plant that runs, looks plausible, and is wired up incorrectly.
//! So it is a typed enum with the correspondence stated on every variant.
//!
//! | Internal | Paper | Stream                                        |
//! |----------|-------|-----------------------------------------------|
//! | 1        | 2     | D feed                                        |
//! | 2        | 3     | E feed                                        |
//! | 3        | 1     | A feed                                        |
//! | 4        | 4     | A and C feed                                  |
//! | 5        | 5     | Stripper overhead vapour to the mixing zone   |
//! | 6        | 6     | Mixing zone outlet to the reactor             |
//! | 7        | 6     | Reactor inlet, an alias of 6                  |
//! | 8        | 7     | Reactor outlet to the condenser and separator |
//! | 9        | 8     | Separator vapour through the compressor       |
//! | 10       | 9     | Purge                                         |
//! | 11       | 10    | Separator liquid underflow to the stripper    |
//! | 12       | none  | Stripper liquid downflow, internal only       |
//! | 13       | 11    | Product                                       |
//!
//! # How the mapping was established
//!
//! Not from the paper, and not from another port. From the Fortran itself:
//!
//! - `teprob.f:564` drives `FTM(1)` from valve 1, and `XMV(1)` is documented as
//!   "D Feed Flow (stream 2)".
//! - `teprob.f:566` gates `FTM(3)` on `IDV(6)`, documented as "A Feed Loss
//!   (Stream 1)".
//! - `teprob.f:567` scales `FTM(4)` by `IDV(7)`, "C Header Pressure Loss
//!   (Stream 4)".
//! - `teprob.f:687` reports `FTM(10)` as `XMEAS(10)`, "Purge Rate (stream 9)".
//! - `teprob.f:682` reports `FTM(9)` as `XMEAS(5)`, "Recycle Flow (stream 8)".
//! - `teprob.f:569` drives `FTM(11)` from valve 7, "Separator Pot Liquid Flow
//!   (stream 10)".
//! - `teprob.f:570` drives `FTM(13)` from valve 8, "Stripper Liquid Product
//!   Flow (stream 11)".

/// A process stream, named by its role rather than by either numbering.
///
/// The discriminant is the Fortran's one-based index, so `Stream::Purge as
/// usize` is 10 and reads `FTM(10)`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(usize)]
pub enum Stream {
    /// `FTM(1)`, paper stream 2. Pure D, driven by `XMV(1)`.
    DFeed = 1,
    /// `FTM(2)`, paper stream 3. Pure E, driven by `XMV(2)`.
    EFeed = 2,
    /// `FTM(3)`, paper stream 1. Pure A, driven by `XMV(3)`. Lost when
    /// `IDV(6)` is active.
    AFeed = 3,
    /// `FTM(4)`, paper stream 4. Mixed A, B and C, driven by `XMV(4)`. The
    /// target of the two composition step disturbances and of `IDV(7)`.
    AcFeed = 4,
    /// `FTM(5)`, paper stream 5. Stripper overhead vapour returning to the
    /// mixing zone.
    StripperOverhead = 5,
    /// `FTM(6)`, paper stream 6. Mixing zone outlet, driven by the pressure
    /// difference against the reactor rather than by a valve.
    MixingZoneOutlet = 6,
    /// `FTM(7)`, paper stream 6. An alias of [`Stream::MixingZoneOutlet`]:
    /// `teprob.f:655-661` copies flow, enthalpy, temperature and composition
    /// across wholesale.
    ReactorInlet = 7,
    /// `FTM(8)`, paper stream 7. Reactor effluent to the condenser and
    /// separator. Throttled by `IDV(20)`.
    ReactorOutlet = 8,
    /// `FTM(9)`, paper stream 8. Separator vapour through the compressor and
    /// back to the mixing zone, less whatever the recycle valve bleeds off.
    Recycle = 9,
    /// `FTM(10)`, paper stream 9. Purge to atmosphere, driven by `XMV(6)`. The
    /// only exit for the inert B.
    Purge = 10,
    /// `FTM(11)`, paper stream 10. Separator liquid to the stripper, driven by
    /// `XMV(7)`.
    SeparatorUnderflow = 11,
    /// `FTM(12)`, no paper number. The liquid that does not strip out, from
    /// `teprob.f:643`. Internal to the stripper model.
    StripperDownflow = 12,
    /// `FTM(13)`, paper stream 11. Product, driven by `XMV(8)`.
    Product = 13,
}

impl Stream {
    /// How many streams the Fortran carries.
    pub const COUNT: usize = 13;

    /// All thirteen, in Fortran order.
    pub const ALL: [Stream; Self::COUNT] = [
        Stream::DFeed,
        Stream::EFeed,
        Stream::AFeed,
        Stream::AcFeed,
        Stream::StripperOverhead,
        Stream::MixingZoneOutlet,
        Stream::ReactorInlet,
        Stream::ReactorOutlet,
        Stream::Recycle,
        Stream::Purge,
        Stream::SeparatorUnderflow,
        Stream::StripperDownflow,
        Stream::Product,
    ];

    /// The Fortran's one-based index, as in `FTM(i)`.
    #[must_use]
    pub const fn fortran_index(self) -> usize {
        self as usize
    }

    /// The zero-based index for Rust arrays.
    #[must_use]
    pub const fn index(self) -> usize {
        self as usize - 1
    }

    /// The stream number used by Downs and Vogel, where one exists.
    ///
    /// [`Stream::StripperDownflow`] has none: it is internal to the stripper
    /// model and does not appear in the paper's flowsheet.
    #[must_use]
    pub const fn paper_number(self) -> Option<u8> {
        Some(match self {
            Stream::DFeed => 2,
            Stream::EFeed => 3,
            Stream::AFeed => 1,
            Stream::AcFeed => 4,
            Stream::StripperOverhead => 5,
            Stream::MixingZoneOutlet | Stream::ReactorInlet => 6,
            Stream::ReactorOutlet => 7,
            Stream::Recycle => 8,
            Stream::Purge => 9,
            Stream::SeparatorUnderflow => 10,
            Stream::StripperDownflow => return None,
            Stream::Product => 11,
        })
    }

    /// A short human-readable name.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Stream::DFeed => "D feed",
            Stream::EFeed => "E feed",
            Stream::AFeed => "A feed",
            Stream::AcFeed => "A and C feed",
            Stream::StripperOverhead => "stripper overhead",
            Stream::MixingZoneOutlet => "mixing zone outlet",
            Stream::ReactorInlet => "reactor inlet",
            Stream::ReactorOutlet => "reactor outlet",
            Stream::Recycle => "recycle",
            Stream::Purge => "purge",
            Stream::SeparatorUnderflow => "separator underflow",
            Stream::StripperDownflow => "stripper downflow",
            Stream::Product => "product",
        }
    }

    /// Recover a stream from its Fortran index.
    #[must_use]
    pub const fn from_fortran_index(index: usize) -> Option<Self> {
        Some(match index {
            1 => Stream::DFeed,
            2 => Stream::EFeed,
            3 => Stream::AFeed,
            4 => Stream::AcFeed,
            5 => Stream::StripperOverhead,
            6 => Stream::MixingZoneOutlet,
            7 => Stream::ReactorInlet,
            8 => Stream::ReactorOutlet,
            9 => Stream::Recycle,
            10 => Stream::Purge,
            11 => Stream::SeparatorUnderflow,
            12 => Stream::StripperDownflow,
            13 => Stream::Product,
            _ => return None,
        })
    }
}

/// Thirteen values, one per stream, indexed by [`Stream`] rather than by a
/// bare integer.
///
/// Layout matches the Fortran array, so this can mirror `FTM`, `XMWS`, `TST`
/// or `HST` directly. The index arithmetic is the whole point: `Stream` is
/// one-based to match the listing, and getting that off by one is the mistake
/// the module documentation exists to prevent.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
#[repr(transparent)]
pub struct ByStream<T>([T; Stream::COUNT]);

impl<T> ByStream<T> {
    /// Wrap thirteen values given in Fortran order.
    pub const fn new(values: [T; Stream::COUNT]) -> Self {
        Self(values)
    }

    /// The underlying array, in Fortran order.
    pub const fn as_array(&self) -> &[T; Stream::COUNT] {
        &self.0
    }

    /// The underlying array, mutably.
    pub const fn as_mut_array(&mut self) -> &mut [T; Stream::COUNT] {
        &mut self.0
    }

    /// Iterate values in Fortran order.
    pub fn iter(&self) -> core::slice::Iter<'_, T> {
        self.0.iter()
    }

    /// Iterate as `(stream, value)` pairs.
    pub fn enumerate(&self) -> impl Iterator<Item = (Stream, &T)> {
        Stream::ALL.into_iter().zip(self.0.iter())
    }
}

impl<T> core::ops::Index<Stream> for ByStream<T> {
    type Output = T;

    fn index(&self, stream: Stream) -> &T {
        &self.0[stream.index()]
    }
}

impl<T> core::ops::IndexMut<Stream> for ByStream<T> {
    fn index_mut(&mut self, stream: Stream) -> &mut T {
        &mut self.0[stream.index()]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fortran_indices_are_one_based_and_dense() {
        for (i, stream) in Stream::ALL.into_iter().enumerate() {
            assert_eq!(stream.fortran_index(), i + 1);
            assert_eq!(stream.index(), i);
            assert_eq!(Stream::from_fortran_index(i + 1), Some(stream));
        }
        assert_eq!(Stream::from_fortran_index(0), None);
        assert_eq!(Stream::from_fortran_index(14), None);
    }

    /// The whole reason this module exists. Three of the four feeds have an
    /// internal index that differs from their paper number, and the two that
    /// are easiest to confuse are A and D.
    #[test]
    fn the_feed_indices_really_do_disagree_with_the_paper() {
        assert_eq!(Stream::DFeed.fortran_index(), 1);
        assert_eq!(Stream::DFeed.paper_number(), Some(2));

        assert_eq!(Stream::AFeed.fortran_index(), 3);
        assert_eq!(Stream::AFeed.paper_number(), Some(1));

        assert_eq!(Stream::EFeed.fortran_index(), 2);
        assert_eq!(Stream::EFeed.paper_number(), Some(3));

        // Only the A and C feed happens to agree.
        assert_eq!(Stream::AcFeed.fortran_index(), 4);
        assert_eq!(Stream::AcFeed.paper_number(), Some(4));
    }

    #[test]
    fn the_reactor_inlet_is_an_alias_of_the_mixing_zone_outlet() {
        assert_eq!(
            Stream::ReactorInlet.paper_number(),
            Stream::MixingZoneOutlet.paper_number()
        );
        assert_eq!(Stream::ReactorInlet.paper_number(), Some(6));
    }

    #[test]
    fn only_the_stripper_downflow_has_no_paper_number() {
        let unnumbered: alloc::vec::Vec<_> = Stream::ALL
            .into_iter()
            .filter(|s| s.paper_number().is_none())
            .collect();
        assert_eq!(unnumbered, alloc::vec![Stream::StripperDownflow]);
    }

    /// `ByStream` must index by the Fortran number, not by the enum's
    /// declaration order. `Stream::Product` is 13, so it must reach slot 12.
    #[test]
    fn by_stream_indexes_by_the_fortran_number() {
        let mut values = ByStream::new([0.0_f64; Stream::COUNT]);
        for stream in Stream::ALL {
            values[stream] = stream.fortran_index() as f64;
        }
        for (slot, value) in values.as_array().iter().enumerate() {
            assert!((*value - (slot + 1) as f64).abs() < f64::EPSILON);
        }
        assert!((values[Stream::Product] - 13.0).abs() < f64::EPSILON);
    }

    /// Paper numbers 1 through 11 must each be covered exactly once, except 6
    /// which the alias makes appear twice.
    #[test]
    fn every_paper_stream_is_accounted_for() {
        let mut seen = [0_u8; 12];
        for stream in Stream::ALL {
            if let Some(n) = stream.paper_number() {
                seen[n as usize] += 1;
            }
        }
        assert_eq!(seen[0], 0, "there is no paper stream 0");
        for (n, count) in seen.iter().enumerate().skip(1).take(11) {
            // Stream 6 appears twice: the mixing zone outlet and its alias.
            let expected = if n == 6 { 2 } else { 1 };
            assert_eq!(
                *count, expected,
                "paper stream {n} is covered {count} times, expected {expected}"
            );
        }
    }
}
