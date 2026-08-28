//! The eight chemical species, and containers indexed by them.
//!
//! The original addresses species by position in bare `DOUBLE PRECISION`
//! arrays, so `XMW(2)` is B's molecular weight only by convention and nothing
//! stops `XLR(I)` being read with the wrong `I`. Here the index is a type.
//!
//! # The species
//!
//! | Component | Role                                                    |
//! |-----------|---------------------------------------------------------|
//! | A         | Gaseous reactant, fed pure in stream 1 and in stream 4   |
//! | B         | Inert. Enters in stream 4, leaves only through the purge |
//! | C         | Gaseous reactant, fed in stream 4                       |
//! | D         | Reactant fed in stream 2                                |
//! | E         | Reactant fed in stream 3                                |
//! | F         | Byproduct                                               |
//! | G         | Product                                                 |
//! | H         | Product                                                 |
//!
//! The four reactions, from `teprob.f:518-524`:
//!
//! ```text
//! A + C + D -> G        A + C + E -> H        A + E -> F        3D -> 2F
//! ```
//!
//! B takes part in nothing, which is why the purge exists: without it, the
//! inert would accumulate until the plant tripped on pressure.

use core::ops::{Index, IndexMut};

/// One of the eight species.
///
/// The discriminants are the zero-based positions in the Fortran arrays, so
/// `Component::B as usize` is 1 and reads `XMW(2)`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(usize)]
pub enum Component {
    /// Gaseous reactant.
    A = 0,
    /// Inert.
    B = 1,
    /// Gaseous reactant.
    C = 2,
    /// Reactant, fed in stream 2.
    D = 3,
    /// Reactant, fed in stream 3.
    E = 4,
    /// Byproduct.
    F = 5,
    /// Product.
    G = 6,
    /// Product.
    H = 7,
}

impl Component {
    /// How many species there are. The original hard-codes 8 everywhere.
    pub const COUNT: usize = 8;

    /// All eight, in Fortran order.
    pub const ALL: [Component; Self::COUNT] = [
        Component::A,
        Component::B,
        Component::C,
        Component::D,
        Component::E,
        Component::F,
        Component::G,
        Component::H,
    ];

    /// The single-letter name used throughout the TEP literature.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Component::A => "A",
            Component::B => "B",
            Component::C => "C",
            Component::D => "D",
            Component::E => "E",
            Component::F => "F",
            Component::G => "G",
            Component::H => "H",
        }
    }

    /// The zero-based array index. `XMW(i)` in Fortran is `index() + 1`.
    #[must_use]
    pub const fn index(self) -> usize {
        self as usize
    }

    /// The one-based index the Fortran uses.
    #[must_use]
    pub const fn fortran_index(self) -> usize {
        self as usize + 1
    }

    /// Whether this species takes part in any reaction.
    ///
    /// Only B does not, which is the entire reason the plant has a purge.
    #[must_use]
    pub const fn is_inert(self) -> bool {
        matches!(self, Component::B)
    }

    /// Recover a component from a zero-based index.
    #[must_use]
    pub const fn from_index(index: usize) -> Option<Self> {
        Some(match index {
            0 => Component::A,
            1 => Component::B,
            2 => Component::C,
            3 => Component::D,
            4 => Component::E,
            5 => Component::F,
            6 => Component::G,
            7 => Component::H,
            _ => return None,
        })
    }
}

/// Eight values, one per species, indexed by [`Component`] rather than by a
/// bare integer.
///
/// Layout is identical to the Fortran array, so this can mirror a `COMMON`
/// block field directly.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
#[repr(transparent)]
pub struct ByComponent<T>([T; Component::COUNT]);

impl<T> ByComponent<T> {
    /// Wrap eight values given in Fortran order.
    pub const fn new(values: [T; Component::COUNT]) -> Self {
        Self(values)
    }

    /// The underlying array, in Fortran order.
    pub const fn as_array(&self) -> &[T; Component::COUNT] {
        &self.0
    }

    /// The underlying array, mutably.
    pub const fn as_mut_array(&mut self) -> &mut [T; Component::COUNT] {
        &mut self.0
    }

    /// Consume into the underlying array.
    pub fn into_array(self) -> [T; Component::COUNT] {
        self.0
    }

    /// Iterate values in Fortran order.
    pub fn iter(&self) -> core::slice::Iter<'_, T> {
        self.0.iter()
    }

    /// Iterate mutably.
    pub fn iter_mut(&mut self) -> core::slice::IterMut<'_, T> {
        self.0.iter_mut()
    }

    /// Iterate as `(component, value)` pairs.
    pub fn enumerate(&self) -> impl Iterator<Item = (Component, &T)> {
        Component::ALL.into_iter().zip(self.0.iter())
    }
}

impl<T> Index<Component> for ByComponent<T> {
    type Output = T;

    fn index(&self, component: Component) -> &T {
        &self.0[component.index()]
    }
}

impl<T> IndexMut<Component> for ByComponent<T> {
    fn index_mut(&mut self, component: Component) -> &mut T {
        &mut self.0[component.index()]
    }
}

impl<T> From<[T; Component::COUNT]> for ByComponent<T> {
    fn from(values: [T; Component::COUNT]) -> Self {
        Self(values)
    }
}

impl<'a, T> IntoIterator for &'a ByComponent<T> {
    type Item = &'a T;
    type IntoIter = core::slice::Iter<'a, T>;

    fn into_iter(self) -> Self::IntoIter {
        self.0.iter()
    }
}

/// Mole fractions across the eight species.
///
/// # Why the sum check is a debug assertion and not an invariant
///
/// The original's own feed compositions do not sum to one exactly. They are
/// written as single-precision literals like `0.9999` and `0.0001` at
/// `teprob.f:1134-1159`, so widening them to double leaves a residual of about
/// 1.7e-8. Enforcing an exact sum would reject the original's own data, and
/// enforcing it loosely in release builds would cost time in the inner loop for
/// no benefit. So the check runs in debug builds with a tolerance chosen to
/// admit what the original actually does.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
#[repr(transparent)]
pub struct Composition(ByComponent<f64>);

impl Composition {
    /// How far from 1.0 a composition may sum in debug builds.
    ///
    /// Set by the original's single-precision feed literals, not by taste.
    pub const SUM_TOLERANCE: f64 = 1e-6;

    /// Build a composition, checking the sum in debug builds.
    #[must_use]
    pub fn new(fractions: [f64; Component::COUNT]) -> Self {
        let composition = Self(ByComponent::new(fractions));
        debug_assert!(
            composition.sums_to_one(),
            "mole fractions must sum to 1 within {}, got {}",
            Self::SUM_TOLERANCE,
            composition.sum()
        );
        composition
    }

    /// Build without checking. For intermediate states that are not yet
    /// normalised, which the original produces routinely.
    #[must_use]
    pub const fn new_unchecked(fractions: [f64; Component::COUNT]) -> Self {
        Self(ByComponent::new(fractions))
    }

    /// The sum of all eight fractions.
    #[must_use]
    pub fn sum(&self) -> f64 {
        // Summed in Fortran order. Reassociating would change the last bits.
        let mut total = 0.0;
        for value in self.0.iter() {
            total += *value;
        }
        total
    }

    /// Whether the sum is within [`Composition::SUM_TOLERANCE`] of one.
    #[must_use]
    pub fn sums_to_one(&self) -> bool {
        (self.sum() - 1.0).abs() <= Self::SUM_TOLERANCE
    }

    /// The underlying per-component values.
    #[must_use]
    pub const fn fractions(&self) -> &ByComponent<f64> {
        &self.0
    }
}

impl Index<Component> for Composition {
    type Output = f64;

    fn index(&self, component: Component) -> &f64 {
        &self.0[component]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testing::assert_exact;

    #[test]
    fn indices_match_the_fortran_positions() {
        assert_eq!(Component::A.index(), 0);
        assert_eq!(Component::A.fortran_index(), 1);
        assert_eq!(Component::H.index(), 7);
        assert_eq!(Component::H.fortran_index(), 8);
        for (i, component) in Component::ALL.into_iter().enumerate() {
            assert_eq!(component.index(), i);
            assert_eq!(Component::from_index(i), Some(component));
        }
        assert_eq!(Component::from_index(8), None);
    }

    #[test]
    fn only_b_is_inert() {
        let inert: alloc::vec::Vec<_> = Component::ALL
            .into_iter()
            .filter(|c| c.is_inert())
            .collect();
        assert_eq!(inert, alloc::vec![Component::B]);
    }

    #[test]
    fn indexing_reaches_the_right_slot() {
        let mut values = ByComponent::new([0.0; 8]);
        values[Component::E] = 42.0;
        assert_exact(values.as_array()[4], 42.0, "raw slot 4");
        assert_exact(values[Component::E], 42.0, "indexed by E");
        assert_exact(values[Component::D], 0.0, "D untouched");
    }

    #[test]
    fn enumerate_pairs_components_with_their_values() {
        let values = ByComponent::new([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]);
        let collected: alloc::vec::Vec<_> = values.enumerate().map(|(c, v)| (c, *v)).collect();
        assert_eq!(collected[0].0, Component::A);
        assert_exact(collected[0].1, 1.0, "first value");
        assert_eq!(collected[7].0, Component::H);
        assert_exact(collected[7].1, 8.0, "last value");
    }

    #[test]
    fn a_normalised_composition_is_accepted() {
        let c = Composition::new([0.125; 8]);
        assert!(c.sums_to_one());
    }

    /// The tolerance exists for the original's own data, so check it admits it.
    #[test]
    fn the_originals_single_precision_feed_composition_is_accepted() {
        // teprob.f:1134-1137, stream 1: 0.0001 of B and 0.9999 of D, both
        // written without a D suffix and so stored at f32 precision.
        let mut fractions = [0.0; 8];
        fractions[Component::B.index()] = f64::from(0.0001_f32);
        fractions[Component::D.index()] = f64::from(0.9999_f32);
        let c = Composition::new_unchecked(fractions);
        assert!(
            c.sums_to_one(),
            "the tolerance must admit the original's own feed data, which sums \
             to {} rather than exactly 1",
            c.sum()
        );
        assert!(
            (c.sum() - 1.0).abs() > 1e-9,
            "and that data really is off by more than a rounding error"
        );
    }

    /// The check is a `debug_assert!`, so it exists only in a debug build and
    /// this test can only run in one.
    ///
    /// Without the `cfg` the test fails in release, where the assertion is
    /// compiled out and `Composition::new` returns normally. That went
    /// unnoticed because `cargo xtask ci` runs the workspace tests in debug;
    /// it surfaced the first time anything ran `cargo test --release` on this
    /// crate.
    #[test]
    #[cfg(debug_assertions)]
    #[should_panic(expected = "must sum to 1")]
    fn a_composition_that_does_not_sum_to_one_is_rejected_in_debug() {
        let _ = Composition::new([0.5; 8]);
    }

    /// And in release it does not panic, which is the other half of the
    /// contract: the check is a development aid and costs nothing shipped.
    #[test]
    #[cfg(not(debug_assertions))]
    fn a_composition_that_does_not_sum_to_one_is_accepted_in_release() {
        let c = Composition::new([0.5; 8]);
        assert!(!c.sums_to_one());
    }

    #[test]
    fn unchecked_construction_skips_the_assertion() {
        let c = Composition::new_unchecked([0.5; 8]);
        assert!(!c.sums_to_one());
    }
}
