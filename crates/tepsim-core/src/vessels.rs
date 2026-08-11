//! Unpacking the state vector into what the balance equations actually need.
//!
//! Ported from `teprob.f:417-472`. This is the first step of every derivative
//! evaluation: the fifty integrated states become per-vessel inventories, mole
//! fractions, specific energies, temperatures, densities and volumes.
//!
//! # The equations
//!
//! For each vessel, writing \\(n_i\\) for the component holdups:
//!
//! \\[
//!   N = \\sum_i n_i, \\qquad x_i = \\frac{n_i}{N}, \\qquad e = \\frac{E}{N}
//! \\]
//!
//! The temperature is then whatever makes the mixture's specific enthalpy equal
//! \\(e\\), solved by [`crate::thermo::temperature_from_enthalpy`], and the
//! liquid density follows from [`crate::thermo::liquid_density`]. The volume is
//! \\(V = N / \\rho\\).
//!
//! # Variables
//!
//! | Fortran | Here | Meaning |
//! |---|---|---|
//! | `UCLR`, `UCLS`, `UCLC` | [`Liquid::moles`] | liquid component holdups |
//! | `UCVR`, `UCVS` | [`Unpacked::reactor_vapour`] etc. | non-condensible vapour holdups |
//! | `UCVV` | [`Vapour::moles`] | mixing zone holdups |
//! | `UTLR`, `UTLS`, `UTLC`, `UTVV` | `total` | summed holdups |
//! | `XLR`, `XLS`, `XLC`, `XVV` | `fractions` | mole fractions |
//! | `ESR`, `ESS`, `ESC`, `ESV` | `specific_energy` | energy per mole |
//! | `TCR`, `TCS`, `TCC`, `TCV` | `celsius` | temperature |
//! | `DLR`, `DLS`, `DLC` | [`Liquid::density`] | molar density |
//! | `VLR`, `VLS`, `VLC` | [`Liquid::volume`] | liquid volume |
//!
//! # The temperatures are warm-started, and that is not a detail
//!
//! `TESUB2` takes its temperature argument as both the initial guess and the
//! result (`teprob.f:1432`, `1438`), and the four call sites at
//! `teprob.f:460-465` pass `TCR`, `TCS`, `TCC` and `TCV` straight out of
//! `COMMON`. Each evaluation therefore starts its Newton solves from the
//! previous evaluation's answers, and since Newton stops on a step below 1e-12
//! the converged value depends on where it started.
//!
//! So the four temperatures are *state*, not derived quantities, and
//! [`TemperatureSeeds`] carries them. B-0015 measured the cost of getting this
//! wrong: seeding them from a different point on the nominal trajectory moves
//! up to 21 of the 50 derivatives. A port that solved from a fixed guess would
//! be tidier and would not be bit-exact.
//!
//! # `273.15` is single precision here too
//!
//! All three occurrences in this range (`teprob.f:461`, `463`, `466`) are
//! written without a `D` suffix, exactly like the one at `teprob.f:1411` that
//! delta D-001 turns on. [`crate::thermo::ABSOLUTE_ZERO_OFFSET`] is therefore
//! the right constant for all of them.

use crate::component::{ByComponent, Component, Composition};
use crate::state::State;
use crate::thermo::{
    ABSOLUTE_ZERO_OFFSET, EnergyBasis, TemperatureError, liquid_density, temperature_from_enthalpy,
};

/// The four Newton warm-start temperatures, in degrees Celsius.
///
/// Carried between evaluations because the original carries them. See the
/// module documentation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TemperatureSeeds {
    /// `TCR`, the reactor liquid temperature.
    pub reactor: f64,
    /// `TCS`, the separator liquid temperature.
    pub separator: f64,
    /// `TCC`, the stripper liquid temperature.
    pub stripper: f64,
    /// `TCV`, the mixing zone vapour temperature.
    pub mixing: f64,
}

impl Default for TemperatureSeeds {
    /// The nominal operating point, which is what `TEINIT` leaves in `COMMON`
    /// before the first evaluation.
    fn default() -> Self {
        Self {
            reactor: 120.4,
            separator: 80.109,
            stripper: 65.731,
            mixing: 86.120,
        }
    }
}

/// A vessel's liquid phase, unpacked.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Liquid {
    /// Component holdups. Slots A, B and C are zero for the reactor and the
    /// separator (`teprob.f:420-421`): the non-condensibles form no liquid.
    pub moles: ByComponent<f64>,
    /// Total moles, summed in Fortran order.
    pub total: f64,
    /// Mole fractions.
    pub fractions: Composition,
    /// Internal energy per mole.
    pub specific_energy: f64,
    /// Temperature in degrees Celsius.
    pub celsius: f64,
    /// Molar density, lbmol per cubic foot.
    pub density: f64,
    /// Volume, cubic feet.
    pub volume: f64,
}

impl Liquid {
    /// The temperature in kelvin, as `teprob.f:461` computes it.
    ///
    /// Uses the single-precision offset; see the module documentation.
    #[must_use]
    pub fn kelvin(&self) -> f64 {
        self.celsius + ABSOLUTE_ZERO_OFFSET
    }
}

/// A vessel's vapour phase, unpacked. Only the mixing zone has one that comes
/// straight from the state.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Vapour {
    /// Component holdups.
    pub moles: ByComponent<f64>,
    /// Total moles.
    pub total: f64,
    /// Mole fractions.
    pub fractions: Composition,
    /// Internal energy per mole.
    pub specific_energy: f64,
    /// Temperature in degrees Celsius.
    pub celsius: f64,
}

impl Vapour {
    /// The temperature in kelvin, as `teprob.f:466` computes it.
    #[must_use]
    pub fn kelvin(&self) -> f64 {
        self.celsius + ABSOLUTE_ZERO_OFFSET
    }
}

/// Everything `teprob.f:417-472` produces from the state vector.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Unpacked {
    /// Reactor liquid.
    pub reactor: Liquid,
    /// Separator liquid.
    pub separator: Liquid,
    /// Stripper liquid.
    pub stripper: Liquid,
    /// Mixing zone vapour.
    pub mixing: Vapour,
    /// `UCVR`: the reactor's non-condensible vapour holdups, A/B/C only.
    ///
    /// The D-H slots are zero here and are filled in from the vapour-liquid
    /// equilibrium later (`teprob.f:500-501`), which is B-0018's job.
    pub reactor_vapour: ByComponent<f64>,
    /// `UCVS`: the separator's, likewise.
    pub separator_vapour: ByComponent<f64>,
    /// The temperatures this evaluation converged on, to seed the next one.
    pub seeds: TemperatureSeeds,
}

/// Sum in Fortran order. `teprob.f:443-448` accumulates from zero in index
/// order, and reassociating would change the last bits.
fn total(moles: &ByComponent<f64>) -> f64 {
    let mut sum = 0.0;
    for value in moles.iter() {
        sum += *value;
    }
    sum
}

/// Divide each holdup by the total, in place, as `teprob.f:450-453` does.
fn fractions(moles: &ByComponent<f64>, total: f64) -> Composition {
    let mut out = [0.0; Component::COUNT];
    for (slot, value) in out.iter_mut().zip(moles.iter()) {
        *slot = value / total;
    }
    // Unchecked: the original produces un-normalised intermediates routinely,
    // and an adversarial state can drive a total to zero. A NaN here is a real
    // answer about a degenerate state, not something to assert away.
    Composition::new_unchecked(out)
}

/// Unpack the state into per-vessel quantities.
///
/// `seeds` are the previous evaluation's converged temperatures; see the module
/// documentation for why they are inputs rather than something to derive.
///
/// # Errors
///
/// [`TemperatureError::DidNotConverge`] if any of the four Newton solves runs
/// out of iterations. The original silently returns the guess instead; that is
/// delta D-001, and it has never fired on the physical domain.
// @port teprob.f:417-472
pub fn unpack(y: &State, seeds: TemperatureSeeds) -> Result<Unpacked, TemperatureError> {
    // teprob.f:417-430. The non-condensibles are vapour only, the heavies
    // liquid only, and the split differs by vessel.
    let mut reactor_vapour = ByComponent::new([0.0; Component::COUNT]);
    let mut separator_vapour = ByComponent::new([0.0; Component::COUNT]);
    let mut reactor_liquid = ByComponent::new([0.0; Component::COUNT]);
    let mut separator_liquid = ByComponent::new([0.0; Component::COUNT]);
    for component in Component::ALL {
        let held = y.reactor.moles[component];
        let separated = y.separator.moles[component];
        if component.index() < 3 {
            reactor_vapour[component] = held;
            separator_vapour[component] = separated;
        } else {
            reactor_liquid[component] = held;
            separator_liquid[component] = separated;
        }
    }

    // teprob.f:441-457, then 460-472.
    let liquid =
        |moles: ByComponent<f64>, energy: f64, seed: f64| -> Result<Liquid, TemperatureError> {
            let total = total(&moles);
            let fractions = fractions(&moles, total);
            let specific_energy = energy / total;
            let celsius = temperature_from_enthalpy(
                &fractions,
                seed,
                specific_energy,
                EnergyBasis::LiquidEnthalpy,
            )?;
            let density = liquid_density(&fractions, celsius);
            Ok(Liquid {
                moles,
                total,
                fractions,
                specific_energy,
                celsius,
                density,
                volume: total / density,
            })
        };

    let reactor = liquid(reactor_liquid, y.reactor.energy, seeds.reactor)?;
    let separator = liquid(separator_liquid, y.separator.energy, seeds.separator)?;
    let stripper = liquid(y.stripper.moles, y.stripper.energy, seeds.stripper)?;

    // teprob.f:465. The mixing zone is vapour, so its temperature comes from
    // the internal-energy correlation rather than the liquid one.
    let mixing_total = total(&y.mixing.moles);
    let mixing_fractions = fractions(&y.mixing.moles, mixing_total);
    let mixing_specific_energy = y.mixing.energy / mixing_total;
    let mixing_celsius = temperature_from_enthalpy(
        &mixing_fractions,
        seeds.mixing,
        mixing_specific_energy,
        EnergyBasis::VapourInternalEnergy,
    )?;

    Ok(Unpacked {
        seeds: TemperatureSeeds {
            reactor: reactor.celsius,
            separator: separator.celsius,
            stripper: stripper.celsius,
            mixing: mixing_celsius,
        },
        reactor,
        separator,
        stripper,
        mixing: Vapour {
            moles: y.mixing.moles,
            total: mixing_total,
            fractions: mixing_fractions,
            specific_energy: mixing_specific_energy,
            celsius: mixing_celsius,
        },
        reactor_vapour,
        separator_vapour,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testing::assert_exact;

    /// The phase split is by component index and differs between vessels, which
    /// is the easiest thing here to get backwards.
    #[test]
    fn the_non_condensibles_are_vapour_and_the_heavies_are_liquid() {
        let mut y = State::default();
        for component in Component::ALL {
            y.reactor.moles[component] = 10.0 + component.index() as f64;
            y.separator.moles[component] = 20.0 + component.index() as f64;
            y.stripper.moles[component] = 1.0;
            y.mixing.moles[component] = 1.0;
        }
        y.reactor.energy = 100.0;
        y.separator.energy = 100.0;
        y.stripper.energy = 8.0;
        y.mixing.energy = 8.0;

        let unpacked = unpack(&y, TemperatureSeeds::default()).expect("converges");

        for component in [Component::A, Component::B, Component::C] {
            assert_exact(
                unpacked.reactor_vapour[component],
                y.reactor.moles[component],
                "reactor A/B/C are vapour",
            );
            assert_exact(
                unpacked.reactor.moles[component],
                0.0,
                "reactor A/B/C hold no liquid",
            );
        }
        for component in [Component::D, Component::H] {
            assert_exact(
                unpacked.reactor.moles[component],
                y.reactor.moles[component],
                "reactor D-H are liquid",
            );
            assert_exact(
                unpacked.reactor_vapour[component],
                0.0,
                "reactor D-H vapour comes from the equilibrium, not the state",
            );
        }
        // The stripper and the mixing zone keep all eight.
        assert_exact(
            unpacked.stripper.moles[Component::A],
            1.0,
            "the stripper is all liquid, including A",
        );
        assert_exact(
            unpacked.mixing.moles[Component::A],
            1.0,
            "the mixing zone is all vapour, including A",
        );
    }

    #[test]
    fn mole_fractions_sum_to_one_and_specific_energy_divides_the_total() {
        let mut y = State::default();
        for component in Component::ALL {
            y.stripper.moles[component] = f64::from(component.index() as u32 + 1);
        }
        y.stripper.energy = 36.0;
        for component in Component::ALL {
            y.reactor.moles[component] = 5.0;
            y.separator.moles[component] = 5.0;
            y.mixing.moles[component] = 5.0;
        }
        y.reactor.energy = 25.0;
        y.separator.energy = 25.0;
        y.mixing.energy = 40.0;

        let unpacked = unpack(&y, TemperatureSeeds::default()).expect("converges");
        assert_exact(unpacked.stripper.total, 36.0, "1 through 8 sums to 36");
        assert!(unpacked.stripper.fractions.sums_to_one());
        assert_exact(
            unpacked.stripper.specific_energy,
            1.0,
            "36 units of energy over 36 moles",
        );
    }

    /// The seeds come out as the converged temperatures, ready for the next
    /// evaluation. Without this the warm start does not carry.
    #[test]
    fn the_converged_temperatures_become_the_next_seeds() {
        let mut y = State::default();
        for component in Component::ALL {
            y.reactor.moles[component] = 10.0;
            y.separator.moles[component] = 10.0;
            y.stripper.moles[component] = 10.0;
            y.mixing.moles[component] = 10.0;
        }
        y.reactor.energy = 60.0;
        y.separator.energy = 50.0;
        y.stripper.energy = 40.0;
        y.mixing.energy = 80.0;

        let unpacked = unpack(&y, TemperatureSeeds::default()).expect("converges");
        assert_exact(
            unpacked.seeds.reactor,
            unpacked.reactor.celsius,
            "the reactor seed is its own answer",
        );
        assert_exact(
            unpacked.seeds.mixing,
            unpacked.mixing.celsius,
            "the mixing seed is its own answer",
        );
        assert!(
            unpacked.seeds != TemperatureSeeds::default(),
            "the seeds did not move at all, so nothing was solved"
        );
    }

    /// Kelvin uses the single-precision offset, like everything else that
    /// touches `273.15` in this file.
    #[test]
    fn kelvin_uses_the_single_precision_offset() {
        let liquid = Liquid {
            moles: ByComponent::new([0.0; 8]),
            total: 0.0,
            fractions: Composition::new_unchecked([0.0; 8]),
            specific_energy: 0.0,
            celsius: 0.0,
            density: 0.0,
            volume: 0.0,
        };
        assert_exact(liquid.kelvin(), ABSOLUTE_ZERO_OFFSET, "0 C in kelvin");
        assert!(
            liquid.kelvin() < 273.15_f64,
            "the widened f32 rounds down, so this must not be the double literal"
        );
    }
}
