//! What a run costs, per hour, from what the plant actually measures.
//!
//! # Measurements, not state
//!
//! Every input here is an `XMEAS` channel. That is deliberate: an operations
//! layer on a real plant has instruments and not a state vector, and a cost
//! function written against internal state would be unusable on the published
//! `d00`-`d21` files, which are measurements. It also means everything here
//! sees the same noise, dead time and analyser staircases an operator sees.
//!
//! # The four terms, and the `teprob.f` line that defines each unit
//!
//! | Term | Channel | Definition | Unit |
//! |---|---|---|---|
//! | Purge loss | `XMEAS(10)` with `XMEAS(29..36)` | `teprob.f:688` | kscmh, mol % |
//! | Product loss | `XMEAS(17)` with `XMEAS(37..41)` | `teprob.f:695` | m³/h, mol % |
//! | Steam | `XMEAS(19)` | `teprob.f:697` | kg/h |
//! | Compressor | `XMEAS(20)` | `teprob.f:699` | kW |
//!
//! Material leaving in the purge and in the product is what the process
//! actually loses money on: the purge vents unreacted A, C and inert B, and
//! product carrying unconverted D and E is reactant that was paid for and not
//! sold.
//!
//! # Converting a flow to a mass rate, which is where the care is needed
//!
//! `FTM` is molar, in lbmol/h, and `teprob.f:688` is the proof rather than an
//! assumption: `XMEAS(10)=FTM(10)*0.359/35.3145` is 0.359 thousand standard
//! cubic feet per lbmol, then cubic feet to cubic metres, giving kscmh. So the
//! purge molar rate is recoverable by inverting that, and a per-component mass
//! rate follows from the analyser and `XMW`.
//!
//! The product is harder and more interesting. `teprob.f:695` is
//! `XMEAS(17)=FTM(13)/DLC/35.3145`, a *volumetric* flow, and the molar density
//! `DLC` is not measured. It is computable: the stripper underflow composition
//! is `XMEAS(37..41)` and its temperature is `XMEAS(18)`, and
//! [`tepsim_core::thermo::liquid_density`] is the same correlation the plant
//! itself used, already validated against the oracle. Closing an unmeasured
//! quantity with the model is exactly what this layer is for.

use tepsim_core::constants::XMW;
use tepsim_core::{Component, Composition, thermo::liquid_density};

/// Thousand standard cubic feet per lbmol, from `teprob.f:683-688`.
///
/// The same literal appears on every gas flow measurement in that block, which
/// is what makes it a unit conversion rather than a per-stream coefficient.
const KSCF_PER_LBMOL: f64 = 0.359;

/// Cubic feet per cubic metre, from the same lines.
const CUBIC_FEET_PER_CUBIC_METRE: f64 = 35.3145;

/// Kilograms per pound.
///
/// The exact international definition, deliberately *not* the `0.454` that
/// appears at `teprob.f:697`. That constant is part of the model, is inside the
/// measurement this crate reads, and reproducing it is the simulator's job.
/// This one converts a mass this crate computed, so it should be right rather
/// than bit-compatible with a 1993 rounding.
const KG_PER_LB: f64 = 0.453_592_37;

/// Prices, per hour of operation.
///
/// # These are inputs, and there is no default
///
/// Downs and Vogel state an operating cost, and its coefficients are in the
/// paper. They are **not** in `reference/`, so nothing in this repository can
/// assert them, and this project's rule is that a constant which cannot be
/// checked against its source is a silent failure waiting to happen. A number
/// typed here off a remembered table would be exactly that.
///
/// It is also the wrong shape for the problem. The use case this crate exists
/// for is planning against prices that *move*, so a price is a time-varying
/// input and not a constant of the process.
///
/// So: supply your own. [`Prices::unit`] exists for tests and for reading the
/// terms in whatever units you put in.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Prices {
    /// Currency per kilogram of each component vented in the purge.
    ///
    /// Indexed by [`Component`]. Inert B has a disposal cost rather than a
    /// value; A and C are unreacted feed.
    pub purge: [f64; Component::COUNT],
    /// Currency per kilogram of each component leaving in the product stream.
    ///
    /// G and H are the products and would normally be negative here, since
    /// selling them is income. D and E leaving unconverted is a loss.
    pub product: [f64; Component::COUNT],
    /// Currency per kilowatt-hour of compressor work.
    pub compressor: f64,
    /// Currency per kilogram of stripper steam.
    pub steam: f64,
}

impl Prices {
    /// Every price 1.0, which turns a [`CostRate`] into a report of the
    /// underlying physical rates: kg/h, kWh/h and kg/h respectively.
    ///
    /// Useful for reading what the plant is actually doing before deciding what
    /// it is worth, and it is what the tests here use, since the structure is
    /// the part this crate can validate.
    #[must_use]
    pub const fn unit() -> Self {
        Self {
            purge: [1.0; Component::COUNT],
            product: [1.0; Component::COUNT],
            compressor: 1.0,
            steam: 1.0,
        }
    }

    /// Every price zero: a run costs nothing whatever it does.
    #[must_use]
    pub const fn free() -> Self {
        Self {
            purge: [0.0; Component::COUNT],
            product: [0.0; Component::COUNT],
            compressor: 0.0,
            steam: 0.0,
        }
    }
}

/// A cost rate, broken into the four terms so a total can be argued with.
///
/// A single number tells you what a run costs and nothing about why, which is
/// useless to a planner. These are per hour, in whatever currency [`Prices`]
/// was denominated in.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CostRate {
    /// Material vented in the purge, stream 9.
    pub purge: f64,
    /// Material leaving in the product, stream 11.
    pub product: f64,
    /// Compressor work.
    pub compressor: f64,
    /// Stripper steam.
    pub steam: f64,
}

impl CostRate {
    /// The four terms added up.
    #[must_use]
    pub fn total(&self) -> f64 {
        self.purge + self.product + self.compressor + self.steam
    }

    /// The terms, with names, largest first.
    ///
    /// For saying *why* a plant is expensive, which is the question a planner
    /// is actually asked.
    #[must_use]
    pub fn ranked(&self) -> [(&'static str, f64); 4] {
        let mut terms = [
            ("purge", self.purge),
            ("product", self.product),
            ("compressor", self.compressor),
            ("steam", self.steam),
        ];
        terms.sort_by(|a, b| b.1.total_cmp(&a.1));
        terms
    }
}

/// Molar flow of the purge, lbmol/h, by inverting `teprob.f:688`.
#[must_use]
pub fn purge_molar_rate(purge_kscmh: f64) -> f64 {
    purge_kscmh * CUBIC_FEET_PER_CUBIC_METRE / KSCF_PER_LBMOL
}

/// Molar flow of the product, lbmol/h, by inverting `teprob.f:695`.
///
/// `DLC` is not measured, so it is computed from the product analyser and the
/// stripper temperature with the plant's own liquid density correlation.
#[must_use]
pub fn product_molar_rate(
    product_m3_per_hour: f64,
    composition: &Composition,
    celsius: f64,
) -> f64 {
    product_m3_per_hour * CUBIC_FEET_PER_CUBIC_METRE * liquid_density(composition, celsius)
}

/// The purge composition, from `XMEAS(29..36)`, normalised.
///
/// # Why normalising is correct rather than convenient
///
/// The analysers are noisy and read in mole percent, so eight readings do not
/// sum to exactly 100 and [`Composition::new`] checks the sum in debug builds.
/// Renormalising is what an operations layer does with a noisy analyser, and
/// the alternative, trusting the sum, is what puts a 0.3% flow error into every
/// downstream mass balance.
#[must_use]
pub fn purge_composition(measurements: &[f64]) -> Option<Composition> {
    let mut fractions = [0.0; Component::COUNT];
    for (i, slot) in fractions.iter_mut().enumerate() {
        *slot = *measurements.get(28 + i)?;
    }
    normalise(fractions)
}

/// The product composition, from `XMEAS(37..41)`, normalised.
///
/// The product analyser reports D through H only, because A, B and C are not
/// present in the stripper underflow in any quantity the original models. The
/// three leading fractions are therefore zero rather than unknown.
#[must_use]
pub fn product_composition(measurements: &[f64]) -> Option<Composition> {
    let mut fractions = [0.0; Component::COUNT];
    for (i, slot) in fractions.iter_mut().enumerate().skip(3) {
        *slot = *measurements.get(36 + i - 3)?;
    }
    normalise(fractions)
}

/// Scale to sum to one, or `None` if there is nothing to scale.
fn normalise(mut fractions: [f64; Component::COUNT]) -> Option<Composition> {
    let sum: f64 = fractions.iter().sum();
    if !sum.is_finite() || sum <= 0.0 {
        return None;
    }
    for slot in &mut fractions {
        *slot /= sum;
    }
    Some(Composition::new(fractions))
}

/// Mass rate of each component in a stream, kg/h.
fn component_mass_rates(molar_rate: f64, composition: &Composition) -> [f64; Component::COUNT] {
    let mut rates = [0.0; Component::COUNT];
    for component in Component::ALL {
        rates[component as usize] =
            molar_rate * composition[component] * XMW[component] * KG_PER_LB;
    }
    rates
}

/// The cost rate at one instant, from one sample's measurements.
///
/// `measurements` is `XMEAS(1..41)` zero-indexed, which is what
/// [`tepsim::Sample::measurements`] gives.
///
/// # Returns
///
/// `None` if the slice is short, or if an analyser reads all zeros so no
/// composition can be formed. Both are conditions of the *input*, not failures
/// of the calculation, and a planner should be able to tell them apart from an
/// expensive plant.
#[must_use]
pub fn cost_rate(measurements: &[f64], prices: &Prices) -> Option<CostRate> {
    let purge_kscmh = *measurements.get(9)?;
    let product_m3 = *measurements.get(16)?;
    let stripper_celsius = *measurements.get(17)?;
    let steam_kg_per_hour = *measurements.get(18)?;
    let compressor_kw = *measurements.get(19)?;

    let purge_x = purge_composition(measurements)?;
    let product_x = product_composition(measurements)?;

    let purge_rates = component_mass_rates(purge_molar_rate(purge_kscmh), &purge_x);
    let product_rates = component_mass_rates(
        product_molar_rate(product_m3, &product_x, stripper_celsius),
        &product_x,
    );

    // `mul_add`, where `tepsim-core` forbids it.
    //
    // The lint that suggests this is allowed crate-wide in `tepsim-core`,
    // because a fused multiply-add rounds once where the Fortran rounds twice
    // and the whole point of that crate is to round the way the Fortran does.
    // None of that applies here. This crate is not bit-matching anything, a
    // single rounding is simply more accurate over an eight-term sum, and
    // `mul_add` is specified as fused rather than left to the target, so it
    // stays deterministic across platforms.
    let mut purge = 0.0_f64;
    let mut product = 0.0_f64;
    for component in Component::ALL {
        let i = component as usize;
        purge = purge_rates[i].mul_add(prices.purge[i], purge);
        product = product_rates[i].mul_add(prices.product[i], product);
    }

    Some(CostRate {
        purge,
        product,
        // kW is already energy per hour, so a per-kWh price needs no conversion.
        compressor: compressor_kw * prices.compressor,
        steam: steam_kg_per_hour * prices.steam,
    })
}

/// The mean cost rate over a run, and the total it integrates to.
///
/// Samples are evenly spaced, so the mean of the rates is the time average and
/// no trapezoid is needed. Samples whose analysers have not reported yet are
/// skipped rather than counted as zero, which would drag the mean down for the
/// first dead time of every run.
#[must_use]
pub fn mean_cost_rate(run: &tepsim::Run, prices: &Prices) -> Option<(CostRate, usize)> {
    let mut sum = CostRate {
        purge: 0.0,
        product: 0.0,
        compressor: 0.0,
        steam: 0.0,
    };
    let mut counted = 0_usize;
    for sample in &run.samples {
        if let Some(rate) = cost_rate(&sample.measurements, prices) {
            sum.purge += rate.purge;
            sum.product += rate.product;
            sum.compressor += rate.compressor;
            sum.steam += rate.steam;
            counted += 1;
        }
    }
    if counted == 0 {
        return None;
    }
    let n = counted as f64;
    Some((
        CostRate {
            purge: sum.purge / n,
            product: sum.product / n,
            compressor: sum.compressor / n,
            steam: sum.steam / n,
        },
        counted,
    ))
}
