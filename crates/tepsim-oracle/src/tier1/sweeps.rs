//! The input pools a Tier 1 comparison sweeps over: compositions,
//! temperatures, and the sampler behind the random ones.
//!
//! Every routine Tier 1 covers takes the same two arguments, a composition and
//! a temperature, so the pools are shared across `TESUB1`, `TESUB2`, `TESUB3`
//! and `TESUB4` rather than rebuilt per routine.

use tepsim_core::{Component, Composition};

/// A deterministic sampler used only to *generate test inputs*.
///
/// This is SplitMix64. It uses integer operations exclusively, so it produces
/// an identical stream on x86-64, aarch64 and wasm32, which is what makes a
/// recorded Tier 1 number reproducible on another machine. Seeding is explicit
/// and there is no entropy source anywhere: a sweep is a pure function of its
/// [`Sweep`] parameters.
///
/// # Why not `TepRng`
///
/// The plant's own generator is itself under test. Drawing Tier 1 inputs from
/// it would couple the coverage to the thing being validated: a bug in the port
/// of `TESUB7` would silently change *which* compositions the sweep visits, so
/// the sweep could stop covering the case that would have exposed it. The
/// measuring instrument has to be independent of the thing measured.
#[derive(Clone, Debug)]
pub struct Sampler {
    state: u64,
}

impl Sampler {
    /// The odd increment SplitMix64 walks its state by, `2^64 / phi`.
    const GAMMA: u64 = 0x9E37_79B9_7F4A_7C15;

    /// Start a stream from `seed`. Every seed gives a full-period stream.
    #[must_use]
    pub const fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    /// The next 64-bit word.
    pub fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(Self::GAMMA);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// A uniform draw from `[0, 1)`.
    ///
    /// Built from the top 53 bits scaled by `2^-53`, so the multiplication is
    /// exact and every representable multiple of `2^-53` in the interval is
    /// reachable with equal probability.
    pub fn unit(&mut self) -> f64 {
        const SCALE: f64 = 1.0 / (1_u64 << 53) as f64;
        (self.next_u64() >> 11) as f64 * SCALE
    }

    /// A uniform draw from `low ..= high`.
    ///
    /// Deliberately two roundings rather than a fused multiply-add. What this
    /// value has to be is *reproducible on every target*, not accurate to the
    /// last bit, and a plain multiply and add is the same everywhere by the
    /// IEEE-754 definition alone. `mul_add` gets there by a different route on
    /// hardware with an FMA instruction than on hardware without, which is more
    /// trust than a test input needs to be worth.
    #[allow(clippy::suboptimal_flops, reason = "determinism over accuracy")]
    pub fn range(&mut self, low: f64, high: f64) -> f64 {
        low + (high - low) * self.unit()
    }

    /// An integer uniform on `0 .. n`, by Lemire's multiply-shift.
    ///
    /// Slightly biased for `n` that does not divide `2^64`, by at most `n /
    /// 2^64`. Rejection sampling would remove it, at the cost of making the
    /// number of words consumed per draw depend on the values drawn, which is a
    /// worse property for a harness whose whole job is reproducibility.
    fn below(&mut self, n: usize) -> usize {
        ((u128::from(self.next_u64()) * n as u128) >> 64) as usize
    }

    /// `parts` non-negative values summing to one, uniform on the
    /// `(parts - 1)`-simplex, returned in the first `parts` slots.
    ///
    /// This is the spacings method: `parts - 1` uniforms are sorted, and the
    /// gaps between `0`, them, and `1` are the values. Those gaps are exactly
    /// `Dirichlet(1, ..., 1)` distributed, and the construction has three
    /// properties a gamma-ratio sampler would not: it never rejects, so the
    /// word count per sample is fixed; it uses no transcendental function, so
    /// it cannot differ by a ULP between platforms; and every value is
    /// non-negative by construction rather than by luck.
    fn spacings(&mut self, parts: usize) -> [f64; Component::COUNT] {
        debug_assert!((1..=Component::COUNT).contains(&parts));

        let mut buffer = [0.0_f64; Component::COUNT - 1];
        let cuts = &mut buffer[..parts - 1];
        for cut in cuts.iter_mut() {
            *cut = self.unit();
        }
        cuts.sort_by(f64::total_cmp);

        let mut values = [0.0_f64; Component::COUNT];
        let mut previous = 0.0;
        for (value, &cut) in values.iter_mut().zip(cuts.iter()) {
            *value = cut - previous;
            previous = cut;
        }
        values[parts - 1] = 1.0 - previous;
        values
    }

    /// A uniform sample from the interior of the 8-component simplex.
    pub fn dirichlet(&mut self) -> Composition {
        Composition::new(self.spacings(Component::COUNT))
    }

    /// A uniform sample from a random *face* of the simplex: a support of
    /// between one and eight species, uniform within that support and exactly
    /// zero outside it.
    ///
    /// The interesting numerics live near the boundary. Real TEP streams carry
    /// several mole fractions at 1e-5 and below, and a sum of products is at
    /// its least accurate when the terms differ by many orders of magnitude.
    /// The simplex grid reaches the boundary only at its own coarse spacing;
    /// this pool reaches it at full resolution, and reaches the vertices and
    /// edges exactly.
    pub fn simplex_face(&mut self) -> Composition {
        let support = 1 + self.below(Component::COUNT);

        // Partial Fisher-Yates: the first `support` entries become a uniformly
        // chosen subset of the eight species, in a uniformly random order.
        let mut species = [0_usize; Component::COUNT];
        for (i, slot) in species.iter_mut().enumerate() {
            *slot = i;
        }
        for i in 0..support {
            species.swap(i, i + self.below(Component::COUNT - i));
        }

        let values = self.spacings(support);
        let mut fractions = [0.0_f64; Component::COUNT];
        for (&index, &value) in species[..support].iter().zip(values.iter()) {
            fractions[index] = value;
        }
        Composition::new(fractions)
    }
}

/// Every composition on a regular grid over the 8-component simplex.
///
/// A grid point is `counts[i] / divisions` for non-negative integers summing to
/// `divisions`, so there are `C(divisions + 7, 7)` of them: 6435 at 8
/// divisions, 245157 at 16.
///
/// # Use a power of two
///
/// With `divisions` a power of two every fraction is a dyadic rational, so each
/// one is exact in `f64` and each grid point sums to exactly 1.0. At
/// `divisions = 10` the points are off by a rounding error instead, which is
/// harmless for the model but makes any later bit-exactness claim about the
/// *inputs* untrue. There is no cost to the power of two, so take it.
#[derive(Clone, Debug)]
pub struct SimplexGrid {
    divisions: u32,
    counts: [u32; Component::COUNT],
    done: bool,
}

impl SimplexGrid {
    /// Start a grid with the given number of divisions per edge.
    ///
    /// # Panics
    ///
    /// If `divisions` is zero, which would describe no simplex at all.
    #[must_use]
    pub fn new(divisions: u32) -> Self {
        assert!(divisions > 0, "a simplex grid needs at least one division");
        let mut counts = [0; Component::COUNT];
        counts[0] = divisions;
        Self {
            divisions,
            counts,
            done: false,
        }
    }

    /// How many points the grid has: `C(divisions + 7, 7)`.
    #[must_use]
    pub fn len(divisions: u32) -> usize {
        let n = u128::from(divisions);
        let mut count: u128 = 1;
        for k in 1..Component::COUNT as u128 {
            count = count * (n + k) / k;
        }
        count as usize
    }

    /// Step the odometer to the next composition of `divisions` into eight
    /// parts, in the standard colexicographic order.
    fn advance(&mut self) {
        const LAST: usize = Component::COUNT - 1;
        if self.counts[LAST] == self.divisions {
            self.done = true;
            return;
        }

        let mut i = LAST;
        while self.counts[i] == 0 {
            i -= 1;
        }

        if i == LAST {
            // The whole remainder has reached the end; carry it back to the
            // next donor on the left and dump it one place further right.
            let mut j = LAST - 1;
            while self.counts[j] == 0 {
                debug_assert!(j > 0, "a non-zero entry must exist to the left");
                j -= 1;
            }
            let carried = self.counts[LAST];
            self.counts[LAST] = 0;
            self.counts[j] -= 1;
            self.counts[j + 1] = carried + 1;
        } else {
            self.counts[i] -= 1;
            self.counts[i + 1] += 1;
        }
    }
}

impl Iterator for SimplexGrid {
    type Item = Composition;

    fn next(&mut self) -> Option<Composition> {
        if self.done {
            return None;
        }
        let counts = self.counts;
        self.advance();

        let divisions = f64::from(self.divisions);
        let mut fractions = [0.0_f64; Component::COUNT];
        for (fraction, count) in fractions.iter_mut().zip(counts) {
            *fraction = f64::from(count) / divisions;
        }
        Some(Composition::new(fractions))
    }
}

/// A closed temperature interval in degrees Celsius.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TemperatureRange {
    /// The low end, included.
    pub low: f64,
    /// The high end, included.
    pub high: f64,
}

impl TemperatureRange {
    /// Everything the plant can reach before it trips, in degrees Celsius.
    ///
    /// The ceiling is the reactor-temperature shutdown limit: `teprob.f:706`
    /// sets `ISD` when `XMEAS(9)`, which is `TCR`, exceeds 175. No stream can
    /// sit above that for more than one step of a run that has not already
    /// ended. The floor is 0 rather than the coldest thing the plant actually
    /// contains, which is the 35 degree cooling water supply at
    /// `teprob.f:1323`, because the stripper temperature factor at
    /// `teprob.f:617` has an explicit branch for `TCC` below 5.292 and a sweep
    /// that never went there would never exercise it.
    ///
    /// # This is the same band for all three `ITY` modes
    ///
    /// Not an approximation. Every `TESUB1` call site passes some `TST(i)`
    /// (`teprob.f:555-564`, `654-655`), and every `TST(i)` is either the 45
    /// degree feed temperature or one of `TCV`, `TCR`, `TCS`, `TCC`
    /// (`teprob.f:411-412`, `549-554`, `652-653`). The vapour correlations are
    /// evaluated at exactly the temperatures the liquid ones are, so splitting
    /// the domain by mode would describe a distinction the source does not
    /// make.
    pub const PLANT: Self = Self {
        low: 0.0,
        high: 175.0,
    };

    /// Whether `celsius` lies within this range.
    #[must_use]
    pub fn contains(&self, celsius: f64) -> bool {
        celsius >= self.low && celsius <= self.high
    }
}

/// A temperature the model treats specially, and the line that makes it so.
///
/// A uniform grid steps over a branch boundary with probability one, so the
/// exact values go in the sweep by hand. This is the Tier 1 analogue of the
/// adversarial state pool `PLAN.org` calls for at Tier 2.
#[derive(Clone, Copy, Debug)]
pub struct Breakpoint {
    /// Degrees Celsius.
    pub celsius: f64,
    /// What changes here, and where in `teprob.f` to read it.
    pub why: &'static str,
}

/// Every temperature in [`TemperatureRange::PLANT`] that the source singles
/// out, in ascending order.
pub const BREAKPOINTS: &[Breakpoint] = &[
    Breakpoint {
        celsius: 0.0,
        why: "TESUB1 enthalpy is identically zero here for ITY 0 and 1, since \
              every term carries a factor of T (teprob.f:1395, 1402)",
    },
    Breakpoint {
        celsius: 5.292,
        why: "below this the stripper temperature factor is pinned at 0.1 \
              instead of following the hyperbola (teprob.f:617-618)",
    },
    Breakpoint {
        celsius: 35.0,
        why: "nominal reactor cooling water supply temperature, SZERO(5) \
              (teprob.f:1323)",
    },
    Breakpoint {
        celsius: 40.0,
        why: "nominal condenser cooling water supply temperature, SZERO(6) \
              (teprob.f:1328)",
    },
    Breakpoint {
        celsius: 45.0,
        why: "the four feed stream temperatures, TST(1..4) (teprob.f:1142, \
              1151, 1160, 1169)",
    },
    Breakpoint {
        celsius: 100.0,
        why: "above this the condenser duty QUC is exactly zero rather than \
              proportional to the driving force (teprob.f:678)",
    },
    Breakpoint {
        celsius: 170.0,
        why: "above this the stripper temperature factor becomes the linear \
              branch (teprob.f:615-616)",
    },
    Breakpoint {
        celsius: 175.0,
        why: "the reactor temperature shutdown limit (teprob.f:706)",
    },
];

/// The parameters of one Tier 1 sweep.
///
/// A `Sweep` is a pure description: the same value always enumerates the same
/// cases in the same order, on any platform. Record it alongside the numbers it
/// produced and the numbers can be reproduced.
#[derive(Clone, Copy, Debug)]
pub struct Sweep {
    /// Divisions per edge of the simplex grid. Take a power of two; see
    /// [`SimplexGrid`].
    pub grid_divisions: u32,
    /// How many uniformly spaced temperatures to cross the grid with. The
    /// [`BREAKPOINTS`] are added to these.
    pub temperature_steps: usize,
    /// How many random interior compositions to draw.
    pub dirichlet_samples: usize,
    /// How many random boundary compositions to draw.
    pub face_samples: usize,
    /// The temperature domain.
    pub range: TemperatureRange,
    /// Seed for the two random pools.
    pub seed: u64,
}

impl Sweep {
    /// Small enough to run in a fraction of a second. For testing the harness
    /// itself and for a quick check while iterating.
    pub const SMOKE: Self = Self {
        grid_divisions: 4,
        temperature_steps: 9,
        dirichlet_samples: 2_000,
        face_samples: 1_000,
        range: TemperatureRange::PLANT,
        seed: 0x050E_15EE,
    };

    /// The volume `PLAN.org` asks for at Tier 1: about ten million cases.
    pub const FULL: Self = Self {
        grid_divisions: 8,
        temperature_steps: 48,
        dirichlet_samples: 8_000_000,
        face_samples: 1_640_000,
        range: TemperatureRange::PLANT,
        seed: 0x07E9_5EED,
    };

    /// The temperatures the grid pool is crossed with: a uniform ladder over
    /// [`Sweep::range`] plus every [`BREAKPOINTS`] entry inside it, ascending,
    /// with duplicates removed.
    #[must_use]
    #[allow(clippy::suboptimal_flops, reason = "see Sampler::range")]
    pub fn temperatures(&self) -> Vec<f64> {
        let mut temperatures = Vec::with_capacity(self.temperature_steps + BREAKPOINTS.len());
        let last = self.temperature_steps.saturating_sub(1);
        for step in 0..self.temperature_steps {
            // The endpoint is assigned rather than computed so that it is the
            // range bound exactly, not the bound plus a rounding error.
            temperatures.push(if step == last {
                self.range.high
            } else {
                let fraction = step as f64 / last as f64;
                self.range.low + (self.range.high - self.range.low) * fraction
            });
        }
        for breakpoint in BREAKPOINTS {
            if self.range.contains(breakpoint.celsius) {
                temperatures.push(breakpoint.celsius);
            }
        }
        temperatures.sort_by(f64::total_cmp);
        temperatures.dedup_by(|a, b| a.to_bits() == b.to_bits());
        temperatures
    }

    /// How many cases [`Sweep::cases`] will yield.
    #[must_use]
    pub fn len(&self) -> usize {
        SimplexGrid::len(self.grid_divisions) * self.temperatures().len()
            + self.dirichlet_samples
            + self.face_samples
    }

    /// Whether the sweep yields nothing at all.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Every case, grid pool first, then the interior samples, then the
    /// boundary samples.
    ///
    /// The grid pool is temperature-major: the whole composition grid at the
    /// first temperature, then the whole grid at the second. Each random case
    /// draws its composition before its temperature. Both orderings are part of
    /// what a seed means, so changing either changes every recorded number.
    pub fn cases(self) -> impl Iterator<Item = Case> {
        let divisions = self.grid_divisions;
        let grid = self
            .temperatures()
            .into_iter()
            .flat_map(move |celsius| {
                SimplexGrid::new(divisions).map(move |composition| (composition, celsius))
            })
            .enumerate()
            .map(|(index, (composition, celsius))| Case {
                pool: Pool::Grid,
                index,
                composition,
                celsius,
            });

        let range = self.range;
        let mut interior = Sampler::new(self.seed);
        let dirichlet = (0..self.dirichlet_samples).map(move |index| {
            let composition = interior.dirichlet();
            Case {
                pool: Pool::Dirichlet,
                index,
                composition,
                celsius: interior.range(range.low, range.high),
            }
        });

        // A separate stream, so that changing one pool's size leaves the other
        // pool's cases untouched and its recorded numbers comparable.
        let mut boundary = Sampler::new(self.seed ^ u64::MAX);
        let faces = (0..self.face_samples).map(move |index| {
            let composition = boundary.simplex_face();
            Case {
                pool: Pool::Face,
                index,
                composition,
                celsius: boundary.range(range.low, range.high),
            }
        });

        grid.chain(dirichlet).chain(faces)
    }
}

/// Which pool a case came from.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Pool {
    /// The regular simplex grid crossed with the temperature ladder.
    Grid,
    /// Uniform random interior compositions.
    Dirichlet,
    /// Uniform random compositions on a random face of the simplex.
    Face,
}

impl core::fmt::Display for Pool {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.write_str(match self {
            Pool::Grid => "grid",
            Pool::Dirichlet => "dirichlet",
            Pool::Face => "face",
        })
    }
}

/// One evaluation point.
#[derive(Clone, Copy, Debug)]
pub struct Case {
    /// Which pool produced it.
    pub pool: Pool,
    /// Its position within that pool, zero based.
    pub index: usize,
    /// The composition.
    pub composition: Composition,
    /// The temperature in degrees Celsius.
    pub celsius: f64,
}

impl Case {
    /// The composition as the raw eight-element array the Fortran expects.
    #[must_use]
    pub fn z(&self) -> [f64; Component::COUNT] {
        *self.composition.fractions().as_array()
    }
}

impl core::fmt::Display for Case {
    /// Terse enough for a one-line report, complete enough to find the case
    /// again: the temperature is printed with Rust's shortest round-tripping
    /// form, so re-reading it recovers the same bits.
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "{}#{} T={}", self.pool, self.index, self.celsius)
    }
}
