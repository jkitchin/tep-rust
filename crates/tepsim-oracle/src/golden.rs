//! The golden trace: a recorded run of the original Fortran, committed to the
//! repository so the fidelity preflight works without a Fortran toolchain.
//!
//! # Why raw bit patterns
//!
//! Every value is stored as the hexadecimal IEEE-754 bit pattern of an `f64`.
//! Decimal text does not round-trip exactly at every precision, and a
//! comparison that is only approximately exact is not much use as the anchor
//! for a validation ladder. Hex bits also make a diff show exactly which value
//! moved.
//!
//! # Why 100 steps
//!
//! Long enough to leave the first-step special cases behind, short enough that
//! the preflight runs in well under a second on any machine. The first step is
//! the interesting one anyway: `TEFUNC` draws nothing at t=0 and exactly 264
//! uniforms on the step after, so a port that gets the ordering wrong diverges
//! immediately rather than subtly.
//!
//! # Format
//!
//! Comment lines start with `#` and carry provenance as `key: value`. Data
//! lines are the step index followed by [`VALUES_PER_STEP`] hex words: 50
//! states as they were *before* the step, 50 derivatives, 41 measurements, and
//! the generator word, in that order.

use core::fmt::Write as _;

/// States, then derivatives, then measurements, then the generator word.
pub const VALUES_PER_STEP: usize = 50 + 50 + 41 + 1;

/// Steps recorded. See the module docs for why this number.
pub const STEPS: usize = 100;

/// One second, in hours: the step size the original's `INTGTR` uses.
pub const DT_HOURS: f64 = 1.0 / 3600.0;

/// The seed compiled into `teprob.f:1187`.
pub const SEED: f64 = 4651207995.0;

/// Where the committed trace lives, relative to the workspace root.
///
/// Deliberately not under `reference/`: that holds vendored upstream material
/// that never changes, whereas this is our artefact and is regenerated on a
/// deliberate re-baseline.
pub const PATH: &str = "golden/nominal-100-steps.trace";

/// One recorded step.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Step {
    /// The 50 states as they were before this step was evaluated.
    pub states: [f64; 50],
    /// The 50 derivatives `TEFUNC` returned.
    pub derivatives: [f64; 50],
    /// `XMEAS(1..41)` after the call, including measurement noise.
    pub measurements: [f64; 41],
    /// `COMMON/RANDSD/ G` after the call.
    pub rng: f64,
}

/// A parsed trace, with the provenance recorded when it was generated.
#[derive(Clone, Debug, PartialEq)]
pub struct Trace {
    /// gfortran version string, as `major.minor.patch`.
    pub gfortran: String,
    /// The exact flag set passed to gfortran.
    pub fflags: String,
    /// The generator seed the run started from.
    pub seed: f64,
    /// Integration step in hours.
    pub dt_hours: f64,
    /// The recorded steps, in order.
    pub steps: Vec<Step>,
}

/// Something wrong with a trace file.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum TraceError {
    /// A required provenance header was absent.
    MissingHeader(&'static str),
    /// A data line did not have `1 + VALUES_PER_STEP` fields.
    WrongFieldCount {
        /// One-based line number in the file.
        line: usize,
        /// How many fields were present.
        got: usize,
    },
    /// A field was not a 16-digit hexadecimal word.
    BadHex {
        /// One-based line number in the file.
        line: usize,
        /// The offending text.
        field: String,
    },
    /// Step indices must run 0, 1, 2, ... with no gaps.
    StepOutOfOrder {
        /// What the index should have been.
        expected: usize,
        /// What it was.
        got: usize,
    },
    /// The file recorded a different number of steps than expected.
    WrongStepCount {
        /// How many steps were expected.
        expected: usize,
        /// How many were present.
        got: usize,
    },
    /// A recorded value was not finite, which no valid run produces.
    NotFinite {
        /// Zero-based step index.
        step: usize,
        /// Zero-based position within the step's values.
        index: usize,
    },
}

impl core::fmt::Display for TraceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::MissingHeader(key) => write!(f, "missing provenance header `{key}`"),
            Self::WrongFieldCount { line, got } => write!(
                f,
                "line {line}: expected {} fields, found {got}",
                VALUES_PER_STEP + 1
            ),
            Self::BadHex { line, field } => {
                write!(f, "line {line}: `{field}` is not a 16-digit hex word")
            }
            Self::StepOutOfOrder { expected, got } => {
                write!(
                    f,
                    "step indices must be consecutive: expected {expected}, found {got}"
                )
            }
            Self::WrongStepCount { expected, got } => {
                write!(f, "expected {expected} steps, found {got}")
            }
            Self::NotFinite { step, index } => write!(
                f,
                "step {step}, value {index} is not finite; no valid run produces that"
            ),
        }
    }
}

impl std::error::Error for TraceError {}

impl Trace {
    /// Render to the on-disk format.
    pub fn to_text(&self) -> String {
        let mut out = String::with_capacity(STEPS * VALUES_PER_STEP * 17 + 1024);
        out.push_str("# tep-rust golden trace, format 1\n");
        out.push_str("#\n");
        out.push_str("# A recorded run of the ORIGINAL Fortran, committed so the fidelity\n");
        out.push_str("# preflight works without a Fortran toolchain. Regenerate with:\n");
        out.push_str("#   cargo run -p tepsim-oracle --features oracle --bin gen-golden-trace\n");
        out.push_str("#\n");
        out.push_str("# Regenerating is a deliberate re-baseline: every validation number in\n");
        out.push_str("# LOG.org was measured against these values. Do not regenerate to make a\n");
        out.push_str("# failing test pass.\n");
        out.push_str("#\n");
        let _ = writeln!(out, "# gfortran: {}", self.gfortran);
        let _ = writeln!(out, "# fflags: {}", self.fflags);
        let _ = writeln!(out, "# seed: {:016x}", self.seed.to_bits());
        let _ = writeln!(out, "# dt_hours: {:016x}", self.dt_hours.to_bits());
        let _ = writeln!(out, "# steps: {}", self.steps.len());
        let _ = writeln!(out, "# values_per_step: {VALUES_PER_STEP}");
        out.push_str("# layout: step, 50 states (before), 50 derivatives, 41 measurements, rng\n");
        out.push_str("# values are hexadecimal IEEE-754 f64 bit patterns\n");

        for (i, step) in self.steps.iter().enumerate() {
            let _ = write!(out, "{i}");
            for value in step
                .states
                .iter()
                .chain(&step.derivatives)
                .chain(&step.measurements)
                .chain(core::iter::once(&step.rng))
            {
                let _ = write!(out, " {:016x}", value.to_bits());
            }
            out.push('\n');
        }
        out
    }

    /// Parse the on-disk format, checking structure as it goes.
    pub fn parse(text: &str) -> Result<Self, TraceError> {
        let mut gfortran = None;
        let mut fflags = None;
        let mut seed = None;
        let mut dt_hours = None;
        let mut steps = Vec::new();

        for (line_no, raw) in text.lines().enumerate() {
            let line_no = line_no + 1;
            let line = raw.trim();
            if line.is_empty() {
                continue;
            }
            if let Some(comment) = line.strip_prefix('#') {
                let comment = comment.trim();
                if let Some((key, value)) = comment.split_once(':') {
                    let value = value.trim();
                    match key.trim() {
                        "gfortran" => gfortran = Some(value.to_string()),
                        "fflags" => fflags = Some(value.to_string()),
                        "seed" => seed = parse_hex(value, line_no).ok(),
                        "dt_hours" => dt_hours = parse_hex(value, line_no).ok(),
                        _ => {}
                    }
                }
                continue;
            }

            let fields: Vec<&str> = line.split_whitespace().collect();
            if fields.len() != VALUES_PER_STEP + 1 {
                return Err(TraceError::WrongFieldCount {
                    line: line_no,
                    got: fields.len(),
                });
            }
            let index: usize = fields[0].parse().map_err(|_| TraceError::BadHex {
                line: line_no,
                field: fields[0].to_string(),
            })?;
            if index != steps.len() {
                return Err(TraceError::StepOutOfOrder {
                    expected: steps.len(),
                    got: index,
                });
            }

            let mut values = [0.0_f64; VALUES_PER_STEP];
            for (slot, field) in values.iter_mut().zip(&fields[1..]) {
                *slot = parse_hex(field, line_no)?;
            }
            for (i, v) in values.iter().enumerate() {
                if !v.is_finite() {
                    return Err(TraceError::NotFinite {
                        step: index,
                        index: i,
                    });
                }
            }

            let mut step = Step {
                states: [0.0; 50],
                derivatives: [0.0; 50],
                measurements: [0.0; 41],
                rng: values[VALUES_PER_STEP - 1],
            };
            step.states.copy_from_slice(&values[0..50]);
            step.derivatives.copy_from_slice(&values[50..100]);
            step.measurements.copy_from_slice(&values[100..141]);
            steps.push(step);
        }

        Ok(Self {
            gfortran: gfortran.ok_or(TraceError::MissingHeader("gfortran"))?,
            fflags: fflags.ok_or(TraceError::MissingHeader("fflags"))?,
            seed: seed.ok_or(TraceError::MissingHeader("seed"))?,
            dt_hours: dt_hours.ok_or(TraceError::MissingHeader("dt_hours"))?,
            steps,
        })
    }

    /// Fail unless the trace holds exactly [`STEPS`] steps.
    pub fn require_full_length(&self) -> Result<(), TraceError> {
        if self.steps.len() == STEPS {
            Ok(())
        } else {
            Err(TraceError::WrongStepCount {
                expected: STEPS,
                got: self.steps.len(),
            })
        }
    }
}

fn parse_hex(field: &str, line: usize) -> Result<f64, TraceError> {
    if field.len() != 16 {
        return Err(TraceError::BadHex {
            line,
            field: field.to_string(),
        });
    }
    u64::from_str_radix(field, 16)
        .map(f64::from_bits)
        .map_err(|_| TraceError::BadHex {
            line,
            field: field.to_string(),
        })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample() -> Trace {
        Trace {
            gfortran: "15.2.0".into(),
            fflags: "-c -O0".into(),
            seed: SEED,
            dt_hours: DT_HOURS,
            steps: vec![
                Step {
                    states: [1.5; 50],
                    derivatives: [-0.25; 50],
                    measurements: [120.4; 41],
                    rng: 12345.0,
                },
                Step {
                    states: [1.75; 50],
                    derivatives: [0.5; 50],
                    measurements: [120.5; 41],
                    rng: 67890.0,
                },
            ],
        }
    }

    #[test]
    fn round_trips_through_the_text_format() {
        let trace = sample();
        assert_eq!(Trace::parse(&trace.to_text()).expect("parse"), trace);
    }

    /// The reason the format stores bits rather than decimals.
    #[test]
    fn awkward_values_survive_the_round_trip_exactly() {
        let mut trace = sample();
        trace.steps[0].states[0] = 1.0 / 3.0;
        trace.steps[0].states[1] = f64::MIN_POSITIVE;
        trace.steps[0].states[2] = -0.0;
        trace.steps[0].derivatives[0] = 1.234_567_890_123_456_7e-300;
        let parsed = Trace::parse(&trace.to_text()).expect("parse");
        assert_eq!(
            parsed.steps[0].states[0].to_bits(),
            (1.0f64 / 3.0).to_bits()
        );
        assert_eq!(parsed.steps[0].states[2].to_bits(), (-0.0f64).to_bits());
        assert_eq!(parsed, trace);
    }

    #[test]
    fn a_missing_provenance_header_is_rejected() {
        let text = sample().to_text().replace("# gfortran: 15.2.0\n", "");
        assert_eq!(
            Trace::parse(&text),
            Err(TraceError::MissingHeader("gfortran"))
        );
    }

    #[test]
    fn a_truncated_data_line_is_rejected() {
        let mut text = sample().to_text();
        // Drop the last value from the final line.
        let cut = text.trim_end().rfind(' ').expect("a space");
        text.truncate(cut);
        text.push('\n');
        assert!(matches!(
            Trace::parse(&text),
            Err(TraceError::WrongFieldCount { .. })
        ));
    }

    #[test]
    fn a_corrupted_hex_word_is_rejected() {
        let text = sample()
            .to_text()
            .replacen("3ff8000000000000", "not_a_hex_word__", 1);
        assert!(matches!(
            Trace::parse(&text),
            Err(TraceError::BadHex { .. })
        ));
    }

    #[test]
    fn out_of_order_steps_are_rejected() {
        let text = sample().to_text().replacen("\n1 ", "\n7 ", 1);
        assert_eq!(
            Trace::parse(&text),
            Err(TraceError::StepOutOfOrder {
                expected: 1,
                got: 7
            })
        );
    }

    #[test]
    fn a_non_finite_value_is_rejected() {
        let mut trace = sample();
        trace.steps[0].derivatives[7] = f64::NAN;
        assert_eq!(
            Trace::parse(&trace.to_text()),
            Err(TraceError::NotFinite { step: 0, index: 57 })
        );
    }

    #[test]
    fn a_short_trace_fails_the_length_requirement() {
        assert_eq!(
            sample().require_full_length(),
            Err(TraceError::WrongStepCount {
                expected: STEPS,
                got: 2
            })
        );
    }

    #[test]
    fn the_layout_arithmetic_is_what_the_docs_claim() {
        assert_eq!(VALUES_PER_STEP, 142);
        assert_eq!(50 + 50 + 41 + 1, VALUES_PER_STEP);
    }
}
