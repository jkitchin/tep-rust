//! `tep`: run the Tennessee Eastman Process from a terminal.
//!
//! Deliberately small. Argument parsing is hand-written rather than pulled from
//! a crate, because the whole surface is eight flags and a subcommand, and
//! `tepsim-cli` is the one place a heavyweight dependency would show up in a
//! user's build. When the scenario engine lands and this needs to parse
//! scenario files, that decision is worth revisiting.

#![forbid(unsafe_code)]

use std::io::{self, BufWriter, Write};
use std::process::ExitCode;

use tepsim::{Csv, Decimating, Integrator, Outcome, Scenario, Simulation};

const USAGE: &str = "\
tep - Tennessee Eastman Process simulator

USAGE:
    tep run [OPTIONS]        Simulate and write CSV to stdout
    tep faults               List the twenty disturbances
    tep help                 This message

RUN OPTIONS:
    --fault <1-20>           Disturbance to inject (default: none)
    --hours <h>              Simulated duration (default: 48)
    --seed <n>               Generator word (default: 4651207995)
    --every <steps>          Sample every N steps (default: 180, i.e. 3 min)
    --decimate <n>           Write only every nth recorded sample
    --integrator <name>      euler (default, matches the original), rk4, dopri5
    --open-loop              Hold the valves instead of controlling
    --no-forced-idv12        Do not switch IDV(12) on at hour eight
    --labels                 Include ground-truth columns

EXAMPLES:
    tep run --hours 8 > normal.csv
    tep run --fault 4 --hours 24 --labels > idv4.csv
    tep run --fault 1 --seed 12345 --every 60 | head
    tep run --hours 4 --integrator rk4 > accurate.csv

NOTE ON --integrator:
    Only `euler` reproduces the original. The published data uses fixed-step
    explicit Euler at one second and carries about 1% of integration error
    against an accurate solution; rk4 and dopri5 remove that error and so give
    different numbers.
";

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    match args.first().map(String::as_str) {
        None | Some("help" | "--help" | "-h") => {
            print!("{USAGE}");
            ExitCode::SUCCESS
        }
        Some("faults") => {
            list_faults();
            ExitCode::SUCCESS
        }
        Some("run") => match parse_run(&args[1..]) {
            Ok(options) => run(options),
            Err(message) => fail(&message),
        },
        Some(other) => fail(&format!("unknown command `{other}`")),
    }
}

fn fail(message: &str) -> ExitCode {
    eprintln!("tep: {message}");
    eprintln!("try `tep help`");
    ExitCode::FAILURE
}

fn list_faults() {
    println!("The twenty disturbances of teprob.f. IDV(21) exists in later");
    println!("versions of the model and not in this one.\n");
    for fault in &tepsim::tepsim_core::FAULTS {
        println!("  IDV({:<2}) {}", fault.index, fault.published);
        println!("          {:?}, {}", fault.shape, fault.line);
        if !fault.effect.is_empty() && fault.effect != fault.published {
            println!("          {}", fault.effect);
        }
    }
}

/// What `tep run` was asked to do.
#[derive(Debug)]
struct Options {
    scenario: Scenario,
    labels: bool,
    decimate: usize,
}

fn parse_run(args: &[String]) -> Result<Options, String> {
    let mut scenario = Scenario::baseline();
    let mut labels = false;
    let mut decimate = 1_usize;
    let mut rest = args.iter();

    while let Some(flag) = rest.next() {
        let mut value = || {
            rest.next()
                .ok_or_else(|| format!("`{flag}` needs a value"))
                .cloned()
        };
        match flag.as_str() {
            "--fault" => {
                let raw = value()?;
                let n: usize = raw
                    .parse()
                    .map_err(|_| format!("`--fault {raw}` is not a number"))?;
                if !(1..=tepsim::DISTURBANCES).contains(&n) {
                    return Err(format!(
                        "`--fault {n}` is out of range: this model has {} disturbances",
                        tepsim::DISTURBANCES
                    ));
                }
                scenario = scenario.with_fault(n);
            }
            "--hours" => {
                let raw = value()?;
                let hours: f64 = raw
                    .parse()
                    .map_err(|_| format!("`--hours {raw}` is not a number"))?;
                // Written so a NaN falls into the guard rather than past it.
                if hours.is_nan() || hours <= 0.0 {
                    return Err(format!("`--hours {hours}` must be positive"));
                }
                scenario = scenario.with_hours(hours);
            }
            "--seed" => {
                let raw = value()?;
                let seed: f64 = raw
                    .parse()
                    .map_err(|_| format!("`--seed {raw}` is not a number"))?;
                if seed.is_nan() || seed <= 0.0 {
                    return Err("`--seed` must be positive".to_string());
                }
                scenario = scenario.with_seed(seed);
            }
            "--every" => {
                let raw = value()?;
                let steps: usize = raw
                    .parse()
                    .map_err(|_| format!("`--every {raw}` is not a number"))?;
                if steps == 0 {
                    return Err("`--every 0` would sample nothing".to_string());
                }
                scenario = scenario.sampling_every(steps);
            }
            "--decimate" => {
                let raw = value()?;
                let factor: usize = raw
                    .parse()
                    .map_err(|_| format!("`--decimate {raw}` is not a number"))?;
                if factor == 0 {
                    return Err("`--decimate 0` would keep nothing".to_string());
                }
                decimate = factor;
            }
            "--integrator" => {
                let raw = value()?;
                let method = Integrator::parse(&raw).ok_or_else(|| {
                    format!("`--integrator {raw}` is not one of euler, rk4, dopri5")
                })?;
                scenario = scenario.with_integrator(method);
            }
            "--open-loop" => scenario = scenario.open_loop(),
            "--no-forced-idv12" => scenario.driver_forces_idv12 = false,
            "--labels" => labels = true,
            other => return Err(format!("unknown option `{other}`")),
        }
    }
    Ok(Options {
        scenario,
        labels,
        decimate,
    })
}

fn run(options: Options) -> ExitCode {
    let Options {
        scenario,
        labels,
        decimate,
    } = options;

    eprintln!(
        "tep: {:.1} h, {} steps, sampling every {} ({} rows), seed {}, {}",
        scenario.hours,
        scenario.steps(),
        scenario.sample_every,
        scenario.samples().div_ceil(decimate),
        scenario.seed,
        if scenario.controlled {
            "closed loop"
        } else {
            "open loop"
        }
    );
    if !scenario.integrator.is_faithful() {
        eprintln!(
            "tep: integrator {}, which does NOT reproduce the original; only \
             euler does",
            scenario.integrator.name()
        );
    }

    // Streamed through the library's own CSV sink rather than collected first.
    // A 48-hour run sampled at every step is 172,800 rows of 53 channels, and
    // there is no reason for the process ever to hold them.
    //
    // `Decimating` wraps the CSV sink and not the writer: it has to see whole
    // samples, and decimating formatted text would drop parts of rows.
    let stdout = io::stdout();
    let mut writer = Adapter::new(BufWriter::new(stdout.lock()));
    let outcome = {
        let mut csv = Csv::new(&mut writer);
        if labels {
            csv = csv.with_labels();
        }
        let mut sink = Decimating::new(csv, decimate);
        Simulation::new(scenario).run_into(&mut sink)
    };

    // A closed pipe is how `| head` ends, and it is not a failure.
    if writer.broken_pipe() {
        return ExitCode::SUCCESS;
    }
    if let Err(error) = writer.finish() {
        if error.kind() == io::ErrorKind::BrokenPipe {
            return ExitCode::SUCCESS;
        }
        return fail(&format!("writing output: {error}"));
    }
    if let Some(error) = writer.error() {
        return fail(&format!("writing output: {error}"));
    }

    match outcome {
        Outcome::Completed => {
            eprintln!("tep: completed");
            ExitCode::SUCCESS
        }
        Outcome::Tripped { step, hours, cause } => {
            // Not a failure exit: a trip is a result, and the frozen samples
            // after it are part of what the original produces.
            eprintln!(
                "tep: the plant tripped at step {step} ({hours:.3} h) on {cause:?}; \
                 the plant is frozen after the trip"
            );
            ExitCode::SUCCESS
        }
        Outcome::SolveFailed { step } => {
            eprintln!("tep: a temperature solve failed to converge at step {step}");
            ExitCode::FAILURE
        }
    }
}

/// Bridges `core::fmt::Write`, which the library's sinks use, to
/// `std::io::Write`, which stdout is.
///
/// The library cannot use `io::Write`: the same code compiles to wasm32 for the
/// browser, where there is no `std`. Holding the first error rather than
/// returning it is forced by `fmt::Write`'s signature, which is why
/// [`Adapter::error`] exists.
struct Adapter<W> {
    out: W,
    error: Option<io::Error>,
}

impl<W: Write> Adapter<W> {
    const fn new(out: W) -> Self {
        Self { out, error: None }
    }

    fn broken_pipe(&self) -> bool {
        self.error
            .as_ref()
            .is_some_and(|e| e.kind() == io::ErrorKind::BrokenPipe)
    }

    const fn error(&self) -> Option<&io::Error> {
        self.error.as_ref()
    }

    fn finish(&mut self) -> io::Result<()> {
        self.out.flush()
    }
}

impl<W: Write> core::fmt::Write for Adapter<W> {
    fn write_str(&mut self, s: &str) -> core::fmt::Result {
        if self.error.is_some() {
            return Err(core::fmt::Error);
        }
        match self.out.write_all(s.as_bytes()) {
            Ok(()) => Ok(()),
            Err(error) => {
                self.error = Some(error);
                Err(core::fmt::Error)
            }
        }
    }
}

#[cfg(test)]
#[allow(
    clippy::float_cmp,
    reason = "parsed flags are compared to the exact values they were given"
)]
mod tests {
    use super::*;

    fn parse(args: &[&str]) -> Result<Options, String> {
        let owned: Vec<String> = args.iter().map(|s| (*s).to_string()).collect();
        parse_run(&owned)
    }

    #[test]
    fn no_arguments_is_the_baseline() {
        let options = parse(&[]).expect("parses");
        assert_eq!(options.scenario, Scenario::baseline());
        assert!(!options.labels);
        assert_eq!(options.decimate, 1);
    }

    #[test]
    fn every_flag_reaches_the_scenario() {
        let options = parse(&[
            "--fault",
            "4",
            "--hours",
            "12",
            "--seed",
            "99",
            "--every",
            "60",
            "--decimate",
            "5",
            "--open-loop",
            "--no-forced-idv12",
            "--labels",
        ])
        .expect("parses");

        assert_eq!(
            options.scenario.active_faults().collect::<Vec<_>>(),
            vec![4]
        );
        assert_eq!(options.scenario.hours, 12.0);
        assert_eq!(options.scenario.seed, 99.0);
        assert_eq!(options.scenario.sample_every, 60);
        assert_eq!(options.decimate, 5);
        assert!(!options.scenario.controlled);
        assert!(!options.scenario.driver_forces_idv12);
        assert!(options.labels);
        // Unspecified, so the faithful default.
        assert_eq!(options.scenario.integrator, Integrator::Euler);

        let rk4 = parse(&["--integrator", "rk4"]).expect("parses");
        assert_eq!(rk4.scenario.integrator, Integrator::Rk4);
        assert!(!rk4.scenario.integrator.is_faithful());
    }

    /// Bad input is refused with a message that says what to do, rather than
    /// panicking or silently doing something else.
    #[test]
    fn bad_input_is_refused_rather_than_guessed_at() {
        for (args, expected) in [
            (vec!["--fault", "0"], "out of range"),
            (vec!["--fault", "21"], "out of range"),
            (vec!["--fault", "x"], "not a number"),
            (vec!["--fault"], "needs a value"),
            (vec!["--hours", "0"], "must be positive"),
            (vec!["--hours", "-3"], "must be positive"),
            (vec!["--seed", "0"], "must be positive"),
            (vec!["--every", "0"], "would sample nothing"),
            (vec!["--decimate", "0"], "would keep nothing"),
            (vec!["--integrator", "heun"], "not one of"),
            (vec!["--nonsense"], "unknown option"),
        ] {
            let error = parse(&args).expect_err(&format!("{args:?} should be refused"));
            assert!(
                error.contains(expected),
                "{args:?} gave {error:?}, which does not mention {expected:?}"
            );
        }
    }

    /// Fault 20 is the last one this model has, and 21 is the first it does
    /// not. Pinned here because the CLI is where a user meets that boundary.
    #[test]
    fn the_fault_range_matches_the_model() {
        assert!(parse(&["--fault", "20"]).is_ok());
        assert!(parse(&["--fault", "21"]).is_err());
        assert_eq!(tepsim::DISTURBANCES, 20);
    }

    #[test]
    fn the_usage_text_mentions_every_flag_that_exists() {
        for flag in [
            "--fault",
            "--hours",
            "--seed",
            "--every",
            "--decimate",
            "--integrator",
            "--open-loop",
            "--no-forced-idv12",
            "--labels",
        ] {
            assert!(USAGE.contains(flag), "`{flag}` is undocumented");
        }
    }
}
