//! `tep`: run the Tennessee Eastman Process from a terminal.
//!
//! Deliberately small. Argument parsing is hand-written rather than pulled from
//! a crate, because the whole surface is four flags and a subcommand, and
//! `tepsim-cli` is the one place a heavyweight dependency would show up in a
//! user's build. When the scenario engine lands in B-0054 and this needs to
//! parse scenario files, that decision is worth revisiting.

#![forbid(unsafe_code)]

use std::io::{self, BufWriter, Write};
use std::process::ExitCode;

use tepsim::{Integrator, Outcome, Run, Scenario, Simulation, channel_names};

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
    --open-loop              Hold the valves instead of controlling
    --no-forced-idv12        Do not switch IDV(12) on at hour eight
    --labels                 Include ground-truth columns
    --integrator <name>      euler (default, matches the original), rk4, dopri5

EXAMPLES:
    tep run --hours 8 > normal.csv
    tep run --fault 4 --hours 24 --labels > idv4.csv
    tep run --fault 1 --seed 12345 --every 60 | head
    tep run --hours 4 --integrator rk4 > accurate.csv

NOTE ON --integrator:
    Only `euler` reproduces the original. The published data uses fixed-step
    explicit Euler at one second and carries about 1% of integration error
    against an accurate solution; rk4 and dopri5 remove that error and so give
    different numbers. See the book's validation chapter.
";

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let command = args.first().map(String::as_str);

    match command {
        None | Some("help" | "--help" | "-h") => {
            print!("{USAGE}");
            ExitCode::SUCCESS
        }
        Some("faults") => {
            list_faults();
            ExitCode::SUCCESS
        }
        Some("run") => match parse_run(&args[1..]) {
            Ok((scenario, labels)) => run(scenario, labels),
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

/// Parse `tep run`'s flags.
///
/// Returns the scenario and whether to emit label columns.
fn parse_run(args: &[String]) -> Result<(Scenario, bool), String> {
    let mut scenario = Scenario::baseline();
    let mut labels = false;
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
    Ok((scenario, labels))
}

fn run(scenario: Scenario, labels: bool) -> ExitCode {
    let steps = scenario.steps();
    eprintln!(
        "tep: {:.1} h, {steps} steps, sampling every {} ({} rows), seed {}, {}",
        scenario.hours,
        scenario.sample_every,
        scenario.samples(),
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

    let finished = Simulation::new(scenario).run();
    if let Err(error) = write_csv(&finished, labels) {
        // A closed pipe is how `| head` ends, and it is not a failure.
        if error.kind() == io::ErrorKind::BrokenPipe {
            return ExitCode::SUCCESS;
        }
        return fail(&format!("writing output: {error}"));
    }

    match finished.outcome {
        Outcome::Completed => {
            eprintln!("tep: completed, {} rows", finished.samples.len());
            ExitCode::SUCCESS
        }
        Outcome::Tripped { step, hours, cause } => {
            // Not a failure exit: a trip is a result, and the frozen samples
            // after it are part of what the original produces.
            eprintln!(
                "tep: the plant tripped at step {step} ({hours:.3} h) on {cause:?}; \
                 {} rows written, the plant frozen after the trip",
                finished.samples.len()
            );
            ExitCode::SUCCESS
        }
        Outcome::SolveFailed { step } => {
            eprintln!("tep: a temperature solve failed to converge at step {step}");
            ExitCode::FAILURE
        }
    }
}

fn write_csv(run: &Run, labels: bool) -> io::Result<()> {
    let stdout = io::stdout();
    let mut out = BufWriter::new(stdout.lock());

    write!(out, "step,hours")?;
    for name in channel_names() {
        write!(out, ",{name}")?;
    }
    if labels {
        write!(out, ",fault,hours_since_onset")?;
    }
    writeln!(out)?;

    for sample in &run.samples {
        write!(out, "{},{:.6}", sample.step, sample.hours)?;
        for value in sample.row() {
            // Seventeen significant digits round-trips an f64 exactly, which
            // is what makes a CSV written here reproducible rather than
            // approximately reproducible.
            write!(out, ",{value:.17e}")?;
        }
        if labels {
            let active: Vec<String> = sample.labels.faults().map(|n| n.to_string()).collect();
            let since = sample
                .labels
                .faults()
                .next()
                .and_then(|n| sample.labels.since_onset[n - 1]);
            write!(out, ",{}", active.join(" "))?;
            match since {
                Some(hours) => write!(out, ",{hours:.6}")?,
                None => write!(out, ",")?,
            }
        }
        writeln!(out)?;
    }
    out.flush()
}

#[cfg(test)]
#[allow(
    clippy::float_cmp,
    reason = "parsed flags are compared to the exact values they were given"
)]
mod tests {
    use super::*;

    fn parse(args: &[&str]) -> Result<(Scenario, bool), String> {
        let owned: Vec<String> = args.iter().map(|s| (*s).to_string()).collect();
        parse_run(&owned)
    }

    #[test]
    fn no_arguments_is_the_baseline() {
        let (scenario, labels) = parse(&[]).expect("parses");
        assert_eq!(scenario, Scenario::baseline());
        assert!(!labels);
    }

    #[test]
    fn every_flag_reaches_the_scenario() {
        let (scenario, labels) = parse(&[
            "--fault",
            "4",
            "--hours",
            "12",
            "--seed",
            "99",
            "--every",
            "60",
            "--open-loop",
            "--no-forced-idv12",
            "--labels",
        ])
        .expect("parses");

        assert_eq!(scenario.active_faults().collect::<Vec<_>>(), vec![4]);
        assert_eq!(scenario.hours, 12.0);
        assert_eq!(scenario.seed, 99.0);
        assert_eq!(scenario.sample_every, 60);
        assert!(!scenario.controlled);
        assert!(!scenario.driver_forces_idv12);
        assert!(labels);
        assert_eq!(scenario.integrator, Integrator::Euler);

        let (rk4, _) = parse(&["--integrator", "rk4"]).expect("parses");
        assert_eq!(rk4.integrator, Integrator::Rk4);
        assert!(!rk4.integrator.is_faithful());
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
            "--integrator",
            "--fault",
            "--hours",
            "--seed",
            "--every",
            "--open-loop",
            "--no-forced-idv12",
            "--labels",
        ] {
            assert!(USAGE.contains(flag), "`{flag}` is undocumented");
        }
    }
}
