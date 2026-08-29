# Injecting a fault, and finding it in the data

A disturbance is one line: `Scenario::fault(n)` is the baseline with `IDV(n)`
switched on, one-based exactly as the Fortran's `IDV(n)` is. This tutorial
switches on `IDV(4)`, the step change in the reactor cooling water inlet
temperature, and then goes looking for it.

The interesting part is that it is nearly invisible where you would first look.

```rust,ignore
use tepsim::{Scenario, Simulation};

/// Mean and the largest absolute departure from a reference series.
fn compare(name: &str, base: &[f64], faulted: &[f64]) {
    let mean = |x: &[f64]| x.iter().sum::<f64>() / x.len() as f64;
    let worst = base
        .iter()
        .zip(faulted)
        .map(|(b, f)| (f - b).abs())
        .fold(0.0_f64, f64::max);
    println!(
        "  {name:<28} base {:>9.3}   fault {:>9.3}   worst gap {:>8.3}",
        mean(base),
        mean(faulted),
        worst
    );
}

fn main() {
    let hours = 8.0;
    let base = Simulation::new(Scenario::baseline().with_hours(hours)).run();
    let faulted = Simulation::new(Scenario::fault(4).with_hours(hours)).run();

    println!("--- IDV(4): reactor cooling water inlet temperature, step ---");
    println!("samples: {} each", base.samples.len());
    compare(
        "XMEAS(9) reactor temp",
        &base.measurement(9),
        &faulted.measurement(9),
    );
    compare(
        "XMEAS(21) reactor cw out",
        &base.measurement(21),
        &faulted.measurement(21),
    );
    compare(
        "XMV(10) reactor cw valve",
        &base.manipulated(10),
        &faulted.manipulated(10),
    );

    println!();
    println!("--- the first hour, sample by sample ---");
    let temp = faulted.measurement(9);
    let valve = faulted.manipulated(10);
    for i in [0, 1, 2, 3, 4, 6, 9, 14, 19] {
        println!(
            "  {:>6.3} h   XMEAS(9) {:>8.4}   XMV(10) {:>7.3}",
            faulted.samples[i].hours, temp[i], valve[i]
        );
    }

    println!();
    println!("--- ground truth ---");
    let first = &faulted.samples[0];
    println!("  faulted:      {}", first.labels.faulted());
    println!(
        "  faults:       {:?}",
        first.labels.faults().collect::<Vec<_>>()
    );
    println!("  since onset:  {:?} h", first.labels.since_onset[3]);

    println!();
    println!("--- the driver's IDV(12), which is no longer the default ---");
    let forced = Simulation::new(Scenario::faithful().with_hours(9.0).with_fault(4)).run();
    for sample in &forced.samples {
        let faults: Vec<usize> = sample.labels.faults().collect();
        if faults.len() > 1 {
            println!(
                "  faithful: first sample with two faults at {:.3} h {faults:?}",
                sample.hours
            );
            break;
        }
    }
    let plain = Simulation::new(Scenario::fault(4).with_hours(9.0)).run();
    let ever: Vec<usize> = plain
        .samples
        .last()
        .expect("a run has samples")
        .labels
        .faults()
        .collect();
    println!("  default:  at 9 h, still only {ever:?}");
}
```

```text
--- IDV(4): reactor cooling water inlet temperature, step ---
samples: 160 each
  XMEAS(9) reactor temp        base   120.401   fault   120.401   worst gap    0.161
  XMEAS(21) reactor cw out     base    94.597   fault    94.597   worst gap    0.214
  XMV(10) reactor cw valve     base    41.126   fault    44.901   worst gap    5.758

--- the first hour, sample by sample ---
   0.050 h   XMEAS(9) 120.5480   XMV(10)  46.995
   0.100 h   XMEAS(9) 120.4603   XMV(10)  44.373
   0.150 h   XMEAS(9) 120.4306   XMV(10)  44.572
   0.200 h   XMEAS(9) 120.4446   XMV(10)  45.563
   0.250 h   XMEAS(9) 120.4261   XMV(10)  45.382
   0.350 h   XMEAS(9) 120.3864   XMV(10)  44.318
   0.500 h   XMEAS(9) 120.3842   XMV(10)  45.026
   0.750 h   XMEAS(9) 120.3771   XMV(10)  44.616
   1.000 h   XMEAS(9) 120.3385   XMV(10)  43.543

--- ground truth ---
  faulted:      true
  faults:       [4]
  since onset:  Some(0.04972222222222225) h

--- the driver's IDV(12), which is no longer the default ---
  faithful: first sample with two faults at 8.000 h [4, 12]
  default:  at 9 h, still only [4]
```

## The fault is in the valve, not in the temperature

Over eight hours the mean reactor temperature is 120.401 degrees with the fault
and 120.401 degrees without it. To three decimal places the disturbance has left
no trace in the controlled variable at all, and the worst single departure over
160 samples is 0.16 degrees, which is inside the ordinary variation of the
fault-free plant.

Meanwhile the reactor cooling water valve has moved from a mean of 41.1% open to
44.9% open, and its worst departure is 5.8 percentage points. That is the whole
story of `IDV(4)`: the cooling water arrives hotter, the temperature controller
notices immediately, and it opens the valve until the temperature comes back.
The controller has converted a disturbance in a measured variable into a
disturbance in a manipulated one, which is what a controller is for.

The lesson generalises past this fault. On a plant under closed-loop control the
evidence of a disturbance often sits in the manipulated variables rather than in
the measurements, and a monitoring scheme that watches only `XMEAS` is looking
in the place the controller has been busy cleaning up. Both are in the row that
`Sample::row` returns, and the detector in the next tutorial uses all 53.

## Ground truth

`Sample::labels` is a record of what was actually true at that instant:
`active` is the twenty flags, `faults()` iterates the one-based indices of the
live ones, and `since_onset[i]` is how long `IDV(i + 1)` has been running.

The original records nothing of the sort. Every detection delay in the
literature is measured against an onset the author knew from the experimental
protocol rather than from the data, which is fine until two papers disagree
about where sample 160 falls. Here the onset is in the file.

Note the value: 0.04972222222222225 hours at the first sample. That is not
"about 0.05". The fault was live from the first step, and the label carries the
same simulated clock the row does, so it is 179 seconds. Precision here is what
lets a delay of one sample be distinguished from a delay of zero.

## The disturbance you did not ask for, and no longer get

Run `Scenario::faithful()` for nine hours instead of eight and a second fault
appears at hour eight without being asked for. That is not a bug in the
scheduler; it is `temain_mod.f:366-368`, which switches `IDV(12)` on at eight
hours whatever the scenario said.

It was the default here until 2026-08-28, on the belief that reproducing `d01`
through `d21` required it. Tier 7 showed the opposite: every `dNN_te` file
except `d12_te` sits at the nominal operating point straight across row 160,
which is hour eight, so the published files were made with that line replaced.
The default now follows the files, and the driver's version is a named
constructor away.

It is delta D-011 in the register, and it used to be the single most common way
a comparison against the published files went quietly wrong, because it was the
default. It is not any more: since 2026-08-28 a request for `IDV(4)` gets
`IDV(4)`, and `Scenario::faithful()` is how you ask for the driver's version, as
the last block above does. Say which one you used when you report numbers. The
labels make the difference visible either way, which is the point of recording
ground truth rather than assuming it.

## From the command line

```console
$ tep run --fault 4 --hours 8 --labels
```

`--labels` adds the `fault` and `hours_since_onset` columns to the CSV, so the
ground truth travels with the data.
