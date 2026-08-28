# A first run, and the 53 channels

The whole library is three types. A `Scenario` says what to simulate, a
`Simulation` does it, and a `Run` holds what came out. This tutorial runs two
hours of the fault-free plant and takes the output apart.

Every program on these four pages was compiled and run against this commit, and
the output shown under each one is the output it actually produced. Numbers are
quoted rather than described, because a number can be compared against the one
your own copy produces and a description cannot.

Put this in a crate that depends on `tepsim`:

```rust,ignore
use tepsim::{Scenario, Simulation, channel_names};

fn main() {
    println!("--- first run ---");
    let run = Simulation::new(Scenario::baseline().with_hours(2.0)).run();
    println!("samples:  {}", run.samples.len());
    println!("outcome:  {:?}", run.outcome);
    println!("first at: {:.4} h", run.samples[0].hours);
    println!("last at:  {:.4} h", run.samples[run.samples.len() - 1].hours);

    println!();
    println!("--- one row ---");
    let names = channel_names();
    let row = run.samples[0].row();
    println!("row length: {}", row.len());
    for channel in [0, 6, 8, 11, 14, 16, 39, 41, 52] {
        println!("  [{channel:>2}] {:<32} {:>10.4}", names[channel], row[channel]);
    }

    println!();
    println!("--- the analysers hold ---");
    let purge_b = run.measurement(30);
    let pressure = run.measurement(7);
    for i in 0..8 {
        println!(
            "  {:>6.3} h   XMEAS(7) {:>8.2}   XMEAS(30) {:>7.4}",
            run.samples[i].hours, pressure[i], purge_b[i]
        );
    }

    println!();
    println!("--- columns ---");
    let columns = run.columns();
    println!("columns: {}, each {} long", columns.len(), columns[0].len());
    let scenario = run.scenario;
    println!(
        "steps {} at {} h each, one sample every {}",
        scenario.steps(),
        scenario.step_hours,
        scenario.sample_every
    );
    println!(
        "digest {}",
        core::str::from_utf8(&scenario.digest_hex()).unwrap()
    );
}
```

```text
--- first run ---
samples:  40
outcome:  Completed
first at: 0.0497 h
last at:  1.9997 h

--- one row ---
row length: 53
  [ 0] XMEAS_1_A_feed                       0.2514
  [ 6] XMEAS_7_reactor_pressure          2706.1771
  [ 8] XMEAS_9_reactor_temperature        120.3866
  [11] XMEAS_12_separator_level            49.6902
  [14] XMEAS_15_stripper_level             50.9667
  [16] XMEAS_17_stripper_underflow         22.3050
  [39] XMEAS_40_product_G                  53.7240
  [41] XMV_1_D_feed_flow                   62.9625
  [52] XMV_12_agitator_speed               50.0000

--- the analysers hold ---
   0.050 h   XMEAS(7)  2706.18   XMEAS(30) 13.8229
   0.100 h   XMEAS(7)  2704.64   XMEAS(30) 13.8229
   0.150 h   XMEAS(7)  2704.19   XMEAS(30) 13.6879
   0.200 h   XMEAS(7)  2704.26   XMEAS(30) 13.6879
   0.250 h   XMEAS(7)  2706.26   XMEAS(30) 13.8244
   0.300 h   XMEAS(7)  2707.33   XMEAS(30) 13.8244
   0.350 h   XMEAS(7)  2705.33   XMEAS(30) 13.7758
   0.400 h   XMEAS(7)  2705.46   XMEAS(30) 13.7758

--- columns ---
columns: 53, each 40 long
steps 7200 at 0.0002777777777777778 h each, one sample every 180
digest b3415d9a395b8c70
```

## Why forty samples and not 7200

The integrator takes one step per simulated second, so two hours is 7200 steps.
A sample is written every 180 of them, which is the three-minute spacing
`temain_mod.f:401` writes at and the spacing of the published `d00` through
`d21` files. Two hours is therefore forty rows.

The first sample is at 0.0497 hours rather than at 0.05, and that is not an
off-by-one. `Sample::step` on that row is 180, and `Sample::hours` is the time
at which step 180 *began*, which is 179 seconds. The simulated clock is advanced
at the end of a step, after the row has been written, because that is the order
`temain_mod.f` writes in and a row that carried the post-step time would be
labelled with a clock the plant had not reached when it was measured. The last
row is at 1.9997 hours for the same reason. If you want the initial condition
itself it is `tepsim_core`'s nominal state, not a row of the run.

## What the 53 channels are

A row is the plant as an operator sees it: 41 measurements and then the 12
valve positions the controllers are holding. `Sample::row` returns exactly that,
measurements first, and `channel_names()` returns names in the same order, so a
CSV header and a matrix column can never disagree about which is which.

The 41 measurements split into two groups that behave quite differently.

`XMEAS(1)` through `XMEAS(22)` are continuous instruments: flows, pressures,
levels, temperatures, the compressor work. They are read every step and carry
Gaussian measurement noise whose standard deviation is the `XNS` table in
`teprob.f`.

`XMEAS(23)` through `XMEAS(41)` are the three gas chromatographs. They are not
instruments in the same sense at all: they sample on a schedule, take time to
run, and then hold the answer until the next result arrives. That is why
`XMEAS(30)`, the purge stream's B composition, repeats in pairs above while
`XMEAS(7)` moves at every row. The reactor feed and purge analysers run every
0.1 hours and the product analyser every 0.25, so at a 0.05-hour output cadence
each answer appears twice or five times. A detector that treats those repeats as
independent observations is counting the same measurement several times, and
that is worth knowing before it produces a p-value.

Both `0.1` literals in `teprob.f` are single precision, so the gas interval is
really 0.10000000149011612 and a step landing on exactly 0.1 does not sample.
That is faithfully reproduced here; see `tepsim_core::analysers` and the delta
register.

`XMV(1)` through `XMV(12)` are the manipulated variables. Eleven of them move.
`XMV(12)`, the agitator speed, sits at 50 forever, because the published control
scheme never touches it. A statistical model fitted to all 53 channels has to
cope with that constant column, and the next tutorial but one shows what
`tepsim-stats` does about it.

## Getting at the data

`Run::samples` is the record itself. Three views on top of it cover most uses.
`Run::measurement(n)` and `Run::manipulated(n)` take the one-based indices of
the original, so `run.measurement(7)` is `XMEAS(7)` and needs no mental
arithmetic. `Run::column(i)` takes the zero-based row position, and
`Run::columns()` returns all 53 at once, which is the shape a covariance matrix
or a chart wants.

## The digest

`Scenario::digest_hex` is a content hash over everything the run's output
depends on: the seed, the duration, the step, the cadence, the disturbances, the
control mode, the quirk flags, the schedule and the integrator. Two scenarios
that describe the same experiment produce the same sixteen characters, and two
that differ in any respect do not.

It is worth putting in a filename or a file header, because it turns "this is
the fault 4 data, I think" into something checkable. The next tutorials use it,
and `Scenario::to_text` writes out the whole scenario in one line that
`Scenario::from_text` reads back, so a dataset can carry its own description
rather than a memory of one.

## The same run from the command line

```console
$ tep run --hours 2
```

writes the same 40 rows as CSV on stdout, with seventeen significant digits, so
the file round-trips an `f64` exactly.
