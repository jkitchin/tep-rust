# tepsim

The Tennessee Eastman Process, as a Python package with no C dependencies.

A pure-Rust port of Downs and Vogel's 1993 challenge problem, taken from the
original Fortran rather than from any later reimplementation, and validated
against that Fortran through a ten-tier ladder. Given the same `exp` and `pow`,
a complete 48-hour closed-loop run of 172,800 integrator steps is bit-identical
in all 41 measurements and all 12 manipulated variables.

```python
import tepsim as tep

run = tep.Simulation(tep.Scenario.baseline(seed=42, hours=48)).run()

run.to_numpy()        # (960, 53) float64, XMEAS(1..41) then XMV(1..12)
run.measurement(7)    # XMEAS(7), reactor pressure
run.columns()         # {'XMEAS_1_A_feed': array([...]), ...}
run.outcome           # 'completed'
```

Faults are the twenty `IDV` disturbances, and the table says what each one
actually does rather than repeating the five the original header calls
"Unknown":

```python
for f in tep.faults():
    print(f.index, f.shape, f.published, "|", f.effect)

run = tep.Simulation(tep.Scenario.fault(1, hours=8)).run()
run.labels()["since_onset"]   # hours since each disturbance came on
```

## Reproducibility

A run is a pure function of its scenario: no clock, no global state, no
randomness outside the seeded generator. The same scenario gives bit-identical
output on x86-64, aarch64 and wasm32, so a dataset is reproducible from its
description rather than from a file. `repr(scenario)` round-trips.

## Threads

`Simulation.run()` releases the GIL for the whole integration, so an ensemble
is a thread pool and nothing has to be pickled.

```python
from concurrent.futures import ThreadPoolExecutor

sims = [tep.Simulation(tep.Scenario.fault(n, hours=8)) for n in range(1, 21)]
with ThreadPoolExecutor() as pool:
    runs = list(pool.map(tep.Simulation.run, sims))
```

## Arrays

Every array a `Run` hands out is a read-only view over one buffer, which is
filled once when the run finishes and moved into NumPy rather than copied.
`to_numpy()` returns the same object every call; `measurement`, `manipulated`,
`column` and `columns` are strided views into it. Call `.copy()` for something
writable, or `numpy.ascontiguousarray` if contiguity matters.

## Licence

New work is MIT. The portions derived from the original Tennessee Eastman
Fortran are under the University of Illinois/NCSA Open Source License, whose
attribution conditions survive into binary distributions. Both texts and
`NOTICE.md` ship inside the wheel.
