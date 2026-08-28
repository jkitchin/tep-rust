"""Tennessee Eastman Process simulator.

A pure-Rust port of Downs and Vogel's 1993 challenge problem, taken from the
original Fortran rather than from any later reimplementation, and validated
against that Fortran through a ten-tier ladder. Given the same ``exp`` and
``pow``, a complete 48-hour closed-loop run of 172,800 integrator steps is
bit-identical in all 41 measurements and all 12 manipulated variables.

Build a :class:`Scenario`, run a :class:`Simulation`, get a :class:`Run`::

    import tepsim as tep

    run = tep.Simulation(tep.Scenario.baseline(seed=42, hours=48)).run()
    matrix = run.to_numpy()          # (960, 53) float64
    pressure = run.measurement(7)    # XMEAS(7), reactor pressure

A run is a pure function of its scenario. There is no clock, no global state,
and no randomness outside the seeded generator, so a dataset is reproducible
from its description rather than from a file.

``run()`` releases the GIL, so an ensemble parallelises with threads and needs
nothing pickled::

    from concurrent.futures import ThreadPoolExecutor

    sims = [tep.Simulation(tep.Scenario.fault(n, hours=8)) for n in range(1, 21)]
    with ThreadPoolExecutor() as pool:
        runs = list(pool.map(tep.Simulation.run, sims))

The arrays a :class:`Run` hands out are read-only views over one buffer, which
is filled once when the run finishes. Call ``.copy()`` for something writable.
"""

from ._tepsim import CHANNELS as CHANNELS
from ._tepsim import DEFAULT_SAMPLE_EVERY as DEFAULT_SAMPLE_EVERY
from ._tepsim import DEFAULT_SEED as DEFAULT_SEED
from ._tepsim import DEFAULT_STEP_HOURS as DEFAULT_STEP_HOURS
from ._tepsim import DISTURBANCES as DISTURBANCES
from ._tepsim import FORCED_DISTURBANCE_STEP as FORCED_DISTURBANCE_STEP
from ._tepsim import MANIPULATED as MANIPULATED
from ._tepsim import MEASUREMENTS as MEASUREMENTS
from ._tepsim import Fault as Fault
from ._tepsim import Run as Run
from ._tepsim import Scenario as Scenario
from ._tepsim import Simulation as Simulation
from ._tepsim import __version__ as __version__
from ._tepsim import channel_names as channel_names
from ._tepsim import faults as faults

__all__ = [
    "CHANNELS",
    "DEFAULT_SAMPLE_EVERY",
    "DEFAULT_SEED",
    "DEFAULT_STEP_HOURS",
    "DISTURBANCES",
    "FORCED_DISTURBANCE_STEP",
    "MANIPULATED",
    "MEASUREMENTS",
    "Fault",
    "Run",
    "Scenario",
    "Simulation",
    "__version__",
    "channel_names",
    "faults",
]
