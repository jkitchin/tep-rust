"""Type stubs for the compiled half of `tepsim`.

Kept beside `__init__.py` rather than as `__init__.pyi`, so the stub describes
the extension module and the package's own re-exports stay checkable Python.
"""

from typing import Any, Dict, Optional, Sequence, Tuple

import numpy as np
from numpy.typing import NDArray

__version__: str

MEASUREMENTS: int
MANIPULATED: int
CHANNELS: int
DISTURBANCES: int

DEFAULT_SEED: float
DEFAULT_STEP_HOURS: float
DEFAULT_SAMPLE_EVERY: int
FORCED_DISTURBANCE_STEP: int

def channel_names() -> Tuple[str, ...]:
    """Short, stable names for the 53 recorded channels, in row order."""
    ...

def faults() -> Tuple[Fault, ...]:
    """The twenty disturbances, in `IDV` order, with what each one does."""
    ...

class Fault:
    """One of the twenty disturbances."""

    @property
    def index(self) -> int: ...
    @property
    def published(self) -> str: ...
    @property
    def effect(self) -> str: ...
    @property
    def shape(self) -> str: ...
    @property
    def channels(self) -> Tuple[int, ...]: ...
    @property
    def spiking(self) -> bool: ...
    @property
    def valves(self) -> Tuple[int, ...]: ...
    @property
    def line(self) -> str: ...
    @property
    def affects_the_plant(self) -> bool: ...
    def __repr__(self) -> str: ...

class Scenario:
    """A complete description of a run."""

    def __init__(
        self,
        *,
        seed: float = ...,
        hours: float = ...,
        step_hours: float = ...,
        sample_every: int = ...,
        faults: Sequence[int] = ...,
        controlled: bool = ...,
        driver_forces_idv12: bool = ...,
        trip_ends_the_run: bool = ...,
    ) -> None: ...
    @staticmethod
    def baseline(
        *,
        seed: float = ...,
        hours: float = ...,
        step_hours: float = ...,
        sample_every: int = ...,
        controlled: bool = ...,
        driver_forces_idv12: bool = ...,
        trip_ends_the_run: bool = ...,
    ) -> Scenario: ...
    @staticmethod
    def fault(
        n: int,
        *,
        seed: float = ...,
        hours: float = ...,
        step_hours: float = ...,
        sample_every: int = ...,
        controlled: bool = ...,
        driver_forces_idv12: bool = ...,
        trip_ends_the_run: bool = ...,
    ) -> Scenario: ...
    @property
    def seed(self) -> float: ...
    @property
    def hours(self) -> float: ...
    @property
    def step_hours(self) -> float: ...
    @property
    def sample_every(self) -> int: ...
    @property
    def controlled(self) -> bool: ...
    @property
    def driver_forces_idv12(self) -> bool: ...
    @property
    def trip_ends_the_run(self) -> bool: ...
    @property
    def faults(self) -> Tuple[int, ...]: ...
    @property
    def steps(self) -> int: ...
    @property
    def samples(self) -> int: ...
    @property
    def digest(self) -> str: ...
    def to_text(self) -> str:
        """This scenario as one line of canonical, versioned text."""
        ...

    @staticmethod
    def from_text(text: str) -> Scenario:
        """Parse `to_text`. Strict: anything wrong is a ValueError that says what."""
        ...

    def with_seed(self, seed: float) -> Scenario: ...
    def with_hours(self, hours: float) -> Scenario: ...
    def with_fault(self, n: int) -> Scenario: ...
    def sampling_every(self, steps: int) -> Scenario: ...
    def open_loop(self) -> Scenario: ...
    def __eq__(self, other: Any) -> bool: ...
    def __repr__(self) -> str: ...

class Simulation:
    """A simulation ready to run."""

    def __init__(self, scenario: Optional[Scenario] = ...) -> None: ...
    @property
    def scenario(self) -> Scenario: ...
    def run(self) -> Run:
        """Run the whole scenario, releasing the GIL for the integration."""
        ...

    def __repr__(self) -> str: ...

class Run:
    """A finished run: its scenario, its samples as arrays, and how it ended."""

    @property
    def scenario(self) -> Scenario: ...
    @property
    def outcome(self) -> str:
        """`'completed'`, `'tripped'` or `'solve_failed'`."""
        ...

    @property
    def tripped_at(self) -> Optional[int]: ...
    @property
    def tripped_hours(self) -> Optional[float]: ...
    @property
    def trip_cause(self) -> Optional[str]: ...
    @property
    def solve_failed_at(self) -> Optional[int]: ...
    @property
    def hours(self) -> NDArray[np.float64]: ...
    @property
    def steps(self) -> NDArray[np.int64]: ...
    def to_numpy(self) -> NDArray[np.float64]:
        """The whole run as one read-only `(n_samples, 53)` float64 array."""
        ...

    def column(self, channel: int) -> NDArray[np.float64]:
        """One channel, zero-based over the 53. A strided view, not a copy."""
        ...

    def measurement(self, n: int) -> NDArray[np.float64]:
        """One measurement, one-based as `XMEAS(n)` is."""
        ...

    def manipulated(self, n: int) -> NDArray[np.float64]:
        """One manipulated variable, one-based as `XMV(n)` is."""
        ...

    def columns(self) -> Dict[str, NDArray[np.float64]]:
        """Every channel by name, in row order."""
        ...

    def labels(self) -> Dict[str, NDArray[Any]]:
        """Ground truth: `'active'` (bool) and `'since_onset'` (float64)."""
        ...

    def __len__(self) -> int: ...
    def __repr__(self) -> str: ...
