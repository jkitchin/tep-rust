"""Tests for the Python bindings.

Run against an installed wheel, not the source tree:

    maturin build --release -m crates/tepsim-py/Cargo.toml -o dist
    pip install --force-reinstall dist/*.whl
    python -m pytest crates/tepsim-py/tests

Or in one step, if the current environment is disposable:

    maturin develop --release -m crates/tepsim-py/Cargo.toml
    python -m pytest crates/tepsim-py/tests

What is tested here is the *binding*, not the plant: that the arrays have the
shape, dtype, writeability and aliasing the documentation claims, that a
scenario determines its output exactly, and that `run()` really does let go of
the GIL. The plant itself is validated against the original Fortran by the Rust
test suite, which is where a numerical claim belongs.
"""

import ast
import os
import time
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pytest

import tepsim as tep


# A short run, for the tests that only care about shapes and identities. Two
# hours is 7200 integrator steps and 40 samples, and takes about 50 ms.
SHORT = tep.Scenario.baseline(hours=2)


def test_module_surface():
    assert isinstance(tep.__version__, str)
    assert (tep.MEASUREMENTS, tep.MANIPULATED, tep.CHANNELS) == (41, 12, 53)
    assert tep.DISTURBANCES == 20
    # `temain_mod.f:366-368` forces IDV(12) on at eight hours, one second per
    # step.
    assert tep.FORCED_DISTURBANCE_STEP == 8 * 3600


def test_the_package_ships_its_types():
    """`py.typed` and the stub have to be installed, not merely committed."""
    root = os.path.dirname(tep.__file__)
    assert os.path.exists(os.path.join(root, "py.typed"))
    assert os.path.exists(os.path.join(root, "_tepsim.pyi"))


def _stub():
    """The installed stub, parsed."""
    path = os.path.join(os.path.dirname(tep.__file__), "_tepsim.pyi")
    with open(path, encoding="utf-8") as handle:
        return ast.parse(handle.read())


def test_the_stub_lists_what_the_module_exports():
    """A stub nobody checks is a stub that quietly goes stale."""
    declared = set()
    for node in _stub().body:
        if isinstance(node, (ast.ClassDef, ast.FunctionDef)):
            declared.add(node.name)
        elif isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name):
            declared.add(node.target.id)

    assert declared == set(tep.__all__)


@pytest.mark.parametrize("name", ["Scenario", "Simulation", "Run", "Fault"])
def test_the_stub_describes_every_public_attribute(name):
    stubbed = {
        item.name
        for node in _stub().body
        if isinstance(node, ast.ClassDef) and node.name == name
        for item in node.body
        if isinstance(item, ast.FunctionDef) and not item.name.startswith("_")
    }
    actual = {
        attribute
        for attribute in dir(getattr(tep, name))
        if not attribute.startswith("_")
    }

    assert stubbed == actual


# ---------------------------------------------------------------------------
# Arrays: shape, dtype, and who owns the memory
# ---------------------------------------------------------------------------


def test_to_numpy_shape_and_dtype():
    run = tep.Simulation(SHORT).run()
    matrix = run.to_numpy()

    assert matrix.shape == (SHORT.samples, tep.CHANNELS) == (40, 53)
    assert matrix.dtype == np.float64
    assert matrix.flags.c_contiguous
    assert len(run) == 40


def test_arrays_are_read_only():
    run = tep.Simulation(SHORT).run()
    for array in (
        run.to_numpy(),
        run.hours,
        run.steps,
        run.measurement(7),
        run.manipulated(1),
        run.labels()["active"],
        run.labels()["since_onset"],
    ):
        assert not array.flags.writeable
    with pytest.raises(ValueError):
        run.to_numpy()[0, 0] = 0.0


def test_to_numpy_does_not_copy():
    """The matrix is built once; every call hands back that same object."""
    run = tep.Simulation(SHORT).run()
    assert run.to_numpy() is run.to_numpy()


def test_columns_are_views_of_the_matrix():
    """A channel is a strided view, not a fresh allocation."""
    run = tep.Simulation(SHORT).run()
    matrix = run.to_numpy()

    for column in (run.column(0), run.measurement(7), run.manipulated(12)):
        assert np.shares_memory(column, matrix)

    for column in run.columns().values():
        assert np.shares_memory(column, matrix)


def test_accessors_agree_with_the_matrix():
    run = tep.Simulation(SHORT).run()
    matrix = run.to_numpy()

    assert np.array_equal(run.measurement(1), matrix[:, 0])
    assert np.array_equal(run.measurement(41), matrix[:, 40])
    assert np.array_equal(run.manipulated(1), matrix[:, 41])
    assert np.array_equal(run.manipulated(12), matrix[:, 52])
    assert np.array_equal(run.column(7), matrix[:, 7])


def test_columns_are_keyed_by_channel_name_in_row_order():
    run = tep.Simulation(SHORT).run()
    columns = run.columns()
    names = tep.channel_names()

    assert len(names) == tep.CHANNELS
    assert tuple(columns) == names
    assert names[0].startswith("XMEAS_1_")
    assert names[41].startswith("XMV_1_")
    for channel, name in enumerate(names):
        assert np.array_equal(columns[name], run.column(channel))


def test_time_and_step_columns():
    run = tep.Simulation(SHORT).run()

    assert run.steps.dtype == np.int64
    # One-based, as `temain_mod.f`'s loop counter is, and one sample every
    # `sample_every` steps.
    assert np.array_equal(run.steps, np.arange(1, 41) * SHORT.sample_every)
    assert run.hours.shape == (40,)
    assert np.all(np.diff(run.hours) > 0)


# ---------------------------------------------------------------------------
# Determinism
# ---------------------------------------------------------------------------


def test_the_same_scenario_gives_identical_arrays():
    """Bit-identical, not close: a run is a pure function of its scenario."""
    scenario = tep.Scenario.baseline(seed=42, hours=4)
    first = tep.Simulation(scenario).run().to_numpy()
    second = tep.Simulation(scenario).run().to_numpy()

    assert np.array_equal(first, second)
    assert first.tobytes() == second.tobytes()


def test_running_twice_gives_the_same_run():
    """`run()` works on a copy, so the simulation stays usable."""
    simulation = tep.Simulation(tep.Scenario.baseline(seed=42, hours=2))
    first = simulation.run()
    second = simulation.run()

    assert np.array_equal(first.to_numpy(), second.to_numpy())
    assert first.outcome == second.outcome


def test_a_different_seed_gives_different_noise():
    """Determinism must not mean the seed is ignored."""
    base = tep.Scenario.baseline(seed=42, hours=4)
    other = base.with_seed(99)

    assert not np.array_equal(
        tep.Simulation(base).run().to_numpy(),
        tep.Simulation(other).run().to_numpy(),
    )


# ---------------------------------------------------------------------------
# The plant does what the scenario asked for
# ---------------------------------------------------------------------------


def test_a_fault_run_differs_from_the_baseline():
    hours = 4
    baseline = tep.Simulation(tep.Scenario.baseline(hours=hours)).run().to_numpy()
    faulted = tep.Simulation(tep.Scenario.fault(1, hours=hours)).run().to_numpy()

    assert baseline.shape == faulted.shape
    assert not np.array_equal(baseline, faulted)
    # IDV(1) steps the mixed feed's A fraction down, which the loops fight for
    # hours. The difference is large, not a rounding artefact.
    assert np.abs(faulted - baseline).max() > 1.0


def test_the_open_loop_plant_trips():
    """The clearest single statement of what the control layer does.

    With the valves held where they started, reactor pressure runs away and
    `teprob.f:703` fires after about three hours.
    """
    run = tep.Simulation(tep.Scenario.baseline(hours=6).open_loop()).run()

    assert run.outcome == "tripped"
    assert run.trip_cause == "reactor pressure high"
    assert 2.5 < run.tripped_hours < 3.5
    # The step is one-based and the recorded time is the time at the *start* of
    # that step, so the two differ by exactly one step.
    assert run.tripped_at == round(run.tripped_hours * 3600) + 1
    # The run ends at the trip, so it is shorter than its scenario planned.
    # Delta D-007, signed off 2026-08-28.
    assert len(run) < run.scenario.samples

    # `teprob.f:807-811` freezes the plant and keeps reporting instead, which
    # is what the published `d06` and `d18` files contain and what any
    # comparison against them needs.
    faithful = tep.Simulation(
        tep.Scenario.baseline(hours=6, trip_ends_the_run=False).open_loop()
    ).run()
    assert faithful.outcome == "tripped"
    assert len(faithful) == faithful.scenario.samples


def test_the_closed_loop_plant_does_not_trip():
    run = tep.Simulation(tep.Scenario.baseline(hours=6)).run()

    assert run.outcome == "completed"
    assert run.tripped_at is None
    assert run.trip_cause is None
    assert run.solve_failed_at is None


def test_a_trip_ends_the_run_by_default():
    """Delta D-007, signed off 2026-08-28."""
    scenario = tep.Scenario(hours=6, controlled=False)
    assert scenario.trip_ends_the_run
    run = tep.Simulation(scenario).run()

    assert run.outcome == "tripped"
    assert len(run) < scenario.samples


def test_the_faithful_configuration_freezes_instead():
    """`teprob.f:807-811`, which is what made the published frozen tails."""
    scenario = tep.Scenario(hours=6, controlled=False, trip_ends_the_run=False)
    run = tep.Simulation(scenario).run()

    assert run.outcome == "tripped"
    assert len(run) == scenario.samples


# ---------------------------------------------------------------------------
# Ground truth
# ---------------------------------------------------------------------------


def test_labels_record_what_was_actually_wrong():
    run = tep.Simulation(tep.Scenario.fault(1, hours=4)).run()
    labels = run.labels()

    assert labels["active"].shape == (len(run), tep.DISTURBANCES)
    assert labels["active"].dtype == np.bool_
    assert labels["active"][:, 0].all()
    # Elapsed time since onset, which for a fault present from the start is
    # simply the sample time.
    assert labels["since_onset"][:, 0] == pytest.approx(run.hours)
    # NaN, not zero, for a disturbance that never came on.
    assert np.isnan(labels["since_onset"][:, 1]).all()


def test_the_driver_forces_idv12_at_eight_hours():
    """`temain_mod.f:366-368`: delta D-011, and it is in every published file."""
    run = tep.Simulation(
        tep.Scenario.baseline(hours=10, driver_forces_idv12=True)
    ).run()
    active = run.labels()["active"][:, 11]

    # By step rather than by `hours`: the recorded time is the time at the start
    # of the step, so the sample at step 28800 reads 7.99972 h and would fall on
    # the wrong side of an 8.0 hour comparison.
    before = run.steps < tep.FORCED_DISTURBANCE_STEP
    assert not active[before].any()
    assert active[~before].all()


def test_idv12_forcing_is_off_by_default():
    assert not tep.Scenario.baseline().driver_forces_idv12
    run = tep.Simulation(tep.Scenario.baseline(hours=10)).run()

    assert not run.labels()["active"].any()


# ---------------------------------------------------------------------------
# The fault table
# ---------------------------------------------------------------------------


def test_the_fault_table():
    table = tep.faults()

    assert len(table) == tep.DISTURBANCES == 20
    assert [f.index for f in table] == list(range(1, 21))
    assert {f.shape for f in table} == {"step", "random", "sticking"}
    for fault in table:
        assert fault.published
        assert fault.effect
        assert fault.line.startswith("teprob.f:")


def test_sticking_faults_do_not_reach_the_plant():
    """A sticking fault widens a valve dead band; it touches no equation."""
    sticking = [f for f in tep.faults() if f.shape == "sticking"]

    assert [f.index for f in sticking] == [14, 15, 19]
    for fault in sticking:
        assert not fault.affects_the_plant
        assert fault.valves
        assert not fault.channels
    assert all(f.affects_the_plant for f in tep.faults() if f.shape != "sticking")


def test_the_five_unknown_faults_say_what_they_do():
    """The header calls IDV(16..20) unknown; the source is explicit."""
    unknown = [f for f in tep.faults() if f.published == "Unknown"]

    assert [f.index for f in unknown] == [16, 17, 18, 19, 20]
    for fault in unknown:
        assert fault.effect != "Unknown"
    # Three of the four non-sticking ones are spike trains, which is why the
    # literature reports them as the hardest to detect.
    assert [f.index for f in tep.faults() if f.spiking] == [17, 18, 20]


# ---------------------------------------------------------------------------
# Scenario
# ---------------------------------------------------------------------------


def test_scenario_defaults_are_the_baseline():
    assert tep.Scenario() == tep.Scenario.baseline()
    assert tep.Scenario(faults=[1]) == tep.Scenario.fault(1)

    baseline = tep.Scenario.baseline()
    assert baseline.seed == tep.DEFAULT_SEED
    assert baseline.hours == 48.0
    assert baseline.step_hours == tep.DEFAULT_STEP_HOURS
    assert baseline.sample_every == tep.DEFAULT_SAMPLE_EVERY
    assert baseline.faults == ()
    assert baseline.controlled
    assert baseline.steps == 172_800
    assert baseline.samples == 960


def test_scenario_builders():
    scenario = (
        tep.Scenario.baseline()
        .with_seed(7)
        .with_hours(3)
        .with_fault(4)
        .with_fault(11)
        .sampling_every(60)
        .open_loop()
    )

    assert scenario.seed == 7.0
    assert scenario.hours == 3.0
    assert scenario.faults == (4, 11)
    assert scenario.sample_every == 60
    assert not scenario.controlled
    assert scenario.samples == 3 * 60

    # Each one returns a new scenario rather than mutating in place.
    assert tep.Scenario.baseline().faults == ()


def test_scenario_repr_round_trips():
    Scenario = tep.Scenario  # noqa: N806  (the name `eval` needs in scope)
    for scenario in (
        Scenario.baseline(),
        Scenario.fault(3, hours=1, seed=5),
        Scenario(faults=[2, 7, 20], controlled=False, trip_ends_the_run=False),
    ):
        # These are all expressible as a constructor call, so that is what
        # `repr` should give: it is what is readable at a prompt.
        assert repr(scenario).startswith("Scenario(")
        assert eval(repr(scenario)) == scenario


def test_repr_round_trips_what_the_constructor_cannot_express():
    """`Scenario(...)` has no argument for a schedule, an extension or an
    integrator, so `repr` falls back to `from_text` for those.

    Without this, `repr` produced valid Python that evaluated to a *different*
    scenario with a different digest, and the class docstring promised it did
    not. Found from a notebook, not from a test, which is why this exists.
    """
    Scenario = tep.Scenario  # noqa: N806  (the name `eval` needs in scope)
    base = Scenario.baseline().to_text()
    for scenario in (
        # A scheduled fault that arrives at hour 6 and clears at hour 12.
        Scenario.from_text(base.replace("events=", "events=6:start:4,12:stop:4")),
        # A non-default integrator.
        Scenario.from_text(base.replace("integrator=euler", "integrator=rk4")),
        # The continuous-disturbance extension.
        Scenario.from_text(base.replace("continuous=0", "continuous=1")),
    ):
        assert repr(scenario).startswith("Scenario.from_text(")
        restored = eval(repr(scenario))
        assert restored == scenario
        # The digest is the thing that actually decides whether two scenarios
        # are the same run, so assert on it as well as on equality.
        assert restored.digest == scenario.digest


# ---------------------------------------------------------------------------
# Serialisation
# ---------------------------------------------------------------------------
#
# `repr` is for a human reading a traceback and only Python can read it back.
# `to_text` is the wire format: the same string the Rust, wasm and browser sides
# read and write, so a scenario built here can be run in a browser and a link
# out of a browser can be run here.


SCENARIOS = [
    tep.Scenario.baseline(),
    tep.Scenario.fault(4, hours=8),
    tep.Scenario(faults=[1, 6, 20], seed=42, hours=6.5, sample_every=60),
    tep.Scenario.baseline(controlled=False, driver_forces_idv12=True),
    tep.Scenario.baseline(trip_ends_the_run=False),
    tep.Scenario.baseline(step_hours=1 / 7200),
]


@pytest.mark.parametrize("scenario", SCENARIOS)
def test_a_scenario_round_trips_through_its_text(scenario):
    text = scenario.to_text()
    assert text.startswith("tepsim.scenario.v1;")
    back = tep.Scenario.from_text(text)
    assert back == scenario
    # The numbers survive exactly, not approximately. 4651207995 is the
    # generator word `teprob.f:1187` compiles in, and a text that rounded it
    # would name a run other than the one it produced.
    assert back.seed == scenario.seed
    assert back.step_hours == scenario.step_hours
    assert back.hours == scenario.hours
    # And the rendering is canonical: one scenario, one text.
    assert back.to_text() == text


@pytest.mark.parametrize("scenario", SCENARIOS)
def test_the_digest_survives_the_text(scenario):
    # The property that makes a serialised scenario worth anything. A dataset
    # labelled with this digest and shipped with the text can be checked.
    assert tep.Scenario.from_text(scenario.to_text()).digest == scenario.digest
    assert len(scenario.digest) == 16


def test_distinct_scenarios_have_distinct_digests_and_texts():
    seen_text = {s.to_text() for s in SCENARIOS}
    seen_digest = {s.digest for s in SCENARIOS}
    assert len(seen_text) == len(SCENARIOS)
    assert len(seen_digest) == len(SCENARIOS)


@pytest.mark.parametrize(
    ("text", "expected"),
    [
        # A field this build does not have. Rejected by name rather than
        # ignored, because ignoring it means running a different scenario.
        (tep.Scenario.baseline().to_text() + ";wobble=3", "wobble"),
        # A version this build does not read.
        (tep.Scenario.baseline().to_text().replace("v1", "v2", 1), "tepsim.scenario.v2"),
        # A field left out. Not defaulted: a text says what it runs.
        (tep.Scenario.baseline().to_text().replace(";trip=1", ""), "trip"),
        # A number that is not one.
        (tep.Scenario.baseline().to_text().replace("hours=48", "hours=soon"), "soon"),
        # A number out of range.
        (tep.Scenario.baseline().to_text().replace("seed=4651207995", "seed=0"), "seed"),
        # An integrator that does not exist.
        (
            tep.Scenario.baseline().to_text().replace("integrator=euler", "integrator=leapfrog"),
            "leapfrog",
        ),
        # Not a scenario text at all.
        ("", "tepsim.scenario.v1"),
        ("hours=1", "tepsim.scenario.v1"),
    ],
)
def test_a_bad_text_raises_and_says_what_was_wrong(text, expected):
    with pytest.raises(ValueError, match=expected):
        tep.Scenario.from_text(text)


def test_a_text_scenario_runs_to_the_same_numbers():
    # The end of the claim: a scenario that travelled as text produces the same
    # arrays as the one it came from.
    scenario = tep.Scenario.fault(1, hours=2)
    direct = tep.Simulation(scenario).run().to_numpy()
    travelled = tep.Simulation(tep.Scenario.from_text(scenario.to_text())).run().to_numpy()
    np.testing.assert_array_equal(direct, travelled)


@pytest.mark.parametrize(
    "make",
    [
        lambda: tep.Scenario.fault(0),
        lambda: tep.Scenario.fault(21),
        lambda: tep.Scenario(faults=[1, 99]),
        lambda: tep.Scenario(hours=-1),
        lambda: tep.Scenario(hours=float("nan")),
        lambda: tep.Scenario(step_hours=0),
        lambda: tep.Scenario(sample_every=0),
        lambda: tep.Scenario.baseline().with_fault(0),
        lambda: tep.Scenario.baseline().with_hours(-1),
        lambda: tep.Scenario.baseline().sampling_every(0),
    ],
)
def test_scenario_rejects_nonsense(make):
    with pytest.raises(ValueError):
        make()


@pytest.mark.parametrize(
    "access",
    [
        lambda run: run.measurement(0),
        lambda run: run.measurement(42),
        lambda run: run.manipulated(0),
        lambda run: run.manipulated(13),
        lambda run: run.column(53),
    ],
)
def test_channel_indices_are_checked(access):
    run = tep.Simulation(SHORT).run()
    with pytest.raises(ValueError):
        access(run)


def test_simulation_carries_its_scenario():
    scenario = tep.Scenario.fault(6, hours=1)
    simulation = tep.Simulation(scenario)

    assert simulation.scenario == scenario
    assert simulation.run().scenario == scenario
    assert tep.Simulation().scenario == tep.Scenario.baseline()


# ---------------------------------------------------------------------------
# Threads
# ---------------------------------------------------------------------------


@pytest.mark.skipif(
    (os.cpu_count() or 1) < 4, reason="needs four cores to show any parallelism"
)
def test_run_releases_the_gil():
    """An ensemble has to parallelise with threads, not just interleave.

    The margin is deliberately loose. Four independent 12-hour runs on four
    cores measured 5.4x on the machine this was written on; anything above 2x
    is impossible without the GIL actually being released, and anything below
    it is a scheduling accident rather than a regression.
    """
    workers = 4
    sims = [
        tep.Simulation(tep.Scenario.baseline(hours=12, seed=42 + 2 * i))
        for i in range(workers)
    ]

    start = time.perf_counter()
    for simulation in sims:
        simulation.run()
    serial = time.perf_counter() - start

    start = time.perf_counter()
    with ThreadPoolExecutor(max_workers=workers) as pool:
        threaded_runs = list(pool.map(tep.Simulation.run, sims))
    threaded = time.perf_counter() - start

    assert serial / threaded > 2.0, f"serial {serial:.3f} s, threaded {threaded:.3f} s"
    # Parallel or not, the answers have to be the ones the scenarios determine.
    for simulation, run in zip(sims, threaded_runs):
        assert np.array_equal(run.to_numpy(), simulation.run().to_numpy())
