"""`book/src/python.md`'s quickstart, as a script that really runs.

`crates/tepsim/tests/book_examples.rs` pins the page's code fence to this file
byte for byte, and re-runs the file to check the transcript quoted beside it.

    .xtask-python/venv/bin/python book/examples/quickstart.py
"""

from concurrent.futures import ThreadPoolExecutor

import numpy as np

import tepsim as tep

print("tepsim %s: XMEAS(1..%d), XMV(1..%d), %d channels, IDV(1..%d)"
      % (tep.__version__, tep.MEASUREMENTS, tep.MANIPULATED, tep.CHANNELS,
         tep.DISTURBANCES))

# A Scenario says what to simulate, a Simulation does it, a Run holds what
# came out. That is the whole API.
run = tep.Simulation(tep.Scenario.baseline(seed=42, hours=48)).run()
print()
print("run:      %r" % run)
print("matrix:   %s %s" % (run.to_numpy().shape, run.to_numpy().dtype))
print("outcome:  %s" % run.outcome)
print("XMEAS(7): mean %.2f kPa over %.0f h"
      % (run.measurement(7).mean(), run.hours[-1]))

# The twenty disturbances say what they do, not only what the original header
# called them. Five of the published descriptions are the word "Unknown".
print()
print("the five the original leaves unexplained")
for fault in tep.faults():
    if fault.published == "Unknown":
        print("  IDV(%2d) %-9s %s" % (fault.index, fault.shape, fault.effect))

# Ground truth travels with the data, which the original records nowhere.
faulted = tep.Simulation(tep.Scenario.fault(1, hours=8)).run()
labels = faulted.labels()
print()
print("IDV(1) over 8 h: active at the last sample %s, %.2f h since onset"
      % (labels["active"][-1, 0], labels["since_onset"][-1, 0]))

# run() releases the GIL, so an ensemble is a thread pool and nothing has to be
# pickled.
sims = [tep.Simulation(tep.Scenario.fault(n, hours=8)) for n in range(1, 21)]
with ThreadPoolExecutor() as pool:
    runs = list(pool.map(tep.Simulation.run, sims))
print()
print("twenty 8-hour faulted runs: %d completed, %d tripped"
      % (sum(r.outcome == "completed" for r in runs),
         sum(r.outcome == "tripped" for r in runs)))

# A run is a pure function of its scenario, and a scenario is one line of text
# that parses back to an equal scenario.
scenario = tep.Scenario.fault(4, hours=8)
print()
print("digest:   %s" % scenario.digest)
print("text:     %s" % scenario.to_text())
print("parses back equal:      %s"
      % (tep.Scenario.from_text(scenario.to_text()) == scenario))
print("two runs bit-identical: %s"
      % np.array_equal(tep.Simulation(scenario).run().to_numpy(),
                       tep.Simulation(scenario).run().to_numpy()))
