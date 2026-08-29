// The worker protocol, driven under Node.
//
// Run after a build:
//
//   apps/studio/build.sh
//   node --test --test-force-exit apps/studio/node/
//
// `--test-force-exit` is not papering over a leak. The worker's yield primitive
// is a `MessageChannel`, and an open `MessagePort` is a live handle that keeps
// Node's event loop from draining, so the process would sit there after the
// last assertion. In a browser the worker is torn down with the page and there
// is nothing to drain.
//
// These are the checks that would otherwise need a browser and a pair of eyes.
// The rendering cannot be tested here and is not pretended to be; what is here
// is every claim the page makes about numbers.

import test from "node:test";
import assert from "node:assert/strict";

import { runToCompletion, request, startWorker } from "./harness.mjs";

// Pinned by `crates/tepsim-wasm/tests/determinism.rs` from a native run of this
// commit: one hour of the baseline plant, 3,600 Euler steps.
const EXPECTED_SELF_CHECK = "c8a26889992f1719";

const worker = await startWorker();
await worker.next("ready");

test("the wasm module reproduces the native build bit for bit", async () => {
  // Import the bindings directly rather than through the worker: this is the
  // wasm half of Tier 9, and it should fail loudly and on its own line if a
  // toolchain change breaks it.
  const { default: init, selfCheckDigest, rowWidth, measurementCount, manipulatedCount } =
    await import("../dist/js/tepsim_wasm.js");
  await init({ module_or_path: globalThis.__tepsimWasmModule });

  assert.equal(
    selfCheckDigest(),
    EXPECTED_SELF_CHECK,
    "this JavaScript runtime does not reproduce the native reference",
  );
  // A packed row is [hours, XMEAS(1..41), XMV(1..12)]. The flowsheet and the
  // trend picker both index into it, so the width is load bearing.
  assert.equal(measurementCount(), 41);
  assert.equal(manipulatedCount(), 12);
  assert.equal(rowWidth(), 54);
});

test("a run streams chunks and finishes with the planned sample count", async () => {
  const { started, last, chunks, rows } = await runToCompletion(worker, { hours: 2 });

  assert.equal(started.rowWidth, 54);
  assert.equal(started.columnLabels.length, 54);
  assert.equal(started.isFaithful, true, "euler is the faithful integrator");

  // Two hours, sampled every 180 seconds, is 40 samples.
  assert.equal(started.totalSamples, 40);
  assert.equal(rows, 40, "every planned sample arrived exactly once");
  assert.equal(last.emitted, 40);
  assert.equal(last.outcome, "completed");
  assert.ok(chunks >= 2, `expected several chunks, got ${chunks}`);

  // The last sample's time is the run length to within one integrator step.
  assert.ok(Math.abs(last.hours - 2) < 1 / 3600 + 1e-9, `ended at ${last.hours} h`);
});

test("a run is a pure function of its scenario", async () => {
  const first = await runToCompletion(worker, { hours: 2 });
  const second = await runToCompletion(worker, { hours: 2 });
  assert.equal(
    first.last.checksum,
    second.last.checksum,
    "the same scenario produced different numbers on a second run",
  );
  assert.equal(first.started.scenarioDigest, second.started.scenarioDigest);
});

test("wall-clock pacing changes when chunks are asked for, never what is in them", async () => {
  // The one place a clock exists anywhere near the simulation. If pacing could
  // reach the numbers this would catch it, which is the whole reason the test
  // is here rather than in a comment.
  // The chunk sizes differ deliberately, so this also says that where the run
  // is cut into chunks changes nothing. Pacing only ever acts between chunks,
  // so a single-chunk run would not exercise it at all.
  const unpaced = await runToCompletion(worker, {
    hours: 0.5,
    speedMultiple: 0,
    chunkSamples: 20,
  });

  const before = Date.now();
  const paced = await runToCompletion(worker, {
    hours: 0.5,
    speedMultiple: 2000,
    chunkSamples: 2,
  });
  const elapsed = Date.now() - before;

  assert.equal(paced.last.checksum, unpaced.last.checksum);
  assert.equal(paced.last.emitted, unpaced.last.emitted);
  // 0.5 simulated hours at 2000x is 900 ms of wall clock, less the last chunk,
  // which is not followed by a wait. Allow generous slack: the point is that
  // pacing did something, not that it is precise.
  assert.ok(elapsed > 400, `paced run took only ${elapsed} ms, pacing did nothing`);
});

test("toggling a disturbance restarts the run and the new run carries it", async () => {
  worker.drain();
  // Long enough that the toggle certainly lands mid-run.
  worker.send(await request({ hours: 24 }));
  const started = await worker.next("started");
  assert.deepEqual([...started.faults], []);

  const before = await worker.next("chunk");
  assert.equal([...before.activeFaults].length, 0, "baseline has no ground truth fault");

  // IDV(6), A feed loss: a step fault that reaches the plant.
  worker.send({ type: "setFault", id: 6, active: true });

  const restart = await worker.next("restart");
  assert.equal(restart.id, 6);
  assert.equal(restart.active, true);

  const restarted = await worker.next("started");
  assert.deepEqual([...restarted.faults], [6]);
  assert.notEqual(
    restarted.scenarioDigest,
    started.scenarioDigest,
    "the digest must change when the scenario does, or a link is a lie",
  );

  const after = await worker.next("chunk");
  assert.equal(after.emittedBefore, 0, "the restarted run begins at step zero");
  assert.ok([...after.activeFaults].includes(6), "ground truth reports IDV(6)");

  worker.send({ type: "stop" });
});

test("a restarted run is identical to one started with the fault set", async () => {
  // The claim the fault panel rests on: rebuilding is not a browser-only path
  // that produces browser-only numbers.
  worker.drain();
  worker.send(await request({ hours: 2 }));
  await worker.next("started");
  await worker.next("chunk");
  worker.send({ type: "setFault", id: 8, active: true });
  await worker.next("restart");
  const restarted = await worker.next("started");

  let rebuilt = null;
  for (;;) {
    const chunk = await worker.next("chunk");
    if (chunk.finished) {
      rebuilt = chunk;
      break;
    }
  }

  const direct = await runToCompletion(worker, { hours: 2, faults: [8] });
  assert.equal(restarted.scenarioDigest, direct.started.scenarioDigest);
  assert.equal(rebuilt.checksum, direct.last.checksum);
});

test("an open-loop run trips, and the protocol says so rather than erroring", async () => {
  // A trip is a result, not a failure mode, so the page has to be able to show
  // one without treating it as an error and the protocol has to carry it.
  const { last, started } = await runToCompletion(worker, {
    hours: 6,
    controlled: false,
  });
  assert.equal(last.outcome, "tripped");
  assert.ok(last.tripHours > 0 && last.tripHours < 6, `tripped at ${last.tripHours} h`);
  assert.ok(typeof last.tripCause === "string" && last.tripCause.length > 0);
  // Delta D-007, signed off 2026-08-28: the run ends at the trip, so it is
  // shorter than the plan. Before that it froze the plant and emitted every
  // planned sample, and this assertion read `equal` instead.
  assert.ok(
    last.emitted < started.totalSamples,
    `the run continued past the trip: ${last.emitted} of ${started.totalSamples}`,
  );
});

test("and `tripEndsTheRun: false` freezes the plant and reports to the end", async () => {
  // The other half of D-007. `teprob.f:807-811` zeroes the derivatives and the
  // plant keeps reporting, which is what the frozen tails in the published
  // `d06` and `d18` files are, so the browser has to be able to reproduce it.
  const { last, started } = await runToCompletion(worker, {
    hours: 6,
    controlled: false,
    tripEndsTheRun: false,
  });
  assert.equal(last.outcome, "tripped");
  assert.equal(last.emitted, started.totalSamples);
});

test("an invalid scenario is reported, not thrown into the void", async () => {
  worker.drain();
  worker.send(await request({ hours: -1 }));
  const error = await worker.next("error");
  assert.match(error.message, /duration|positive|finite/i);
});

test("throughput, measured", async () => {
  // PLAN.org budgets at least 100x real time for the plant in wasm. This is
  // Node rather than a browser, but both are V8 on the same wasm engine, so it
  // is the right order of magnitude and it is a number rather than a hope.
  //
  // A large chunk is used deliberately: the yield between chunks is a latency
  // control, and measuring throughput at a 20-sample chunk measures the event
  // loop as much as the plant.
  const hours = 24;
  const before = process.hrtime.bigint();
  const { last } = await runToCompletion(worker, { hours, chunkSamples: 200 });
  const elapsedMs = Number(process.hrtime.bigint() - before) / 1e6;

  const multiple = (last.hours * 3600_000) / elapsedMs;
  const stepsPerSecond = (hours * 3600) / (elapsedMs / 1000);
  console.log(
    `    throughput: ${hours} simulated hours in ${elapsedMs.toFixed(0)} ms, ` +
      `${Math.round(multiple).toLocaleString()}x real time, ` +
      `${Math.round(stepsPerSecond).toLocaleString()} Euler steps per second`,
  );

  assert.ok(
    multiple >= 100,
    `PLAN.org budgets at least 100x real time; measured ${multiple.toFixed(0)}x`,
  );
});
