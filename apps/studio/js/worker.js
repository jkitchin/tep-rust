// TEP Studio's simulation worker: one `Sim`, advanced in chunks, forever.
//
// The whole reason this file exists on a second thread is that the plant is a
// tight numeric loop and the user interface is not. A 48-hour run is 172,800
// Euler steps. On the main thread that is a frozen tab; here it is a thread
// that nothing is waiting on.
//
// Three details carry the design, and each of them is a decision rather than a
// convenience.
//
// Transfer, not clone. `postMessage(msg, [values.buffer])` moves the chunk's
// `ArrayBuffer` across threads as a pointer. Without the transfer list the
// structured clone algorithm copies every byte a second time, and at 54 values
// a row that is the dominant cost of a fast run.
//
// Yield between chunks. A worker in a `while (true)` loop is exactly as
// unresponsive as a blocked main thread, it just fails somewhere less visible.
// Returning to the event loop between chunks is what lets "stop", "pace" and
// "toggle IDV(6)" be heard at all. The chunk size is therefore a latency
// control: one sample costs `sampleEvery` integrator steps, so 20 samples at
// the default cadence is 3,600 steps of arithmetic before this thread can hear
// anything.
//
// No `SharedArrayBuffer`. It would remove even the pointer move, and it needs
// COOP and COEP response headers. Neither GitHub Pages nor a Hugging Face
// Static Space can set those, and free static hosting is the reason TEP Studio
// exists. See PLAN.org, "The browser application".

import init, { Scenario, Sim } from "./tepsim_wasm.js";

// The wasm-bindgen `--target web` glue fetches the `.wasm` beside itself, which
// is right in a browser and impossible in Node, whose `fetch` refuses `file:`
// URLs. `apps/studio/node/` drives this exact file under a shim of the Worker
// globals and sets the bytes here first, so the protocol below is tested rather
// than assumed. In a browser the value is `undefined` and `init` does what it
// always does.
const wasmOverride = globalThis.__tepsimWasmModule;
await init(wasmOverride === undefined ? undefined : { module_or_path: wasmOverride });

let sim = null;
let running = false;
let chunkSamples = 20;

// Wall-clock pacing. Zero means "as fast as this machine goes"; any other value
// is a multiple of real time to hold the run down to, so a human can watch the
// flowsheet move. This gates *when* the next chunk is asked for and never what
// it computes, so a paced run and an unpaced one emit identical bytes. The
// determinism check in `apps/studio/node/protocol.test.mjs` asserts exactly
// that, because a clock anywhere near a simulation deserves a test.
let speedMultiple = 0;
let runStartedAt = 0;

// Yield to the event loop so queued messages run. `setTimeout(0)` is clamped to
// about 4 ms after a few nested calls, which would cap throughput at 250 chunks
// a second; a MessageChannel round trip is not clamped.
const channel = new MessageChannel();
const pending = [];
channel.port1.onmessage = () => pending.shift()?.();
function yieldToEventLoop() {
  return new Promise((resolve) => {
    pending.push(resolve);
    channel.port2.postMessage(null);
  });
}

function now() {
  // `performance` exists in workers and in Node 16 and later. The fallback is
  // for the shim, and never reaches the simulation either way.
  return globalThis.performance ? globalThis.performance.now() : 0;
}

// Sleep until the wall clock has caught up with simulated time. Only ever
// called between chunks.
function sleep(ms) {
  return new Promise((resolve) => setTimeout(resolve, ms));
}

self.postMessage({ type: "ready" });

function report(error) {
  running = false;
  self.postMessage({ type: "error", message: String(error) });
}

self.onmessage = ({ data }) => {
  try {
    // `start` and `setFault` reach an async function, so a throw inside them
    // arrives as a rejected promise this `try` cannot see: the `Sim`
    // constructor rejecting a bad scenario would become an unhandled rejection
    // and the page would wait forever for a chunk. Both get their own catch.
    if (data.type === "start") {
      start(data).catch(report);
      return;
    }
    if (data.type === "stop") {
      running = false;
      return;
    }
    if (data.type === "setFault") {
      setFault(data).catch(report);
      return;
    }
    if (data.type === "speed") {
      // Takes effect at the next chunk boundary. Deliberately does not restart:
      // pacing is not part of the scenario and changes no number.
      speedMultiple = Number(data.multiple) || 0;
      return;
    }
    if (data.type === "chunkSamples") {
      chunkSamples = Math.max(1, Number(data.samples) | 0);
      return;
    }
  } catch (error) {
    report(error);
  }
};

// Every field a run's output depends on, copied across the thread boundary as
// plain numbers and booleans. Nothing here is optional: a scenario that left a
// field to a default on one side and set it on the other would produce a run
// whose digest did not match its link.
function buildScenario(request) {
  const scenario = new Scenario();
  scenario.seed = request.seed;
  scenario.hours = request.hours;
  scenario.stepHours = request.stepHours;
  scenario.sampleEvery = request.sampleEvery;
  scenario.controlled = request.controlled;
  scenario.driverForcesIdv12 = request.driverForcesIdv12;
  scenario.tripEndsTheRun = request.tripEndsTheRun;
  scenario.setIntegrator(request.integrator);
  scenario.clearFaults();
  for (const id of request.faults ?? []) scenario.setFault(id, true);
  return scenario;
}

// A disturbance cannot be switched on inside a run. `tepsim::Simulation` takes
// its scenario at construction and hands the disturbance vector to the driver
// there; the bindings' author found no seam in the native API and declined to
// invent one, because a run reachable only through the browser is a run nobody
// can reproduce, and reproducibility is the entire point of this project.
//
// So the toggle is recorded, `pendingRestart` goes true, and the run is rebuilt
// from step zero with the new scenario. At the throughput this thing actually
// gets, that is fast enough to feel immediate, and what comes out is a run a
// native caller can reproduce exactly.
async function setFault(data) {
  if (!sim) return;
  sim.setFault(data.id, data.active);
  if (!sim.pendingRestart) return;
  const scenario = sim.requestedScenario;
  self.postMessage({
    type: "restart",
    id: data.id,
    active: data.active,
    digest: scenario.digest(),
  });
  await startWith(scenario);
}

async function start(request) {
  running = false;
  // Let any in-flight loop see the flag and finish before `sim` is replaced.
  await yieldToEventLoop();
  if (request.chunkSamples !== undefined) {
    chunkSamples = Math.max(1, Number(request.chunkSamples) | 0);
  }
  if (request.speedMultiple !== undefined) {
    speedMultiple = Number(request.speedMultiple) || 0;
  }
  await startWith(buildScenario(request));
}

async function startWith(scenario) {
  running = false;
  await yieldToEventLoop();

  // Throws if the scenario cannot produce a well-defined run. The page sees it
  // as an "error" message rather than an unhandled rejection.
  const previous = sim;
  try {
    sim = new Sim(scenario);
  } catch (error) {
    // Free before rethrowing. A settings panel that fails validation on every
    // keystroke would otherwise leak a scenario handle per attempt.
    scenario.free();
    throw error;
  }

  self.postMessage({
    type: "started",
    rowWidth: sim.rowWidth,
    columnIds: sim.columnIds,
    columnLabels: sim.columnLabels,
    totalSamples: sim.totalSamples,
    totalSteps: sim.totalSteps,
    hours: scenario.hours,
    scenarioDigest: scenario.digest(),
    isFaithful: scenario.isFaithful,
    faults: [...scenario.activeFaults],
  });

  // Every wasm-bindgen object is a handle into wasm memory that JavaScript's
  // collector knows nothing about. The `Sim` constructor borrows the scenario
  // rather than taking it, so this side still owns it, and a page that restarts
  // on every fault toggle would otherwise grow without limit.
  scenario.free();
  previous?.free();

  running = true;
  runStartedAt = now();
  const mine = sim;

  while (running && sim === mine && !sim.isFinished) {
    const emittedBefore = sim.emittedSamples;
    const values = sim.stepChunk(chunkSamples);

    self.postMessage(
      {
        type: "chunk",
        values,
        emittedBefore,
        emitted: sim.emittedSamples,
        total: sim.totalSamples,
        hours: sim.hours,
        checksum: sim.checksum(),
        activeFaults: sim.activeFaults,
        finished: sim.isFinished,
        outcome: sim.outcome,
        tripHours: sim.tripHours,
        tripCause: sim.tripCause,
      },
      // The transfer list. Afterwards this thread's view of the buffer is
      // detached, which is correct: this thread is done with it.
      [values.buffer],
    );

    await yieldToEventLoop();

    if (speedMultiple > 0 && running && sim === mine) {
      // Simulated milliseconds the run has covered, divided by the multiple, is
      // where the wall clock should be. Ahead of that, wait.
      const target = (sim.hours * 3600_000) / speedMultiple;
      const behind = target - (now() - runStartedAt);
      if (behind > 1) await sleep(behind);
    }
  }

  if (sim === mine) running = false;
}
