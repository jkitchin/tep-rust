// The simulation worker.
//
// One Sim per worker. The loop calls stepChunk, posts the Float64Array with its
// ArrayBuffer in the transfer list, then yields to the event loop before asking
// for the next chunk.
//
// Both halves of that matter.
//
// Transferring rather than cloning means a chunk crosses threads as a pointer
// move. postMessage otherwise structured-clones the buffer, which copies every
// byte a second time. SharedArrayBuffer would remove the copy entirely and is
// deliberately not used: it requires COOP and COEP response headers, which
// neither GitHub Pages nor a Hugging Face Static Space can set, and free static
// hosting is the reason the browser app exists. See PLAN.org, "The browser
// application".
//
// Yielding is what makes the worker answer messages at all. A worker in a tight
// loop is exactly as unresponsive as a blocked main thread; it just fails
// somewhere less visible. Between chunks, control returns to the event loop and
// a "setFault" or a "stop" gets processed. The chunk size is therefore a
// latency control: one sample costs sampleEvery integrator steps, so a chunk of
// 20 samples at the default cadence is 3,600 steps of work before this worker
// can hear anything.

import init, { Scenario, Sim } from "./pkg/tepsim_wasm.js";

let sim = null;
let running = false;

// Yield to the event loop so queued messages run. setTimeout(0) is clamped to
// about 4 ms after a few nested calls, which would cap throughput for small
// chunks; a MessageChannel round trip is not clamped.
const channel = new MessageChannel();
const pending = [];
channel.port1.onmessage = () => pending.shift()?.();
function yieldToEventLoop() {
  return new Promise((resolve) => {
    pending.push(resolve);
    channel.port2.postMessage(null);
  });
}

await init();
self.postMessage({ type: "ready" });

function report(error) {
  running = false;
  self.postMessage({ type: "error", message: String(error) });
}

self.onmessage = ({ data }) => {
  try {
    // `start` and `setFault` reach an async function, so a throw inside them
    // arrives as a rejected promise that this `try` cannot see: the Sim
    // constructor rejecting a bad scenario would become an unhandled rejection
    // and the caller would wait forever for a chunk. Both paths get their own
    // catch. Found by running this file under a shim of the Worker globals.
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
  } catch (error) {
    report(error);
  }
};

function buildScenario(request) {
  const scenario = new Scenario();
  scenario.seed = request.seed;
  scenario.hours = request.hours;
  scenario.stepHours = request.stepHours;
  scenario.sampleEvery = request.sampleEvery;
  scenario.controlled = request.controlled;
  scenario.setIntegrator(request.integrator);
  scenario.clearFaults();
  for (const id of request.faults ?? []) scenario.setFault(id, true);
  return scenario;
}

// A disturbance cannot be switched on inside a run: tepsim::Simulation takes
// its scenario at construction, and these bindings will not invent a side door
// the native API does not have, because a run reachable only through the
// browser is a run nobody can reproduce. So the request is recorded and the run
// is rebuilt from the start with the new scenario. At a hundred times real time
// that is fast enough to feel immediate, and the result is a run a native
// caller can reproduce exactly.
async function setFault(data) {
  if (!sim) return;
  sim.setFault(data.id, data.active);
  if (!sim.pendingRestart) return;
  const scenario = sim.requestedScenario;
  self.postMessage({ type: "restart", digest: scenario.digest() });
  await startWith(scenario);
}

async function start(request) {
  running = false;
  // Let any in-flight loop see the flag and finish before sim is replaced.
  await yieldToEventLoop();
  await startWith(buildScenario(request), request.chunkSamples);
}

let chunkSamples = 20;

async function startWith(scenario, requestedChunk) {
  running = false;
  await yieldToEventLoop();

  if (requestedChunk !== undefined) {
    chunkSamples = Math.max(1, requestedChunk | 0);
  }

  // Throws if the scenario cannot produce a well-defined run. The caller sees
  // it as an "error" message rather than an unhandled rejection; see the catch
  // in onmessage.
  const previous = sim;
  try {
    sim = new Sim(scenario);
  } catch (error) {
    // Free before rethrowing: a tab that fails validation on every keystroke
    // would otherwise leak a scenario handle per attempt.
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
    scenarioDigest: scenario.digest(),
    isFaithful: scenario.isFaithful,
  });

  // Every wasm-bindgen object is a handle into wasm memory that JavaScript's
  // collector knows nothing about. The Sim constructor borrows the scenario
  // rather than taking it, so this side still owns it, and a worker that
  // restarts on every fault toggle would otherwise grow without limit.
  scenario.free();
  previous?.free();

  running = true;
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
      // The transfer list. After this the worker's view of the buffer is
      // detached, which is correct: the worker is done with it.
      [values.buffer],
    );

    await yieldToEventLoop();
  }

  if (sim === mine) running = false;
}
