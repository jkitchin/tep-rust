// Drive `js/worker.js` under Node, with no browser anywhere.
//
// The worker protocol is the part of TEP Studio most likely to break silently.
// A chart that draws the wrong colour is obvious; a restart that quietly
// carries the old fault set, or a paced run that emits different numbers from
// an unpaced one, is not. Neither can be caught by reading the file.
//
// So the Worker globals get a shim. `worker.js` needs three things from its
// environment: a `self` with `postMessage` and a settable `onmessage`, a
// `MessageChannel` (Node has had one since v15), and some way to load the wasm.
// That last one is why `worker.js` reads `globalThis.__tepsimWasmModule`: the
// wasm-bindgen `--target web` glue fetches the module beside itself, and Node's
// `fetch` refuses `file:` URLs. Handing it the bytes is a two-line seam and it
// is what makes the file below the same file the browser runs, rather than a
// copy of it that drifts.
//
// One honest limitation. The shim ignores the transfer list, so chunk buffers
// are not detached the way a real `postMessage` detaches them, and a bug that
// depended on that would not show up here. Everything else about the protocol
// does.

import { readFile } from "node:fs/promises";
import { fileURLToPath } from "node:url";
import path from "node:path";

const here = path.dirname(fileURLToPath(import.meta.url));
export const DIST = path.join(here, "..", "dist");

/** Read the built wasm module, with a useful message when it is not there. */
export async function wasmBytes() {
  const file = path.join(DIST, "js", "tepsim_wasm_bg.wasm");
  try {
    return await readFile(file);
  } catch {
    throw new Error(
      `${file} is missing. Build the app first: apps/studio/build.sh`,
    );
  }
}

/**
 * Install the shim and import the worker.
 *
 * ES modules are cached per process, so `worker.js` is evaluated once and its
 * `sim`, `running` and pacing state are shared by everything in the file that
 * calls this. That is exactly the situation in a browser, where there is also
 * one worker, so the tests are written to run in sequence against one of them.
 *
 * @returns {Promise<{send: Function, next: Function, drain: Function}>}
 */
export async function startWorker() {
  const bytes = await wasmBytes();
  globalThis.__tepsimWasmModule = bytes;

  /** @type {{type: string, [k: string]: unknown}[]} */
  const inbox = [];
  /** @type {{match: Function, resolve: Function}[]} */
  const waiters = [];

  // Take the message at `index` and discard everything queued before it.
  //
  // Ordered-stream semantics, which is what a page actually has: messages
  // arrive in order and waiting for a "started" means the chunks of the run it
  // replaced are gone. Without the discard a test that skipped past a stale
  // chunk would find it again later and assert against the wrong run, which is
  // exactly the false failure this harness produced before.
  const take = (index) => inbox.splice(0, index + 1).pop();

  const deliver = (message) => {
    inbox.push(message);
    // One waiter at a time, in order, so a match cannot jump the queue.
    const waiter = waiters[0];
    if (!waiter || !waiter.match(message)) return;
    waiters.shift();
    waiter.resolve(take(inbox.length - 1));
  };

  globalThis.self = {
    // The second argument is the transfer list, ignored here. See the note
    // above about what that costs.
    postMessage: (message) => deliver(message),
    onmessage: null,
  };

  await import(path.join(DIST, "js", "worker.js"));

  return {
    /** Post a message to the worker, as the page would. */
    send(message) {
      globalThis.self.onmessage({ data: message });
    },

    /**
     * Resolve with the next message satisfying `match`, searching messages
     * already delivered first so a fast worker cannot outrun the test.
     * Anything queued ahead of the match is discarded.
     *
     * @param {(m: object) => boolean} match
     * @param {number} timeoutMs
     */
    next(match, timeoutMs = 30_000) {
      const predicate = typeof match === "string" ? (m) => m.type === match : match;
      const index = inbox.findIndex(predicate);
      if (index >= 0) return Promise.resolve(take(index));
      return new Promise((resolve, reject) => {
        const timer = setTimeout(() => {
          reject(new Error(`timed out waiting for ${match}`));
        }, timeoutMs);
        waiters.push({
          match: predicate,
          resolve: (m) => {
            clearTimeout(timer);
            resolve(m);
          },
        });
      });
    },

    /** Throw away everything received so far. */
    drain() {
      inbox.length = 0;
    },

    /** Everything received so far, oldest first. */
    get received() {
      return inbox;
    },
  };
}

/**
 * The baseline request, with overrides. Mirrors what `app.js` posts.
 *
 * The scenario travels as the one string the bindings serialise it to, so this
 * builds a real `Scenario`, applies the overrides through its setters, and
 * sends `scenario.text`. Writing the text out by hand here would put a fifth
 * copy of the field list in the repository, which is the thing B-0054a removed.
 *
 * `chunkSamples` and `speedMultiple` are not part of a scenario: they decide
 * when the next chunk is asked for and never what is in it.
 */
export async function request(overrides = {}) {
  const { hours, faults, integrator, chunkSamples, speedMultiple, ...fields } = overrides;
  const { Scenario } = await bindings();
  const scenario = new Scenario();
  scenario.hours = hours ?? 1;
  if (integrator !== undefined) scenario.setIntegrator(integrator);
  for (const [key, value] of Object.entries(fields)) scenario[key] = value;
  for (const id of faults ?? []) scenario.setFault(id, true);
  const text = scenario.text;
  scenario.free();

  return {
    type: "start",
    scenario: text,
    chunkSamples: chunkSamples ?? 20,
    speedMultiple: speedMultiple ?? 0,
  };
}

/**
 * The wasm bindings, initialised once.
 *
 * Imported directly rather than through the worker, because building a request
 * needs a `Scenario` on this side of the boundary. ES modules are cached per
 * process, so `init` is the only thing that has to be guarded.
 */
let ready = null;
export function bindings() {
  ready ??= (async () => {
    const module = await import("../dist/js/tepsim_wasm.js");
    globalThis.__tepsimWasmModule ??= await wasmBytes();
    await module.default({ module_or_path: globalThis.__tepsimWasmModule });
    return module;
  })();
  return ready;
}

/**
 * Run one scenario to completion and report what came out.
 *
 * @param {{send: Function, next: Function, drain: Function}} worker
 * @param {object} overrides
 */
export async function runToCompletion(worker, overrides = {}) {
  worker.drain();
  worker.send(await request(overrides));
  const started = await worker.next("started");

  let chunks = 0;
  let rows = 0;
  let last = null;
  for (;;) {
    const chunk = await worker.next("chunk");
    chunks += 1;
    rows += chunk.values.length / started.rowWidth;
    last = chunk;
    if (chunk.finished) break;
  }
  return { started, last, chunks, rows };
}
