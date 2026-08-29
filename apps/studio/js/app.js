// TEP Studio, main thread.
//
// This file owns the user interface and owns no arithmetic. The plant runs in
// `worker.js` on its own thread and sends back packed rows; everything here is
// forms, panels, a flowsheet and a repaint loop.
//
// Two rules keep it that way.
//
// The main thread never simulates. Not even one step, not even to preview a
// scenario. It builds a `Scenario` to ask the bindings what it would cost
// (`steps`, `sampleCount`, `validationError`, `digest`), which is arithmetic on
// four numbers, and everything else goes over `postMessage`.
//
// Repaints are driven by the display, not by the data. Chunks arrive far faster
// than a monitor refreshes when the run is unpaced. The message handler stores
// the chunk and sets a flag; a `requestAnimationFrame` loop does the drawing.
// Without that split a fast run spends all its time in `fillText`.
//
// # Determinism
//
// Nothing on this page may reach the simulation except the scenario. There is
// no `Math.random` here and no `Date`; `performance.now` appears twice, both
// times to divide simulated time by wall time for the throughput readout, which
// is a number *about* the run and never an input to it. The pacing control in
// the worker is the same: it decides when the next chunk is asked for and never
// what is in it.

import init, {
  Scenario,
  columnLabels,
  columnUnits,
  faults as faultTable,
  integrators,
  rowWidth as bindingRowWidth,
  sampledColumns,
  selfCheckDigest,
  version,
} from "./tepsim_wasm.js";

import { History } from "./history.js";
import { Recorder, csvParts, downloadName, exportMeta, saveCsv } from "./csv.js";
import { TrendGrid } from "./chart.js";
import { buildFlowsheet } from "./flowsheet.js";
import { decodeLink, encodeLink } from "./share.js";
import { formatValue } from "./format.js";

// Produced by a native run of the same commit; `crates/tepsim-wasm/tests/
// determinism.rs` pins it. One hour of the baseline plant, 3,600 Euler steps.
// If this browser disagrees, the cross-platform determinism invariant that the
// whole validation ladder rests on is broken here, and the page says so instead
// of drawing pretty lines over a lie.
const EXPECTED_SELF_CHECK = "c8a26889992f1719";

const $ = (id) => document.getElementById(id);

// ---------------------------------------------------------------------------
// Module load and the determinism check.
// ---------------------------------------------------------------------------

const banner = $("load-banner");
try {
  await init();
} catch (error) {
  banner.className = "banner bad";
  banner.textContent =
    `Could not load ./tepsim_wasm.js (${error}). Build the app first: ` +
    `apps/studio/build.sh, then serve apps/studio/dist over HTTP.`;
  throw error;
}

const digest = selfCheckDigest();
const deterministic = digest === EXPECTED_SELF_CHECK;
banner.className = deterministic ? "banner good" : "banner bad";
banner.innerHTML = deterministic
  ? `<strong>Determinism check passed.</strong> Self-check digest ` +
    `<span class="mono">${digest}</span> matches the native build: one hour of the ` +
    `baseline plant, 3,600 Euler steps, bit for bit. tepsim-wasm ${version()}.`
  : `<strong>Determinism check FAILED.</strong> This browser produced ` +
    `<span class="mono">${digest}</span> where the native build produces ` +
    `<span class="mono">${EXPECTED_SELF_CHECK}</span>. It is not reproducing the ` +
    `reference bit for bit. Do not trust any number below.`;

// ---------------------------------------------------------------------------
// Column metadata and defaults, taken from the bindings rather than retyped.
// ---------------------------------------------------------------------------

const LABELS = columnLabels();
const UNITS = columnUnits();
const ROW_WIDTH = bindingRowWidth();
const STEPPED = sampledColumns();

// The one place a `Scenario` becomes the plain object the form binds to.
// `pushStateIntoScenario` is its inverse, and the two are meant to be read
// together. This is the only field list left on this page: a form has one input
// per field and there is no way around that. What used to be beside it, and is
// now gone, is a second list in `startRequest`, a third in `worker.js` and a
// fourth in `share.js`, none of which any test could keep in step.
function scenarioFields(handle) {
  return {
    seed: handle.seed,
    hours: handle.hours,
    stepHours: handle.stepHours,
    sampleEvery: handle.sampleEvery,
    integrator: handle.integrator,
    controlled: handle.controlled,
    driverForcesIdv12: handle.driverForcesIdv12,
    tripEndsTheRun: handle.tripEndsTheRun,
    faults: [...handle.activeFaults],
  };
}

// A fresh `Scenario` is the definition of the baseline, so the page's defaults
// are read off it instead of being written down a second time and drifting.
const probe = new Scenario();
const DEFAULTS = {
  ...scenarioFields(probe),
  // Reactor pressure and temperature, both separator and stripper level, the
  // product rate and the product G composition: the six a first look at a
  // disturbance actually wants.
  channels: [7, 9, 12, 15, 17, 40],
  chunkSamples: 20,
  // Zero means "as fast as this machine goes", which is the headline this app
  // is here to demonstrate. The select offers real-time multiples for watching
  // the flowsheet.
  speedMultiple: 0,
};
// What a link is written against: the scenario as one serialised token, and the
// view fields, which reach no arithmetic and are therefore still spelled out.
const LINK_DEFAULTS = {
  scenario: probe.text,
  channels: DEFAULTS.channels,
  chunkSamples: DEFAULTS.chunkSamples,
  speedMultiple: DEFAULTS.speedMultiple,
};
probe.free();

let state = { ...DEFAULTS };

// Apply the scenario half of a link. The bindings parse it and they are strict,
// so a link from a build with a field this one does not have is refused by name
// rather than quietly opening a different run. A refusal is reported and the
// baseline is kept, because a mistyped link should still give a working page.
function applyLink(link) {
  state = {
    ...state,
    channels: [...link.channels],
    chunkSamples: link.chunkSamples,
    speedMultiple: link.speedMultiple,
  };
  if (link.scenario === LINK_DEFAULTS.scenario) return null;
  let parsed = null;
  try {
    parsed = Scenario.fromText(link.scenario);
  } catch (error) {
    return String(error);
  }
  Object.assign(state, scenarioFields(parsed));
  parsed.free();
  return null;
}

let linkProblem = applyLink(decodeLink(globalThis.location.hash, LINK_DEFAULTS, ROW_WIDTH));

// ---------------------------------------------------------------------------
// The scenario the panel is editing. One long-lived wasm handle, mutated in
// place, so a keystroke costs a few stores rather than an allocation.
// ---------------------------------------------------------------------------

const scenario = new Scenario();

function pushStateIntoScenario() {
  scenario.seed = state.seed;
  scenario.hours = state.hours;
  scenario.stepHours = state.stepHours;
  scenario.sampleEvery = state.sampleEvery;
  scenario.controlled = state.controlled;
  scenario.driverForcesIdv12 = state.driverForcesIdv12;
  scenario.tripEndsTheRun = state.tripEndsTheRun;
  try {
    scenario.setIntegrator(state.integrator);
  } catch {
    // `state.integrator` comes either from the select, which offers only names
    // the bindings gave it, or from a `Scenario` the bindings parsed, so this
    // is unreachable. Swallowing it stops a hostile link from killing the page.
    state.integrator = DEFAULTS.integrator;
    scenario.setIntegrator(DEFAULTS.integrator);
  }
  scenario.clearFaults();
  for (const id of state.faults) scenario.setFault(id, true);
}

// ---------------------------------------------------------------------------
// Panels.
// ---------------------------------------------------------------------------

const trends = new TrendGrid($("trends"));
trends.setColumns(LABELS, UNITS, STEPPED);

const flowsheet = buildFlowsheet($("pfd"), LABELS, UNITS);

// The twenty disturbances, each with its published description. Toggling one
// restarts the run, because a simulation takes its scenario at construction and
// the bindings decline to invent a side door the native API does not have. The
// panel says so rather than pretending the toggle is live.
{
  const list = faultTable();
  const table = $("faults");
  table.innerHTML = list
    .map(
      (f) => `<tr>
        <td><input type="checkbox" id="idv-${f.id}" data-fault="${f.id}" /></td>
        <td class="mono nowrap"><label for="idv-${f.id}">IDV(${f.id})</label></td>
        <td>${f.label}</td>
        <td class="muted nowrap">${f.kind}${f.affectsThePlant ? "" : ", valves only"}</td>
      </tr>`,
    )
    .join("");
  // Every `Fault` is a handle into wasm memory that JavaScript's collector
  // knows nothing about. The markup above has copied every string it needs.
  for (const f of list) f.free();
}

// The integrator choices come from the bindings, so the select cannot offer a
// name the module would reject.
$("integrator").innerHTML = integrators()
  .map((name) => `<option value="${name}">${name}</option>`)
  .join("");

// The trend picker: every channel except the time column.
{
  const picker = $("channel-picker");
  picker.innerHTML = LABELS.slice(1)
    .map((label, i) => {
      const col = i + 1;
      return `<label class="pick"><input type="checkbox" data-channel="${col}" />
        <span>${label}</span></label>`;
    })
    .join("");
}

// ---------------------------------------------------------------------------
// Reading the form into `state`, and writing `state` back into the form.
// ---------------------------------------------------------------------------

const NUMERIC_FIELDS = ["seed", "hours", "stepHours", "sampleEvery", "chunkSamples"];
const BOOLEAN_FIELDS = ["controlled", "driverForcesIdv12", "tripEndsTheRun"];

// Assigning a `<select>` a value it does not offer leaves it blank, and the
// next `readForm` would then write that blank back into the state. A link can
// legitimately carry such a value: the bindings accept `dormand-prince` as an
// alias for `dopri5`, and the pacing menu offers five multiples out of a
// continuum. So an unofferable value is normalised to the fallback here, once,
// where the state can be corrected with it.
function setSelect(id, value, fallback) {
  const select = $(id);
  select.value = String(value);
  if (select.value === String(value)) return value;
  select.value = String(fallback);
  return fallback;
}

function writeForm() {
  for (const field of NUMERIC_FIELDS) $(field).value = String(state[field]);
  for (const field of BOOLEAN_FIELDS) $(field).checked = state[field];
  state.integrator = setSelect("integrator", state.integrator, DEFAULTS.integrator);
  state.speedMultiple = setSelect(
    "speedMultiple",
    state.speedMultiple,
    DEFAULTS.speedMultiple,
  );
  for (const box of document.querySelectorAll("[data-fault]")) {
    box.checked = state.faults.includes(Number(box.dataset.fault));
  }
  for (const box of document.querySelectorAll("[data-channel]")) {
    box.checked = state.channels.includes(Number(box.dataset.channel));
  }
}

function readForm() {
  for (const field of NUMERIC_FIELDS) {
    const value = Number($(field).value);
    if (Number.isFinite(value)) state[field] = value;
  }
  state.sampleEvery = Math.max(1, Math.round(state.sampleEvery));
  state.chunkSamples = Math.max(1, Math.round(state.chunkSamples));
  for (const field of BOOLEAN_FIELDS) state[field] = $(field).checked;
  state.integrator = $("integrator").value;
  state.speedMultiple = Number($("speedMultiple").value) || 0;
  state.faults = [...document.querySelectorAll("[data-fault]")]
    .filter((b) => b.checked)
    .map((b) => Number(b.dataset.fault));
  state.channels = [...document.querySelectorAll("[data-channel]")]
    .filter((b) => b.checked)
    .map((b) => Number(b.dataset.channel));
}

function syncFragment() {
  // The link carries the scenario as the bindings serialise it, so the handle
  // has to be current first. `refreshScenarioReadout` pushes as well; both are
  // a few dozen stores and neither allocates.
  pushStateIntoScenario();
  const fragment = encodeLink(
    {
      scenario: scenario.text,
      channels: state.channels,
      chunkSamples: state.chunkSamples,
      speedMultiple: state.speedMultiple,
    },
    LINK_DEFAULTS,
  );
  // `replaceState` rather than assigning `location.hash`: assigning pushes a
  // history entry per keystroke and turns the back button into an undo log for
  // a spinner.
  globalThis.history.replaceState(null, "", fragment === "" ? "#" : `#${fragment}`);
  $("share").value = `${globalThis.location.origin}${globalThis.location.pathname}#${fragment}`;
}

function refreshScenarioReadout() {
  pushStateIntoScenario();
  const error = scenario.validationError();
  // A link this build cannot honour is reported and takes the message slot,
  // because it is the news: the panel is showing the baseline, and without
  // saying so the page would look as though that is what the link asked for.
  // It does not block the run. The panel itself is valid, and a page that
  // refused to start because of a link the user has already been told about
  // would be a dead end.
  const message = linkProblem ?? error;
  $("scenario-error").textContent = message ?? "";
  $("scenario-error").hidden = !message;
  // `setButtons` owns both buttons, so that editing a field mid-run cannot
  // re-enable Start behind the running flag's back.
  setButtons();
  $("plan").textContent = error
    ? "-"
    : `${scenario.steps.toLocaleString()} steps, ` +
      `${scenario.sampleCount.toLocaleString()} samples`;
  $("stat-digest").textContent = error ? "-" : scenario.digest();

  const warn = $("integrator-banner");
  warn.hidden = scenario.isFaithful;
  if (!scenario.isFaithful) {
    warn.className = "banner warn";
    warn.innerHTML =
      `<strong>${scenario.integrator} is not the original.</strong> Only Euler ` +
      `reproduces the Fortran, and every claim the validation ladder makes is a claim ` +
      `about Euler. This is a better integration of the same equations and a different ` +
      `set of numbers.`;
  }
  return !error;
}

// ---------------------------------------------------------------------------
// The run.
// ---------------------------------------------------------------------------

let history = new History(ROW_WIDTH);
// The same chunks, kept whole. `History` decimates so the charts stay cheap;
// this keeps every sample so the download is the run rather than a picture of
// it. It holds the chunk arrays the worker transferred, so the second store
// costs references and not copies. See `js/csv.js`.
let recorder = new Recorder(ROW_WIDTH);
// The "started" message of the run being recorded, which is what a saved file
// is labelled with. Not `state`: the panel can be edited mid-run, and then the
// panel is describing a run that is not this one.
let run = null;
let worker = null;
let workerReady = false;
let queuedStart = false;
let dirty = false;
let running = false;
let startedAt = 0;
let chunks = 0;
let lastMessage = null;

function statsReset() {
  chunks = 0;
  lastMessage = null;
  $("stat-samples").textContent = "0";
  $("stat-hours").textContent = "0";
  $("stat-wall").textContent = "0 ms";
  $("stat-speed").textContent = "-";
  $("stat-checksum").textContent = "-";
  $("stat-outcome").textContent = "-";
  $("stat-truth").textContent = "-";
  $("stat-recorded").textContent = "-";
  // Nothing has been recorded of this run yet, so there is nothing to save.
  $("download").disabled = true;
}

function ensureWorker() {
  if (worker) return worker;
  // `new Worker` resolves a bare string against the *document* base URL, not
  // against this module, so a relative path would look for the worker beside
  // index.html and fail. `import.meta.url` is the module's own address and is
  // what makes the layout below the page's own business.
  worker = new Worker(new URL("./worker.js", import.meta.url), { type: "module" });
  worker.onmessage = ({ data }) => onWorkerMessage(data);
  worker.onerror = (event) => {
    $("scenario-error").hidden = false;
    $("scenario-error").textContent = `worker failed: ${event.message ?? event}`;
    running = false;
    setButtons();
  };
  return worker;
}

// The scenario crosses to the worker as the one string the bindings serialise
// it to, and nothing else about it does. `chunkSamples` and `speedMultiple` are
// not part of the scenario: they decide when the next chunk is asked for and
// never what is in it, so a paced run and an unpaced one emit identical bytes.
function startRequest() {
  // The caller has just run `refreshScenarioReadout`, which pushed the form
  // into the handle, so `text` is the scenario the panel is showing.
  return {
    type: "start",
    scenario: scenario.text,
    chunkSamples: state.chunkSamples,
    speedMultiple: state.speedMultiple,
  };
}

function start() {
  if (!refreshScenarioReadout()) return;
  ensureWorker();
  running = true;
  setButtons();
  if (!workerReady) {
    // The module is still instantiating on the worker thread. Remember the
    // intent and send it from the "ready" handler.
    queuedStart = true;
    return;
  }
  worker.postMessage(startRequest());
}

function stop() {
  running = false;
  queuedStart = false;
  worker?.postMessage({ type: "stop" });
  setButtons();
}

function setButtons() {
  $("start").disabled = running || Boolean(scenario.validationError());
  $("stop").disabled = !running;
}

function onWorkerMessage(data) {
  if (data.type === "ready") {
    workerReady = true;
    if (queuedStart) {
      queuedStart = false;
      worker.postMessage(startRequest());
    }
    return;
  }

  if (data.type === "started") {
    run = data;
    history = new History(data.rowWidth);
    recorder = new Recorder(data.rowWidth);
    trends.setHours(data.hours);
    trends.setSelection(state.channels);
    flowsheet.update(null);
    statsReset();
    startedAt = performance.now();
    running = true;
    dirty = true;
    setButtons();
    return;
  }

  if (data.type === "restart") {
    $("stat-outcome").textContent =
      `restarted from step zero: IDV(${data.id}) ${data.active ? "on" : "off"}`;
    return;
  }

  if (data.type === "error") {
    $("scenario-error").hidden = false;
    $("scenario-error").textContent = data.message;
    running = false;
    setButtons();
    return;
  }

  if (data.type === "chunk") {
    // `data.values` arrived by transfer: this thread now owns the buffer and
    // the worker's view of it is detached. `History` copies out of it and the
    // recorder keeps it, so the order of these two does not matter and neither
    // of them costs a second buffer.
    history.append(data.values);
    recorder.append(data.values, data.emittedBefore, data.activeFaults, data.faultAges);
    chunks += 1;
    lastMessage = data;
    dirty = true;
    if (data.finished) {
      running = false;
      setButtons();
    }
  }
}

// The repaint loop. One pass per displayed frame, no matter how many chunks
// landed in between.
function frame() {
  if (dirty) {
    dirty = false;
    trends.draw(history);
    flowsheet.update(history.lastRow());
    if (lastMessage) drawStats(lastMessage);
  }
  requestAnimationFrame(frame);
}

function drawStats(data) {
  const wall = performance.now() - startedAt;
  $("stat-samples").textContent = `${data.emitted.toLocaleString()} of ${data.total.toLocaleString()}`;
  $("stat-hours").textContent = formatValue(data.hours);
  $("stat-wall").textContent = `${wall.toFixed(0)} ms`;
  // Simulated milliseconds over wall milliseconds. Meaningless while pacing is
  // on, because pacing is a cap on exactly this number, so say so instead of
  // reporting the cap as a measurement.
  $("stat-speed").textContent =
    state.speedMultiple > 0
      ? `paced to ${state.speedMultiple}x`
      : wall > 0
        ? `${Math.round((data.hours * 3600_000) / wall).toLocaleString()}x real time`
        : "-";
  $("stat-checksum").textContent = data.checksum;

  const active = [...data.activeFaults];
  $("stat-truth").textContent = active.length
    ? `IDV ${active.join(", ")} active`
    : "no disturbance active";

  if (data.outcome === "tripped") {
    $("stat-outcome").textContent =
      `tripped at ${data.tripHours.toFixed(3)} h: ${data.tripCause}` +
      (data.finished ? "" : " (the plant freezes and keeps reporting, as teprob.f does)");
  } else if (data.finished) {
    $("stat-outcome").textContent = data.outcome ?? "finished";
  } else {
    $("stat-outcome").textContent = "running";
  }

  // Decimation is invisible unless the page admits to it.
  $("stat-retained").textContent =
    history.decimation === 1
      ? `${history.count.toLocaleString()} rows, full resolution`
      : `${history.count.toLocaleString()} rows, every ${history.decimation}th kept`;

  // And so is a record that stopped short. The two readouts differ on purpose:
  // "retained" is what the charts are drawing and "recorded" is what the
  // download would contain, and after 20,000 samples they are not the same
  // thing.
  $("stat-recorded").textContent = recorder.truncated
    ? `${recorder.count.toLocaleString()} rows, then the ceiling; ` +
      `${recorder.dropped.toLocaleString()} not recorded`
    : `${recorder.count.toLocaleString()} rows, full resolution`;
  $("download").disabled = recorder.count === 0;
}

// ---------------------------------------------------------------------------
// Wiring.
// ---------------------------------------------------------------------------

function onSettingChanged({ restart = false, fromUser = true } = {}) {
  // Touching the panel answers the complaint about the link: whatever the link
  // asked for, this is now what the user asked for.
  if (fromUser) linkProblem = null;
  readForm();
  syncFragment();
  refreshScenarioReadout();
  if (restart && running) start();
}

for (const field of NUMERIC_FIELDS) {
  $(field).addEventListener("input", () => onSettingChanged());
}
for (const field of [...BOOLEAN_FIELDS, "integrator"]) {
  $(field).addEventListener("change", () => onSettingChanged());
}

// Pacing is not part of the scenario and changes no number, so it takes effect
// without a restart.
$("speedMultiple").addEventListener("change", () => {
  readForm();
  syncFragment();
  worker?.postMessage({ type: "speed", multiple: state.speedMultiple });
});

// A fault toggle mid-run has to restart, and the worker is the one that knows
// whether a run is in flight, so the request goes there and comes back as a
// "restart" message. With no run in flight this only edits the scenario.
$("faults").addEventListener("change", (event) => {
  const id = Number(event.target.dataset?.fault);
  if (!id) return;
  readForm();
  syncFragment();
  refreshScenarioReadout();
  if (running) {
    worker.postMessage({ type: "setFault", id, active: event.target.checked });
  }
});

$("channel-picker").addEventListener("change", (event) => {
  if (!event.target.dataset?.channel) return;
  readForm();
  syncFragment();
  trends.setSelection(state.channels);
  dirty = true;
});

$("start").addEventListener("click", () => start());
$("stop").addEventListener("click", () => stop());
$("reset").addEventListener("click", () => {
  stop();
  state = { ...DEFAULTS, faults: [], channels: [...DEFAULTS.channels] };
  writeForm();
  onSettingChanged();
  start();
});
$("download").addEventListener("click", () => {
  // The button is disabled in both of these cases; the guard is here because a
  // disabled button is a display property and this is the invariant.
  if (!run || recorder.count === 0) return;
  const meta = exportMeta(run, lastMessage, recorder, version());
  saveCsv(csvParts(recorder, meta), downloadName(meta));
});

$("copy-share").addEventListener("click", async () => {
  const button = $("copy-share");
  try {
    await navigator.clipboard.writeText($("share").value);
    button.textContent = "copied";
  } catch {
    // Clipboard access is permission gated and refused outright in some
    // contexts. The field beside the button is selectable, so this is a
    // downgrade rather than a failure.
    $("share").select();
    button.textContent = "select and copy";
  }
  setTimeout(() => {
    button.textContent = "copy link";
  }, 1500);
});

globalThis.addEventListener("hashchange", () => {
  // Somebody pasted a different link into the same tab. Start from the
  // baseline, so a link that omits a view field means the default rather than
  // whatever the panel happened to be showing.
  state = { ...DEFAULTS };
  linkProblem = applyLink(decodeLink(globalThis.location.hash, LINK_DEFAULTS, ROW_WIDTH));
  writeForm();
  onSettingChanged({ fromUser: false });
  if (running) start();
});

globalThis.addEventListener("resize", () => {
  trends.resize();
  dirty = true;
});

writeForm();
onSettingChanged({ fromUser: false });
trends.setHours(state.hours);
trends.setSelection(state.channels);
requestAnimationFrame(frame);

// Start straight away. A simulator that opens on an empty chart and waits to be
// asked has buried its own point.
start();
