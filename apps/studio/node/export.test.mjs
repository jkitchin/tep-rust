// Saving a run as a file.
//
// A download is the one thing on this page that leaves it, so it is the one
// thing nobody can check afterwards by looking at the page. Everything that
// makes a file trustworthy is here: that it is the whole run and not the
// decimated view the charts draw, that a value survives the round trip to text
// exactly, that the step numbers are the ones `Sample::step` would give, that
// the ground truth is resolved to the sample rather than to the chunk it
// arrived in, and that the scenario in the header reads back as the run that
// produced the numbers under it.
//
// The last few drive the real worker, so what is exported is what a browser
// would export. Only `a.click()` is out of reach, and `saveCsv` takes its host
// as an argument so that everything up to and including revoking the object
// URL is exercised here rather than asserted in a comment.

import test from "node:test";
import assert from "node:assert/strict";

import {
  MAX_RECORDED_ROWS,
  Recorder,
  csvParts,
  downloadName,
  exportMeta,
  saveCsv,
} from "../dist/js/csv.js";
import { History } from "../dist/js/history.js";
import { bindings, request, startWorker } from "./harness.mjs";

const ROW_WIDTH = 54;

// Started before anything is registered. `node:test` begins running as soon as
// the first `test` call is made, and a `test` registered after a top-level
// await that came after one is never seen.
const worker = await startWorker();
await worker.next("ready");
const { version } = await bindings();

/** `count` packed rows whose time column counts up and whose column 1 is `n`. */
function rows(from, count, width = ROW_WIDTH) {
  const out = new Float64Array(count * width);
  for (let i = 0; i < count; i += 1) {
    out[i * width] = (from + i) / 20;
    out[i * width + 1] = from + i;
  }
  return out;
}

/** A plausible `meta`, with whatever the test cares about overridden. */
function meta(overrides = {}) {
  return {
    version: "0.1.0",
    columnIds: ["time_hours", ...Array.from({ length: 53 }, (_, i) => `ch${i + 1}`)],
    sampleEvery: 180,
    scenario: "tepsim.scenario.v1;seed=4651207995;hours=48",
    scenarioDigest: "1122334455667788",
    faults: [],
    planned: 960,
    checksum: "6ee4409dc3b1fe0f",
    outcome: "completed",
    finished: true,
    tripHours: null,
    tripCause: null,
    hours: 2,
    rows: 40,
    truncated: false,
    ...overrides,
  };
}

/** Split a rendered file into its comment block, its header and its rows. */
function parse(parts) {
  const lines = parts.join("").split("\n");
  assert.equal(lines.pop(), "", "the file must end in a newline");
  const comments = lines.filter((line) => line.startsWith("#"));
  const data = lines.filter((line) => !line.startsWith("#"));
  return {
    comments,
    columns: data[0].split(","),
    rows: data.slice(1).map((line) => line.split(",")),
  };
}

test("the recorder keeps every sample where the trend store decimates", () => {
  // The two stores see the same chunks and answer different questions. If this
  // ever stops being true the download quietly becomes a picture of the run.
  const recorder = new Recorder(ROW_WIDTH);
  const history = new History(ROW_WIDTH, 64);
  for (let c = 0; c < 50; c += 1) {
    const values = rows(c * 100, 100);
    recorder.append(values, c * 100, [], []);
    history.append(values);
  }

  assert.equal(recorder.count, 5000);
  assert.equal(recorder.truncated, false);
  assert.ok(history.decimation > 1, "the history store should have decimated");

  const file = parse(csvParts(recorder, meta()));
  assert.equal(file.rows.length, 5000);
  // Column 1 of the packed row counts the samples offered, so this says every
  // one of them is in the file, in order, exactly once.
  file.rows.forEach((row, i) => assert.equal(Number(row[2]), i));
});

test("the ceiling stops recording rather than thinning what it has", () => {
  const recorder = new Recorder(ROW_WIDTH, 25);
  for (let c = 0; c < 3; c += 1) recorder.append(rows(c * 10, 10), c * 10, [], []);

  assert.equal(recorder.count, 25);
  assert.equal(recorder.dropped, 5);
  assert.equal(recorder.truncated, true);

  const file = parse(csvParts(recorder, meta({ truncated: true, rows: 25 })));
  assert.equal(file.rows.length, 25);
  // A prefix of the run at full resolution, not a sample of all of it. That is
  // the whole point of the ceiling behaving this way.
  file.rows.forEach((row, i) => assert.equal(Number(row[2]), i));
  assert.ok(
    file.comments.some((line) => line.includes("Recording stopped at the 25 row ceiling")),
    `the header must own up to the ceiling: ${file.comments.join("\n")}`,
  );
});

test("the default ceiling covers the longest run the panel can ask for", () => {
  // 48 hours at `sampleEvery = 1` is 172,800 samples, and that is the run
  // somebody who wants a file at full cadence is going to ask for.
  assert.ok(
    MAX_RECORDED_ROWS >= 172_800,
    `${MAX_RECORDED_ROWS} rows would truncate a 48-hour run at full cadence`,
  );
});

test("the header is the bindings' column list, framed by the step and the labels", async () => {
  const { columnIds } = await bindings();
  const ids = columnIds();
  const recorder = new Recorder(ROW_WIDTH);
  recorder.append(rows(0, 2), 0, [], []);

  const file = parse(csvParts(recorder, meta({ columnIds: ids })));
  assert.equal(file.columns.length, ROW_WIDTH + 3);
  assert.equal(file.columns[0], "step");
  assert.deepEqual(file.columns.slice(1, 1 + ROW_WIDTH), ids);
  assert.equal(file.columns[1], "time_hours");
  assert.equal(file.columns[2], "XMEAS_1_A_feed");
  assert.equal(file.columns[43], "XMV_1_D_feed_flow");
  // `tepsim::recorder::Csv` names its ground-truth columns the same two things,
  // so a file from the browser and a file from a native run read alike.
  assert.deepEqual(file.columns.slice(-2), ["fault", "hours_since_onset"]);
  // Every row has a cell for every column, including the empty label cells.
  for (const row of file.rows) assert.equal(row.length, file.columns.length);
});

test("a value survives the round trip to text exactly", () => {
  // A deterministic simulator whose export is approximate is a simulator whose
  // export cannot be checked against it. `String` is the shortest text that
  // reads back as the same f64; the awkward cases are the ones that say so.
  const awkward = [
    0.1,
    1 / 3,
    2705.1999999999998,
    -273.15,
    Number.MIN_VALUE,
    Number.MAX_VALUE,
    Number.EPSILON,
    1e-300,
    2 ** 53 + 2,
  ];
  const values = new Float64Array(ROW_WIDTH);
  for (let c = 0; c < ROW_WIDTH; c += 1) values[c] = awkward[c % awkward.length];

  const recorder = new Recorder(ROW_WIDTH);
  recorder.append(values, 0, [], []);
  const file = parse(csvParts(recorder, meta()));

  for (let c = 0; c < ROW_WIDTH; c += 1) {
    assert.equal(
      Number(file.rows[0][c + 1]),
      values[c],
      `column ${c} did not survive: ${file.rows[0][c + 1]}`,
    );
  }
});

test("the step column is the one Sample::step would give", () => {
  // The packed row carries the time and not the step, so this is the single
  // reconstructed column in the file and the only one that can be wrong on its
  // own. A sample lands on every `sampleEvery`th step and steps count from one.
  const recorder = new Recorder(ROW_WIDTH);
  recorder.append(rows(0, 3), 0, [], []);
  recorder.append(rows(3, 2), 3, [], []);

  const file = parse(csvParts(recorder, meta({ sampleEvery: 180 })));
  assert.deepEqual(
    file.rows.map((row) => Number(row[0])),
    [180, 360, 540, 720, 900],
  );

  const every = parse(csvParts(recorder, meta({ sampleEvery: 1 })));
  assert.deepEqual(
    every.rows.map((row) => Number(row[0])),
    [1, 2, 3, 4, 5],
  );
});

test("ground truth is resolved to the row, not to the chunk it arrived in", () => {
  // The chunk's flags were read at its last sample. Handing them to every row
  // in the chunk would move the onset back to the start of the chunk, which is
  // the exact quantity a detection-delay figure measures.
  const recorder = new Recorder(ROW_WIDTH);
  // Ten rows at 0.05 h apart, ending at 0.45 h; IDV(6) had been on for 0.2 h by
  // then, so it came on at 0.25 h and the first six rows are clean.
  recorder.append(rows(0, 10), 0, [6], [0.2]);

  const file = parse(csvParts(recorder, meta()));
  const labels = file.rows.map((row) => row[ROW_WIDTH + 1]);
  assert.deepEqual(labels, ["", "", "", "", "", "6", "6", "6", "6", "6"]);

  const ages = file.rows.map((row) => row[ROW_WIDTH + 2]);
  assert.deepEqual(ages.slice(0, 5), ["", "", "", "", ""]);
  assert.equal(Number(ages[5]), 0, "the onset row's age must be exactly zero");
  assert.equal(Number(ages[9]), 0.2, "and the last row's age is what was reported");
});

test("two disturbances with different onsets are labelled separately", () => {
  const recorder = new Recorder(ROW_WIDTH);
  // IDV(1) from the start of the run, IDV(12) from 0.3 h.
  recorder.append(rows(0, 10), 0, [1, 12], [0.45, 0.15]);
  const file = parse(csvParts(recorder, meta()));
  assert.deepEqual(
    file.rows.map((row) => row[ROW_WIDTH + 1]),
    ["1", "1", "1", "1", "1", "1", "1", "1 12", "1 12", "1 12"],
  );
});

test("the file name carries the fault, the size, and something that changes", () => {
  const baseline = downloadName(meta({ hours: 48, rows: 960 }));
  assert.equal(baseline, "tep-baseline-48h-960rows-6ee4409dc3b1fe0f.csv");

  const faulted = downloadName(meta({ faults: [6], hours: 2.5, rows: 50 }));
  assert.match(faulted, /^tep-idv6-2\.5h-50rows-/);
  assert.equal(downloadName(meta({ faults: [1, 12] })).startsWith("tep-idv1+12-"), true);

  // A trip and a truncated record are both things somebody should see before
  // opening the file.
  assert.match(
    downloadName(meta({ faults: [6], outcome: "tripped" })),
    /-tripped-6ee4409dc3b1fe0f\.csv$/,
  );
  assert.match(downloadName(meta({ truncated: true })), /-truncated-/);

  // Two saves of the same run at different lengths must not collide: the
  // checksum covers every value emitted, so it moves with every extra sample.
  assert.notEqual(
    downloadName(meta({ rows: 40, checksum: "aaaa" })),
    downloadName(meta({ rows: 60, checksum: "bbbb" })),
  );

  // Nothing in the name needs escaping by a file system or a shell.
  assert.match(baseline, /^[A-Za-z0-9.+-]+$/);
});

test("saveCsv builds a blob, names the download, and revokes the URL", async () => {
  // An object URL that is never revoked pins its blob for the life of the
  // document, and one of these is a hundred megabytes.
  const created = [];
  const revoked = [];
  const clicks = [];
  const timers = [];
  const anchor = {
    click() {
      clicks.push({ href: anchor.href, download: anchor.download });
    },
  };
  const host = {
    document: { createElement: (tag) => (tag === "a" ? anchor : null) },
    URL: {
      createObjectURL(blob) {
        created.push(blob);
        return `blob:fake/${created.length}`;
      },
      revokeObjectURL: (href) => revoked.push(href),
    },
    setTimeout: (fn) => timers.push(fn),
  };

  const parts = ["# header\n", "1,2\n", "3,4\n"];
  const href = saveCsv(parts, "tep-baseline-1h-20rows-abcd.csv", host);

  assert.equal(created.length, 1);
  assert.equal(created[0].type, "text/csv;charset=utf-8");
  assert.equal(await created[0].text(), parts.join(""), "the blob is the parts, in order");
  assert.deepEqual(clicks, [{ href, download: "tep-baseline-1h-20rows-abcd.csv" }]);

  // Deferred by one task, because `click` only starts the download.
  assert.deepEqual(revoked, [], "revoked in the same task as the click");
  assert.equal(timers.length, 1);
  timers[0]();
  assert.deepEqual(revoked, [href]);
});

// The rest drive the real worker, so what is saved is what the page would save.

/** Run a scenario to the end, feeding every chunk into a recorder as `app.js` does. */
async function record(overrides = {}) {
  worker.drain();
  worker.send(await request(overrides));
  const started = await worker.next("started");
  const recorder = new Recorder(started.rowWidth);
  let last = null;
  for (;;) {
    const chunk = await worker.next("chunk");
    recorder.append(chunk.values, chunk.emittedBefore, chunk.activeFaults, chunk.faultAges);
    last = chunk;
    if (chunk.finished) break;
  }
  return { started, last, recorder, meta: exportMeta(started, last, recorder, version()) };
}

test("a saved run is the run: every planned sample, in order, at full resolution", async () => {
  const { started, recorder, meta: info } = await record({ hours: 2 });
  assert.equal(started.totalSamples, 40);
  assert.equal(recorder.count, 40);
  assert.equal(recorder.truncated, false);

  const file = parse(csvParts(recorder, info));
  assert.equal(file.rows.length, 40);
  assert.equal(file.columns.length, started.rowWidth + 3);

  // The times are the sampling cadence, and the last is the run length.
  const hours = file.rows.map((row) => Number(row[1]));
  assert.ok(Math.abs(hours.at(-1) - 2) < 1 / 3600 + 1e-9, `ended at ${hours.at(-1)} h`);
  for (let i = 1; i < hours.length; i += 1) {
    assert.ok(hours[i] > hours[i - 1], "the time column must be increasing");
  }
  // 3,600 steps an hour, sampled every 180: the last sample is step 7,200.
  assert.equal(Number(file.rows.at(-1)[0]), 7200);
  assert.equal(Number(file.rows[0][0]), 180);
});

test("the scenario in the header reads back as the run that produced the rows", async () => {
  // The claim the whole file rests on: it does not describe its run, it names
  // it, in the one format `Scenario.fromText` reads.
  const { started, recorder, meta: info } = await record({ hours: 1, faults: [6] });
  const file = parse(csvParts(recorder, info));

  const line = file.comments.find((c) => c.startsWith("# scenario: "));
  assert.ok(line, `no scenario line in:\n${file.comments.join("\n")}`);
  const token = line.slice("# scenario: ".length);

  const { Scenario } = await bindings();
  const back = Scenario.fromText(token);
  assert.equal(back.digest(), started.scenarioDigest, "the digest did not survive the file");
  assert.deepEqual([...back.activeFaults], [6]);
  back.free();

  // And the header says which run it was without anybody having to parse it.
  assert.ok(file.comments.some((c) => c.includes(started.scenarioDigest)));
  assert.match(downloadName(info), /^tep-idv6-/);
});

test("the ground truth column turns on at the sample the driver fires, not the chunk", async () => {
  // `temain_mod.f:366-368` switches IDV(12) on at eight hours whatever the
  // scenario asked for, and eight hours is a sampling instant. The chunk here
  // is three hours wide, so a label taken per chunk instead of per row would be
  // wrong by up to sixty samples and this bound would fail by a factor of 60.
  const { started, recorder, meta: info } = await record({
    hours: 10,
    driverForcesIdv12: true,
    chunkSamples: 60,
  });
  assert.deepEqual([...started.faults], [], "the scenario itself asks for no fault");

  const file = parse(csvParts(recorder, info));
  const first = file.rows.findIndex((row) => row[started.rowWidth + 1].includes("12"));
  assert.ok(first > 0, "IDV(12) never appeared in the label column");

  const sample = 180 / 3600;
  const onset = Number(file.rows[first][1]);
  assert.ok(
    Math.abs(onset - 8) < sample,
    `IDV(12) is labelled from ${onset} h, which is not within one sample of 8 h`,
  );
  assert.equal(file.rows[first - 1][started.rowWidth + 1], "", "the row before must be clean");
  assert.equal(Number(file.rows[first][started.rowWidth + 2]), 0, "the onset row's age is zero");

  // And the age counts up by the sampling interval from there.
  const next = Number(file.rows[first + 1][started.rowWidth + 2]);
  assert.ok(Math.abs(next - sample) < 1e-9, `the next row's age is ${next}`);
});

test("a run that tripped, and one that was stopped, both say so in the file", async () => {
  const tripped = await record({ hours: 6, controlled: false });
  assert.equal(tripped.last.outcome, "tripped");
  const trippedFile = parse(csvParts(tripped.recorder, tripped.meta));
  assert.ok(
    trippedFile.comments.some((c) => c.startsWith("# outcome: tripped at ")),
    `the header hid the trip:\n${trippedFile.comments.join("\n")}`,
  );
  assert.match(downloadName(tripped.meta), /-tripped-/);

  // A run stopped part way is the paused case: the button works and the file
  // says how far it got rather than claiming to be the plan.
  worker.drain();
  worker.send(await request({ hours: 48, chunkSamples: 5 }));
  const started = await worker.next("started");
  const recorder = new Recorder(started.rowWidth);
  let last = null;
  for (let i = 0; i < 3; i += 1) {
    last = await worker.next("chunk");
    recorder.append(last.values, last.emittedBefore, last.activeFaults, last.faultAges);
  }
  worker.send({ type: "stop" });

  const info = exportMeta(started, last, recorder, version());
  const file = parse(csvParts(recorder, info));
  assert.equal(recorder.count, 15);
  assert.equal(file.rows.length, 15);
  assert.ok(
    file.comments.some((c) => c.includes("15 of 960 planned samples")),
    `the header did not say how far it got:\n${file.comments.join("\n")}`,
  );
  assert.ok(file.comments.some((c) => c.includes("still running")));
});
