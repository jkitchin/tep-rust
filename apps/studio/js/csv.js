// The run, as a file.
//
// A link reproduces a run; this saves what the run produced. The two are not
// alternatives. A link is exact and weighs nothing, but it needs this page and
// a browser to become numbers again, and somebody who wants to fit a detector
// to a fault wants the numbers now, in a file, in the tool they already use.
//
// # The full record, not the decimated one
//
// `js/history.js` keeps 20,000 rows and halves its stride when it fills,
// because a chart cannot show more points than it has pixels and a browser tab
// should not allocate 432 MB because someone typed a 1 in a box. That trade is
// right for a display and wrong for a file: a decimated CSV silently answers a
// different question from the one it was asked, and nothing downstream can tell.
//
// So the recorder here is a second, independent store, and it keeps every
// sample. It costs nothing extra to do so: the chunk arrives on this thread as
// a transferred `Float64Array` that the page already owns, `History.append`
// copies what it wants out of it, and holding on to the same array afterwards
// is a reference rather than a copy.
//
// What it does have is a ceiling, and the ceiling *stops recording* rather than
// thinning what it has. 200,000 rows is 86 MB retained and covers the longest
// run the panel can ask for at full cadence (48 hours at `sampleEvery = 1` is
// 172,800 samples). Past it the file is a prefix of the run, at full
// resolution, and both the header of the file and the readout on the page say
// so. A prefix is honest and a thinned record labelled as the run is not.
//
// # The scenario travels in the file
//
// `tepsim.scenario.v1;...` is the whole run: `Scenario.fromText` turns it back
// into bit-identical numbers on any platform. It goes in a comment line rather
// than in the file name, because a name is the wrong place for it. The token
// contains `;` and `:`, which are illegal in a Windows file name and which
// every browser's download path rewrites, so a token in the name would arrive
// mangled and stop round-tripping without saying it had. A comment also
// survives the file being renamed, which is the first thing that happens to
// anything in a downloads folder.
//
// Comment lines are `#`, which `pandas.read_csv(comment="#")`,
// `numpy.loadtxt`, and `read.csv(comment.char = "#")` all skip by default or
// with one argument.
//
// # Numbers
//
// `String(value)` is JavaScript's shortest representation that reads back as
// the same `f64`, so the file is exact and about a third the size of the
// seventeen-digit form. `tepsim::recorder::Csv` writes `{:.17e}` instead, which
// is the same guarantee reached from the other side; both round-trip, and
// neither is an approximation of the other.

/** Rows to record before recording stops. See the note above. */
export const MAX_RECORDED_ROWS = 200_000;

/**
 * Every sample of a run, held as the chunks it arrived in.
 *
 * The mirror of `History`: same input, opposite trade. `History` bounds its
 * memory by losing resolution; this bounds its memory by stopping, and reports
 * that it stopped.
 */
export class Recorder {
  /**
   * @param {number} rowWidth values per packed row, `Sim.rowWidth`
   * @param {number} maxRows rows to record before recording stops
   */
  constructor(rowWidth, maxRows = MAX_RECORDED_ROWS) {
    this.rowWidth = rowWidth;
    this.maxRows = Math.max(1, maxRows | 0);
    /**
     * One entry per chunk kept: the packed values, how many of its rows are
     * recorded, the emission index of its first row, the simulated hours of
     * its last recorded row, and its ground truth.
     * @type {{values: Float64Array, rows: number, first: number, last: number,
     *         faults: number[], ages: number[]}[]}
     */
    this.chunks = [];
    /** Rows recorded. */
    this.count = 0;
    /** Rows the ceiling refused. Non-zero means the file is a prefix. */
    this.dropped = 0;
    /** Rows in the largest chunk kept, which bounds the ground truth's age. */
    this.largestChunk = 0;
    /** Simulated hours of the last recorded row. */
    this.lastHours = 0;
  }

  /** Whether the ceiling has cut the record short. */
  get truncated() {
    return this.dropped > 0;
  }

  /**
   * Record a chunk, ground truth and all.
   *
   * The array is kept, not copied. The page owns it once the worker has
   * transferred it and nothing else writes to it.
   *
   * @param {Float64Array} values `values.length / rowWidth` packed rows
   * @param {number} first emission index of the chunk's first row
   * @param {number[]} faults one-based IDV flags on at the chunk's last sample
   * @param {number[]} ages how long each had been on by then, in hours
   */
  append(values, first, faults, ages) {
    const offered = Math.floor(values.length / this.rowWidth);
    const rows = Math.min(offered, this.maxRows - this.count);
    this.dropped += offered - rows;
    if (rows <= 0) return;
    // The ages the worker measured are relative to the chunk's last sample. If
    // the ceiling cut the chunk short that row is no longer in the file, but it
    // is still the row the ages are relative to, so the age of the row the
    // labels were taken at is what matters and not which rows were kept.
    const last = values[(offered - 1) * this.rowWidth];
    this.chunks.push({ values, rows, first, last, faults, ages });
    this.count += rows;
    this.largestChunk = Math.max(this.largestChunk, rows);
    this.lastHours = values[(rows - 1) * this.rowWidth];
  }
}

/**
 * What a saved file says it is, from the run that produced the rows.
 *
 * Everything comes from the worker's own messages rather than from the panel,
 * because the panel can be edited while a run is in flight and then it is
 * describing a different run. Here rather than in `app.js` so that the tests in
 * `apps/studio/node/` label a file exactly as the page does, instead of
 * carrying a second copy of this list that nothing keeps in step.
 *
 * @param {object} started the run's "started" message
 * @param {object | null} chunk its most recent "chunk" message
 * @param {Recorder} recorder
 * @param {string} version `tepsim_wasm.version()`
 */
export function exportMeta(started, chunk, recorder, version) {
  return {
    version,
    columnIds: started.columnIds,
    sampleEvery: started.sampleEvery,
    scenario: started.scenarioText,
    scenarioDigest: started.scenarioDigest,
    faults: started.faults,
    planned: started.totalSamples,
    checksum: chunk?.checksum ?? "0",
    outcome: chunk?.outcome ?? null,
    finished: chunk?.finished ?? false,
    tripHours: chunk?.tripHours ?? null,
    tripCause: chunk?.tripCause ?? null,
    // The last row actually recorded, which is behind the simulation's clock
    // once the ceiling has cut in, and behind the plan whenever a run was
    // stopped part way.
    hours: recorder.lastHours,
    rows: recorder.count,
    truncated: recorder.truncated,
  };
}

/**
 * Which of a chunk's disturbances were on at one of its rows, and for how long.
 *
 * The worker reports ground truth once per chunk, as of the chunk's last
 * sample, and a chunk can be hours of plant. Attributing that to every row in
 * the chunk would move an onset earlier by up to a chunk, which is exactly the
 * number a detection-delay figure is measuring. So the age of each active
 * disturbance travels with the chunk and every row is decided against its own
 * clock.
 *
 * Everything below is phrased as `behind`, the row's distance back from the
 * sample the ages were measured at, and never as an absolute onset time. A
 * disturbance's onset is recorded at the sample it first acted on, so that
 * sample is faulted and its age there is zero: for it `behind` and `age` are
 * the same subtraction of the same two numbers, so the boundary is exact
 * rather than exact to within an ulp. The driver forces IDV(12) at eight
 * hours, which is a sampling instant, so that boundary is reached on an
 * ordinary run and not only in principle.
 *
 * The one case this cannot resolve is a disturbance switched *off* part way
 * through a chunk by a scheduled `stop` event: it has already left the chunk's
 * list and the rows before the stop read as clean. Nothing but the scenario
 * line records that, and the scenario line does.
 *
 * @returns {{faults: number[], sinceOnset: number | null}}
 */
function labelAt(chunk, hours) {
  const behind = chunk.last - hours;
  const faults = [];
  let sinceOnset = null;
  for (let i = 0; i < chunk.faults.length; i += 1) {
    if (behind > chunk.ages[i]) continue;
    faults.push(chunk.faults[i]);
    // `tepsim::recorder::Csv` reports the age of the first listed fault, and
    // this is the same column.
    if (sinceOnset === null) sinceOnset = chunk.ages[i] - behind;
  }
  return { faults, sinceOnset };
}

/**
 * The header block, comments and column names, ending in a newline.
 *
 * The scenario token is what the bindings serialised, so it is one line of
 * `A-Z a-z 0-9 - . _ ~ ; = , : +` and cannot break out of a comment. Nothing
 * else in here came from anywhere a user can type.
 */
function header(recorder, meta) {
  const lines = [
    `TEP Studio export, tepsim-wasm ${meta.version}.`,
    `scenario: ${meta.scenario}`,
    `scenario digest ${meta.scenarioDigest}, output checksum ${meta.checksum}`,
    ...describeRows(recorder, meta),
    "fault and hours_since_onset are ground truth. A disturbance is listed from",
    "the moment it came on; a scheduled stop is resolved only to the chunk it",
    `landed in, at most ${recorder.largestChunk} samples. The scenario line above is exact.`,
  ];
  const columns = ["step", ...meta.columnIds, "fault", "hours_since_onset"];
  return `${lines.map((line) => `# ${line}`).join("\n")}\n${columns.join(",")}\n`;
}

/** What was recorded and how the run ended, as comment lines. */
function describeRows(recorder, meta) {
  const rows = recorder.count.toLocaleString("en-US");
  const planned = meta.planned.toLocaleString("en-US");
  const resolution =
    "at full resolution: every sample the run emitted, not the decimated " +
    "view the trend charts draw.";
  const what = recorder.truncated
    ? `${rows} samples ${resolution} Recording stopped at the ` +
      `${recorder.maxRows.toLocaleString("en-US")} row ceiling and the run went on, ` +
      `so this is the beginning of it, not all of it.`
    : `${rows} of ${planned} planned samples ${resolution}`;
  return [what, `outcome: ${describeOutcome(meta)}`];
}

function describeOutcome(meta) {
  const at = `${trim(meta.hours, 4)} h`;
  if (meta.outcome === "tripped") {
    const cause = meta.tripCause ?? "shutdown";
    return `tripped at ${trim(meta.tripHours, 4)} h (${cause}), recorded to ${at}`;
  }
  if (meta.finished) return `${meta.outcome ?? "finished"}, recorded to ${at}`;
  return `still running at ${at}; the recording was taken mid-run`;
}

/** A number for prose or a file name: fixed decimals, trailing zeros gone. */
function trim(value, decimals) {
  if (!Number.isFinite(value)) return "0";
  return String(Number(value.toFixed(decimals)));
}

/**
 * The recorded run as an array of strings for a `Blob`.
 *
 * An array rather than one string on purpose. A 48-hour run at full cadence is
 * 172,800 rows and about 100 MB of text, and `Blob` concatenates its parts
 * itself without ever materialising that as a JavaScript string.
 *
 * @param {Recorder} recorder
 * @param {object} meta see `header`
 * @returns {string[]}
 */
export function csvParts(recorder, meta) {
  const parts = [header(recorder, meta)];
  for (const chunk of recorder.chunks) parts.push(chunkText(chunk, recorder.rowWidth, meta));
  return parts;
}

/** One chunk's rows, newline terminated. */
function chunkText(chunk, width, meta) {
  const lines = [];
  const cells = new Array(width + 3);
  for (let r = 0; r < chunk.rows; r += 1) {
    const base = r * width;
    const hours = chunk.values[base];
    // The packed row is `[hours, XMEAS(1..41), XMV(1..12)]` and carries no step
    // number, so this is the one reconstructed column. `Simulation::step`
    // records a sample on every `sampleEvery`th step and counts from one, so
    // the nth sample emitted is step `n * sampleEvery`, which is exactly
    // `Sample::step`.
    cells[0] = String((chunk.first + r + 1) * meta.sampleEvery);
    for (let c = 0; c < width; c += 1) cells[c + 1] = String(chunk.values[base + c]);
    const { faults, sinceOnset } = labelAt(chunk, hours);
    cells[width + 1] = faults.join(" ");
    cells[width + 2] = sinceOnset === null ? "" : String(sinceOnset);
    lines.push(cells.join(","));
  }
  return `${lines.join("\n")}\n`;
}

/**
 * A file name that says what is in the file and does not collide with the last
 * one.
 *
 * The fault is first, because that is what somebody scanning a folder of these
 * is looking for. The checksum is last and is what keeps two downloads apart:
 * it is a digest of every value emitted so far, so it moves with every extra
 * sample, and two downloads that do collide are two downloads of the same
 * numbers. A clock would be the usual answer and is not available here, because
 * nothing in this app may read one: see the head of `js/app.js`.
 *
 * @param {object} meta
 * @returns {string}
 */
export function downloadName(meta) {
  const fault = meta.faults.length > 0 ? `idv${[...meta.faults].join("+")}` : "baseline";
  const notes = [
    `${trim(meta.hours, 3)}h`,
    `${meta.rows}rows`,
    ...(meta.outcome === "tripped" ? ["tripped"] : []),
    ...(meta.truncated ? ["truncated"] : []),
    meta.checksum,
  ];
  return `tep-${fault}-${notes.join("-")}.csv`;
}

/**
 * Hand `parts` to the browser as a download called `name`.
 *
 * `host` is a seam for `apps/studio/node/`, which has no DOM: an object URL
 * that is never revoked pins its blob for the life of the document, and for a
 * hundred-megabyte export that is the whole tab, so it is worth a test rather
 * than a comment claiming it happens.
 *
 * @param {string[]} parts
 * @param {string} name
 * @param {object} host something with `document`, `URL` and `setTimeout`
 * @returns {string} the object URL, for a test to check against
 */
export function saveCsv(parts, name, host = globalThis) {
  const blob = new Blob(parts, { type: "text/csv;charset=utf-8" });
  const href = host.URL.createObjectURL(blob);
  const anchor = host.document.createElement("a");
  anchor.href = href;
  anchor.download = name;
  // Not added to the document. A detached anchor's activation behaviour is the
  // same one a click on a link runs, and appending it would put a stray
  // element in the layout for as long as this takes.
  anchor.click();
  // Revoked on the next task, not here. `click` starts the download and does
  // not wait for it, and a URL revoked in the same task has been seen to abort
  // it in WebKit.
  host.setTimeout(() => host.URL.revokeObjectURL(href), 0);
  return href;
}
