// The pure parts of the page: the history store, the formatting, the link
// encoding, and the flowsheet's binding table.
//
// None of these touch the DOM, so Node can run them exactly as the browser
// does. They are the parts most likely to be wrong in a way that looks right:
// a decimating store that quietly drops the newest row, an axis that divides by
// zero on a flat trace, a link that round-trips every field except the one that
// mattered, a readout wired to the wrong Fortran index.

import test from "node:test";
import assert from "node:assert/strict";

import { History } from "../dist/js/history.js";
import { decimateMinMax, formatHours, formatValue, niceRange } from "../dist/js/format.js";
import { decodeLink, encodeLink } from "../dist/js/share.js";
import { READOUTS, columnOf } from "../dist/js/flowsheet.js";
import { bindings } from "./harness.mjs";

const ROW_WIDTH = 54;

/** `rows` packed rows whose time column counts up and whose column 1 is `n`. */
function rows(from, count, width = ROW_WIDTH) {
  const out = new Float64Array(count * width);
  for (let i = 0; i < count; i += 1) {
    out[i * width] = (from + i) / 20;
    out[i * width + 1] = from + i;
  }
  return out;
}

test("History keeps every row until it is full", () => {
  const h = new History(ROW_WIDTH, 100);
  h.append(rows(0, 40));
  h.append(rows(40, 40));
  assert.equal(h.count, 80);
  assert.equal(h.decimation, 1);
  assert.equal(h.at(0, 1), 0);
  assert.equal(h.at(79, 1), 79);
  assert.equal(h.lastRow()[1], 79);
});

test("History decimates instead of growing, and still spans the whole run", () => {
  const capacity = 64;
  const h = new History(ROW_WIDTH, capacity);
  // Ten times the capacity. A run that outgrows the store must lose
  // resolution, never its ends.
  h.append(rows(0, capacity * 10));

  assert.ok(h.count <= capacity, `kept ${h.count} rows in a store of ${capacity}`);
  assert.ok(h.decimation > 1, "the store should have decimated");
  assert.equal(h.at(0, 1), 0, "the first sample survives");

  // Every kept row is a real row, in order, at the stated stride.
  for (let i = 1; i < h.count; i += 1) {
    assert.equal(h.at(i, 1) - h.at(i - 1, 1), h.decimation);
  }
  // The last kept row is within one stride of the newest sample offered.
  assert.ok(capacity * 10 - 1 - h.at(h.count - 1, 1) < h.decimation);
});

test("History extents ignore nothing and report nothing when empty", () => {
  const h = new History(ROW_WIDTH, 100);
  assert.equal(h.extent(1), null);
  h.append(rows(0, 10));
  assert.deepEqual(h.extent(1), { lo: 0, hi: 9 });
});

test("niceRange gives a flat trace an axis to sit on", () => {
  // Reactor level sits at exactly 75.0 for the first hour of a baseline run.
  // A zero-height axis divides by zero and paints nothing.
  const flat = niceRange(75, 75);
  assert.ok(flat.hi > flat.lo, "a flat trace must still get a band");
  assert.ok(flat.lo < 75 && flat.hi > 75);

  const zero = niceRange(0, 0);
  assert.ok(zero.hi > zero.lo);

  const normal = niceRange(0, 10, 0.1);
  assert.equal(normal.lo, -1);
  assert.equal(normal.hi, 11);
});

test("formatValue keeps a column of readings scannable", () => {
  assert.equal(formatValue(2705.2), "2705.2");
  assert.equal(formatValue(120.41), "120.41");
  assert.equal(formatValue(0.25052), "0.25052");
  assert.equal(formatValue(0.0012), "0.00120");
  // Below a thousandth, fixed decimals stop distinguishing anything and the
  // exponent is the only readable form.
  assert.equal(formatValue(0.00034), "3.40e-4");
  assert.equal(formatValue(0), "0");
  assert.equal(formatValue(Number.NaN), "-");
  assert.equal(formatValue(Number.POSITIVE_INFINITY), "-");
});

test("formatHours carries 59.7 minutes into the next hour", () => {
  assert.equal(formatHours(0), "0:00");
  assert.equal(formatHours(1.5), "1:30");
  assert.equal(formatHours(0.995), "1:00", "59.7 minutes must not read as 0:60");
  assert.equal(formatHours(48), "48:00");
});

test("decimateMinMax is bounded by the canvas and keeps the extremes", () => {
  const h = new History(ROW_WIDTH, 10_000);
  const count = 5000;
  const values = new Float64Array(count * ROW_WIDTH);
  for (let i = 0; i < count; i += 1) {
    values[i * ROW_WIDTH] = i / 100;
    values[i * ROW_WIDTH + 1] = i;
  }
  // One spike, in the middle, which must survive being drawn on 300 pixels.
  values[2500 * ROW_WIDTH + 1] = 99_999;
  h.append(values);

  const points = decimateMinMax(h, 1, 300);
  assert.equal(points.length, 300 * 3, "one triple per pixel column");

  let highest = Number.NEGATIVE_INFINITY;
  for (let c = 0; c < 300; c += 1) highest = Math.max(highest, points[c * 3 + 2]);
  assert.equal(highest, 99_999, "a spike must survive decimation");

  // Fewer samples than pixels: one column per sample, nothing invented.
  const small = new History(ROW_WIDTH, 100);
  small.append(rows(0, 7));
  assert.equal(decimateMinMax(small, 1, 300).length, 7 * 3);

  // Nothing at all: no columns, and no read past the end of the store.
  assert.equal(decimateMinMax(new History(ROW_WIDTH, 100), 1, 300).length, 0);
});

// The link encoding.
//
// Since B-0054a the scenario is one opaque token here: `Scenario.text` produces
// it, `Scenario.fromText` consumes it, and this file neither reads nor writes a
// field of it. That is the point. The three hand-maintained field lists that
// used to keep a link, a worker request and a settings panel in step could not
// be kept in step, and a field added to `Scenario` reached none of them.
//
// What is tested here is what is left: that the token survives the fragment
// byte for byte, and that the view fields, which reach no arithmetic, fall back
// rather than throw.
const BASELINE =
  "tepsim.scenario.v1;seed=4651207995;hours=48;step=2.777777777777778e-4;" +
  "every=180;faults=;controlled=1;idv12=1;trip=0;continuous=0;integrator=euler;events=";

const DEFAULTS = {
  scenario: BASELINE,
  channels: [7, 9, 12, 15, 17, 40],
  chunkSamples: 20,
  speedMultiple: 0,
};

test("the baseline link is empty", () => {
  assert.equal(encodeLink({ ...DEFAULTS }, DEFAULTS), "");
  assert.deepEqual(decodeLink("", DEFAULTS, 54), DEFAULTS);
  assert.deepEqual(decodeLink("#", DEFAULTS, 54), DEFAULTS);
});

test("the scenario token survives the fragment byte for byte", () => {
  // Semicolons, equals signs, commas and colons all appear inside the token and
  // all are fragment-safe. `decodeLink` splits a pair at its *first* `=`, which
  // is the whole reason a scenario full of them can be a value at all.
  const scenario =
    "tepsim.scenario.v1;seed=4651207995.25;hours=6.5;step=1.388888888888889e-4;" +
    "every=60;faults=1,6,20;controlled=0;idv12=0;trip=1;continuous=1;" +
    "integrator=rk4;events=8:start:6,12:magnitude:13:0.5,20:setpoint:9:0.25";
  const link = { ...DEFAULTS, scenario };
  const fragment = encodeLink(link, DEFAULTS);
  assert.ok(fragment.includes(scenario), "the token was mangled on the way out");
  assert.equal(decodeLink(`#${fragment}`, DEFAULTS, 54).scenario, scenario);
});

test("every view field survives the round trip", () => {
  const link = {
    scenario: BASELINE,
    channels: [7, 41, 53],
    chunkSamples: 5,
    speedMultiple: 600,
  };
  assert.deepEqual(decodeLink(`#${encodeLink(link, DEFAULTS)}`, DEFAULTS, 54), link);
});

test("a hostile link opens the baseline rather than a broken page", () => {
  // The view fields fall back. The scenario token is handed back as written and
  // rejected by the bindings, which is where a rejection can carry a reason.
  const nonsense = "#q=leapfrog&p=0,54,9&k=-3&x=-1";
  const decoded = decodeLink(nonsense, DEFAULTS, 54);
  assert.equal(decoded.scenario, "leapfrog", "the token is not interpreted here");
  assert.deepEqual(decoded.channels, [9], "column 0 is time and 54 is past the end");
  assert.equal(decoded.chunkSamples, DEFAULTS.chunkSamples);
  assert.equal(decoded.speedMultiple, DEFAULTS.speedMultiple);

  // A malformed percent escape must not throw: `decodeURIComponent` does.
  assert.doesNotThrow(() => decodeLink("#q=%", DEFAULTS, 54));
});

test("an empty channel list is distinguishable from an absent one", () => {
  // `#p=` must mean "plot nothing", not "use the default", or a link that
  // cleared the picker would be indistinguishable from one that never touched
  // it.
  assert.deepEqual(decodeLink("#p=", DEFAULTS, 54).channels, []);
  assert.deepEqual(decodeLink("#", DEFAULTS, 54).channels, DEFAULTS.channels);
});

test("a link and the bindings agree on the token", async () => {
  // The claim the whole design rests on: what `encodeLink` puts in a fragment
  // is exactly what `Scenario.fromText` reads back, with no encoding step in
  // between that either side could get wrong.
  const { Scenario } = await bindings();
  const built = new Scenario();
  built.hours = 6.5;
  built.setFault(6, true);
  built.setIntegrator("rk4");
  const token = built.text;
  const digest = built.digest();
  built.free();

  const fragment = encodeLink({ ...DEFAULTS, scenario: token }, DEFAULTS);
  const back = Scenario.fromText(decodeLink(`#${fragment}`, DEFAULTS, 54).scenario);
  assert.equal(back.digest(), digest, "the digest did not survive the link");
  assert.equal(back.text, token);
  back.free();
});

test("the baseline token in this file is the one the bindings produce", async () => {
  // A golden, so a change to the format shows up here rather than as a link
  // that silently stops working.
  const { Scenario } = await bindings();
  const probe = new Scenario();
  assert.equal(probe.text, BASELINE);
  probe.free();
});

test("every flowsheet readout points at a real column", () => {
  const seen = new Set();
  for (const readout of READOUTS) {
    const named = (readout.meas === undefined) !== (readout.mv === undefined);
    assert.ok(named, `readout ${JSON.stringify(readout)} must name exactly one index`);
    if (readout.meas !== undefined) {
      assert.ok(readout.meas >= 1 && readout.meas <= 41, `XMEAS(${readout.meas})`);
    } else {
      assert.ok(readout.mv >= 1 && readout.mv <= 12, `XMV(${readout.mv})`);
    }
    const col = columnOf(readout);
    assert.ok(col >= 1 && col < ROW_WIDTH, `column ${col} is off the end of a row`);
    assert.ok(!seen.has(col), `column ${col} appears on the diagram twice`);
    seen.add(col);
    assert.ok(readout.caption.length > 0);
  }
});

test("the flowsheet shows every manipulated variable and every continuous measurement", () => {
  // XMEAS(1..22) are the continuous plant readings and XMV(1..12) the valves:
  // together they are the plant as an operator sees it, and a flowsheet that
  // dropped one would be quietly incomplete. XMEAS(23..41) are the analysers,
  // 19 compositions that belong on trends rather than on a diagram; two are
  // shown where they are the point of the stream (purge B, product G and H).
  const shown = new Set(READOUTS.map(columnOf));
  for (let i = 1; i <= 22; i += 1) {
    assert.ok(shown.has(i), `XMEAS(${i}) is missing from the flowsheet`);
  }
  for (let i = 1; i <= 12; i += 1) {
    assert.ok(shown.has(41 + i), `XMV(${i}) is missing from the flowsheet`);
  }
});
