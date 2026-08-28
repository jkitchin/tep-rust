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
import { decodeState, encodeState } from "../dist/js/share.js";
import { READOUTS, columnOf } from "../dist/js/flowsheet.js";

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

// The link encoding. Every field a run's output depends on has to survive the
// round trip exactly, because the claim a link makes is that it reproduces a
// run and not that it approximates one.
const DEFAULTS = {
  seed: 4651207995,
  hours: 48,
  stepHours: 1 / 3600,
  sampleEvery: 180,
  integrator: "euler",
  controlled: true,
  driverForcesIdv12: true,
  tripEndsTheRun: false,
  faults: [],
  channels: [7, 9, 12, 15, 17, 40],
  chunkSamples: 20,
  speedMultiple: 0,
};

test("the baseline link is empty", () => {
  assert.equal(encodeState({ ...DEFAULTS }, DEFAULTS), "");
  assert.deepEqual(decodeState("", DEFAULTS, 54), DEFAULTS);
  assert.deepEqual(decodeState("#", DEFAULTS, 54), DEFAULTS);
});

test("every field survives the round trip", () => {
  const state = {
    seed: 1234567891,
    hours: 6.5,
    stepHours: 1 / 7200,
    sampleEvery: 60,
    integrator: "rk4",
    controlled: false,
    driverForcesIdv12: false,
    tripEndsTheRun: true,
    faults: [1, 6, 20],
    channels: [7, 41, 53],
    chunkSamples: 5,
    speedMultiple: 600,
  };
  const decoded = decodeState(`#${encodeState(state, DEFAULTS)}`, DEFAULTS, 54);
  assert.deepEqual(decoded, state);
});

test("the seed survives exactly", () => {
  // 4651207995 is what `teprob.f:1187` compiles in. A link that rounded it
  // would silently produce a different run from the one it claims to be.
  const state = { ...DEFAULTS, seed: 4651207995.25 };
  const decoded = decodeState(`#${encodeState(state, DEFAULTS)}`, DEFAULTS, 54);
  assert.equal(decoded.seed, 4651207995.25);
});

test("a hostile link opens the baseline rather than a broken page", () => {
  const nonsense =
    "#s=abc&h=-1&dt=0&n=0&i=leapfrog&c=maybe&f=0,21,99,7&p=0,54,9&k=-3&x=-1";
  const decoded = decodeState(nonsense, DEFAULTS, 54);
  assert.equal(decoded.seed, DEFAULTS.seed, "a non-numeric seed falls back");
  assert.equal(decoded.hours, DEFAULTS.hours, "a negative duration falls back");
  assert.equal(decoded.stepHours, DEFAULTS.stepHours);
  assert.equal(decoded.sampleEvery, DEFAULTS.sampleEvery);
  assert.equal(decoded.integrator, "euler", "an unknown integrator falls back");
  assert.equal(decoded.controlled, DEFAULTS.controlled);
  // Twenty disturbances, not twenty-one: `teprob.f:340` loops DO 500 I=1,20.
  assert.deepEqual(decoded.faults, [7], "only 1..20 survive");
  assert.deepEqual(decoded.channels, [9], "column 0 is time and 54 is past the end");
  assert.equal(decoded.chunkSamples, DEFAULTS.chunkSamples);
  assert.equal(decoded.speedMultiple, DEFAULTS.speedMultiple);
});

test("an empty fault list is distinguishable from an absent one", () => {
  // `#f=` must mean "no faults", not "use the default", or a link that clears
  // the panel would be indistinguishable from one that never touched it.
  const withFaults = { ...DEFAULTS, faults: [3] };
  const cleared = decodeState("#f=", withFaults, 54);
  assert.deepEqual(cleared.faults, []);
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
