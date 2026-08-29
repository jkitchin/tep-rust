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
import {
  decimateMinMax,
  formatHours,
  formatReadout,
  formatValue,
  niceRange,
  shortUnit,
} from "../dist/js/format.js";
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

test("formatReadout is four significant figures and stays inside its box", () => {
  // The nominal steady state, which is what the diagram shows for most of a
  // fault-free run. Four figures each, and no unit is amputated to fit.
  assert.equal(formatReadout(2705.0), "2705");
  assert.equal(formatReadout(120.41), "120.4");
  assert.equal(formatReadout(75.0), "75.00");
  assert.equal(formatReadout(9.3477), "9.348");
  assert.equal(formatReadout(0.25052), "0.2505");
  // Below a tenth, four fixed decimals rather than four figures: the absolute
  // size is the readable fact about a flow that has collapsed.
  assert.equal(formatReadout(0.0012), "0.0012");
  assert.equal(formatReadout(0.00034), "3.4e-4");
  assert.equal(formatReadout(0), "0");
  assert.equal(formatReadout(Number.NaN), "-");
  assert.equal(formatReadout(Number.POSITIVE_INFINITY), "-");

  // The width budget the layout is placed against. Six characters, plus a sign
  // for the noisy near-zero readings a lost feed produces.
  const worst = [
    3200.5, 2705.0, 999.99, 120.41, 75.0, 26.902, 9.3477, 1.0, 0.25052, 0.0012,
    0.00034, -0.0035, -120.41, 1e-9, 5e5,
  ];
  for (const v of worst) {
    assert.ok(
      formatReadout(v).length <= 7,
      `formatReadout(${v}) is "${formatReadout(v)}", wider than the diagram budgets`,
    );
  }

  // The table's formatter is untouched: it is right for a column of readings
  // and this is a separate rule for a picture.
  assert.equal(formatValue(0.25052), "0.25052");
});

test("shortUnit abbreviates what the bindings spell out, and nothing else", () => {
  assert.equal(shortUnit("kPa gauge"), "kPa");
  assert.equal(shortUnit("Deg C"), "°C");
  assert.equal(shortUnit("Mole %"), "mol%");
  assert.equal(shortUnit("kg/hr"), "kg/h");
  assert.equal(shortUnit("m3/hr"), "m3/h");
  // Already short, or not ours: passed through verbatim rather than mapped to
  // a guess. This is what stops it from being a second table of units.
  assert.equal(shortUnit("kscmh"), "kscmh");
  assert.equal(shortUnit("%"), "%");
  assert.equal(shortUnit("kW"), "kW");
  assert.equal(shortUnit("furlongs per fortnight"), "furlongs per fortnight");
});

test("every unit the bindings produce fits the diagram once shortened", async () => {
  // Against the bindings themselves, not against a list here: a unit added or
  // respelled in `tepsim_core::variables` shows up as a failure rather than as
  // a readout running off the edge of the picture.
  const { columnUnits } = await bindings();
  const units = columnUnits();
  assert.ok(units.length > 50, `only ${units.length} units, the call is wrong`);
  for (const unit of units) {
    assert.ok(
      shortUnit(unit).length <= 5,
      `unit "${unit}" is "${shortUnit(unit)}" on the diagram, too wide for a readout`,
    );
  }
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
  "every=180;faults=;controlled=1;idv12=0;trip=1;continuous=0;integrator=euler;events=";

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

test("the flowsheet shows every continuous measurement", () => {
  // XMEAS(1..22) are the continuous plant readings: the plant as an operator
  // sees it, and a flowsheet that dropped one would be quietly incomplete.
  // XMEAS(23..41) are the analysers, 19 compositions that belong on trends
  // rather than on a diagram; three are shown where the composition is the
  // point of the stream (purge B, product G and H).
  const shown = new Set(READOUTS.map(columnOf));
  for (let i = 1; i <= 22; i += 1) {
    assert.ok(shown.has(i), `XMEAS(${i}) is missing from the flowsheet`);
  }
  for (const composition of [30, 40, 41]) {
    assert.ok(shown.has(composition), `XMEAS(${composition}) is missing`);
  }
});

// A valve earns a place on the diagram only where nothing else reports what it
// does. Nine of the twelve sit directly upstream of a measurement that is
// already drawn, so drawing both says the same thing twice and costs the space
// the readable labels needed; the agitator moves no material and holds at 50%
// in every published run. All twelve are on the trends and in the CSV.
//
// The table is the justification, written down. Adding a drop means naming
// what covers it, and the test then checks that the cover is itself on the
// diagram, so this cannot decay into a list of things quietly left out.
const VALVE_COVERED_BY = new Map([
  [1, 2], // D feed flow -> XMEAS(2) D feed
  [2, 3], // E feed flow -> XMEAS(3) E feed
  [3, 1], // A feed flow -> XMEAS(1) A feed
  [4, 4], // A and C feed flow -> XMEAS(4) A and C feed
  [6, 10], // purge valve -> XMEAS(10) purge rate
  [7, 14], // separator pot liquid -> XMEAS(14) underflow
  [8, 17], // stripper liquid product -> XMEAS(17) product rate
  [9, 19], // stripper steam valve -> XMEAS(19) steam flow
]);
const VALVE_MOVES_NOTHING = new Set([12]); // agitator speed

test("the flowsheet shows the valves nothing else reports, and covers the rest", () => {
  const shown = new Set(READOUTS.map(columnOf));

  for (let mv = 1; mv <= 12; mv += 1) {
    const onDiagram = shown.has(41 + mv);
    const covered = VALVE_COVERED_BY.get(mv);
    if (covered !== undefined) {
      assert.ok(!onDiagram, `XMV(${mv}) is drawn as well as XMEAS(${covered})`);
      assert.ok(
        shown.has(covered),
        `XMV(${mv}) was dropped in favour of XMEAS(${covered}), which is not drawn either`,
      );
    } else if (VALVE_MOVES_NOTHING.has(mv)) {
      assert.ok(!onDiagram, `XMV(${mv}) is drawn but was justified as a constant`);
    } else {
      assert.ok(onDiagram, `XMV(${mv}) is the only view of its flow and must be drawn`);
    }
  }

  // The three that survive the rule, named so that removing one is a visible
  // decision rather than an accident: the compressor recycle and the two
  // cooling water flows.
  assert.deepEqual(
    [5, 10, 11].filter((mv) => shown.has(41 + mv)),
    [5, 10, 11],
  );
});
