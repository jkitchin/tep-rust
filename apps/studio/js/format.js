// Pure helpers shared by the charts and the flowsheet.
//
// They live apart from anything that touches the DOM so that
// `apps/studio/node/` can import and exercise them in Node, where there is no
// canvas to draw on. Axis padding and significant-figure rules are the kind of
// thing that looks obviously right and is off by one at the boundary.

/**
 * Pad a data range into an axis range.
 *
 * A flat trace is the case worth thinking about. Reactor level sits at 75.0 for
 * the first hour of a baseline run, so `lo === hi` and a naive axis has zero
 * height and divides by zero. Give it a band proportional to the value, or a
 * unit band when the value is zero.
 *
 * @param {number} lo
 * @param {number} hi
 * @param {number} fraction how much of the span to add at each end
 * @returns {{lo: number, hi: number}}
 */
export function niceRange(lo, hi, fraction = 0.08) {
  if (!Number.isFinite(lo) || !Number.isFinite(hi)) return { lo: 0, hi: 1 };
  const span = hi - lo;
  if (span <= 0) {
    const pad = Math.abs(hi) * 0.05 || 0.5;
    return { lo: hi - pad, hi: hi + pad };
  }
  const pad = span * fraction;
  return { lo: lo - pad, hi: hi + pad };
}

/**
 * A number for a label: enough digits to distinguish neighbouring values, few
 * enough to read.
 *
 * Fixed decimals rather than significant figures, because a column of readings
 * that share a decimal point is scannable and a column of `2.7052e+3` is not.
 * The magnitude picks the count.
 *
 * @param {number} value
 * @returns {string}
 */
export function formatValue(value) {
  if (!Number.isFinite(value)) return "-";
  const magnitude = Math.abs(value);
  if (magnitude === 0) return "0";
  if (magnitude >= 1000) return value.toFixed(1);
  if (magnitude >= 100) return value.toFixed(2);
  if (magnitude >= 1) return value.toFixed(3);
  if (magnitude >= 0.001) return value.toFixed(5);
  return value.toExponential(2);
}

/** Simulated hours as `h:mm`, for a time axis. */
export function formatHours(hours) {
  if (!Number.isFinite(hours)) return "-";
  const whole = Math.floor(hours);
  const minutes = Math.round((hours - whole) * 60);
  // 59.7 minutes rounds to 60, which should read as the next hour, not as
  // ":60". Carry it.
  if (minutes === 60) return `${whole + 1}:00`;
  return `${whole}:${String(minutes).padStart(2, "0")}`;
}

/**
 * Reduce a column to at most two points per pixel, preserving extremes.
 *
 * A 20,000-point trace on a 900-pixel canvas is 20,000 line segments to draw
 * 900 pixels, twenty times a second, on every one of a dozen charts. Drawing
 * the minimum and maximum in each pixel column instead is visually identical
 * (a spike still shows, because a spike is an extreme) and is bounded by the
 * canvas width rather than by the run length.
 *
 * Returns a flat array of `[x, yMin, yMax]` triples in pixel-column order.
 *
 * @param {import("./history.js").History} history
 * @param {number} col column index into a packed row
 * @param {number} pixels width of the plot area
 * @returns {Float64Array}
 */
export function decimateMinMax(history, col, pixels) {
  // An empty history has no extremes to preserve, and asking for one column of
  // nothing would read past the end of the store.
  if (history.count === 0) return new Float64Array(0);
  const columns = Math.max(1, Math.min(pixels | 0, history.count));
  const out = new Float64Array(columns * 3);
  const perColumn = history.count / columns;
  for (let c = 0; c < columns; c += 1) {
    const from = Math.floor(c * perColumn);
    const to = Math.max(from + 1, Math.floor((c + 1) * perColumn));
    let lo = Number.POSITIVE_INFINITY;
    let hi = Number.NEGATIVE_INFINITY;
    for (let i = from; i < to && i < history.count; i += 1) {
      const v = history.at(i, col);
      if (v < lo) lo = v;
      if (v > hi) hi = v;
    }
    // The x of the column is the time of its first row, so the trace lines up
    // with the time axis rather than with the sample index.
    out[c * 3] = history.timeAt(Math.min(from, history.count - 1));
    out[c * 3 + 1] = lo;
    out[c * 3 + 2] = hi;
  }
  return out;
}
