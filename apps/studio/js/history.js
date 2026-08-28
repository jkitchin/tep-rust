// A bounded store for the rows the worker sends back.
//
// The charts and the flowsheet both need history: a chart has to rescale when a
// trace leaves its axis, and rescaling means looking at every point again, not
// just the new ones. So the page keeps the run.
//
// The run can be large. `sampleEvery` is 180 by default, which puts a 48-hour
// run at 960 samples and 414 kB, but the bindings will happily accept
// `sampleEvery = 1`, and then the same run is 172,800 samples and 74 MB. The
// binding-level ceiling is a million samples, which is 432 MB. A browser tab
// should not decide to allocate that because someone typed a 1 in a box.
//
// So the store has a fixed capacity and, when it fills, throws away every other
// row and doubles its stride. The time span it covers is always the whole run;
// what degrades is resolution, gracefully and by a factor of two at a time.
// That is the right trade for a display: nobody can see 172,800 points on a
// 900-pixel canvas anyway. The full-fidelity path off this page is the CSV
// export, which streams straight from the chunks and never comes here.
//
// Halving is amortized O(1) per row: each halve costs `capacity / 2` copies and
// buys `capacity / 2` more appends.

export class History {
  /**
   * @param {number} rowWidth values per packed row, `Sim.rowWidth`
   * @param {number} capacity rows to keep before decimating
   */
  constructor(rowWidth, capacity = 20_000) {
    this.rowWidth = rowWidth;
    this.capacity = Math.max(2, capacity | 0);
    this.values = new Float64Array(this.capacity * rowWidth);
    // Rows kept.
    this.count = 0;
    // Rows offered, kept or not. A row is kept when `seen % decimation === 0`,
    // which is what makes halving consistent: the survivors of a halve are
    // exactly the rows whose index divides the doubled stride.
    this.seen = 0;
    this.decimation = 1;
  }

  /** Forget everything. Capacity and row width are unchanged. */
  clear() {
    this.count = 0;
    this.seen = 0;
    this.decimation = 1;
  }

  /**
   * Append a chunk: `values.length / rowWidth` rows, row-major.
   *
   * @param {Float64Array} values
   */
  append(values) {
    const width = this.rowWidth;
    const rows = Math.floor(values.length / width);
    for (let r = 0; r < rows; r += 1) {
      if (this.seen % this.decimation === 0) {
        if (this.count === this.capacity) this.halve();
        // The halve may have moved this row out of the kept set: after
        // doubling, a row that was `seen % d === 0` is only kept if it is also
        // `seen % 2d === 0`.
        if (this.seen % this.decimation === 0) {
          this.values.set(
            values.subarray(r * width, r * width + width),
            this.count * width,
          );
          this.count += 1;
        }
      }
      this.seen += 1;
    }
  }

  /** Drop every other kept row and double the stride. */
  halve() {
    const width = this.rowWidth;
    const kept = Math.ceil(this.count / 2);
    for (let i = 1; i < kept; i += 1) {
      this.values.copyWithin(i * width, 2 * i * width, (2 * i + 1) * width);
    }
    this.count = kept;
    this.decimation *= 2;
  }

  /** Value of column `col` in kept row `row`. */
  at(row, col) {
    return this.values[row * this.rowWidth + col];
  }

  /** Simulated hours of kept row `row`. Column 0 is the time. */
  timeAt(row) {
    return this.values[row * this.rowWidth];
  }

  /** The most recent kept row as a subarray, or `null` when empty. */
  lastRow() {
    if (this.count === 0) return null;
    const width = this.rowWidth;
    return this.values.subarray((this.count - 1) * width, this.count * width);
  }

  /**
   * Extremes of one column over the kept rows.
   *
   * Returns `null` if there is nothing finite to report, which happens before
   * the first chunk and would otherwise hand the charts an infinite axis.
   *
   * @param {number} col
   * @returns {{lo: number, hi: number} | null}
   */
  extent(col) {
    let lo = Number.POSITIVE_INFINITY;
    let hi = Number.NEGATIVE_INFINITY;
    const width = this.rowWidth;
    for (let i = 0; i < this.count; i += 1) {
      const v = this.values[i * width + col];
      if (!Number.isFinite(v)) continue;
      if (v < lo) lo = v;
      if (v > hi) hi = v;
    }
    return lo <= hi ? { lo, hi } : null;
  }
}
