// The trend grid: one small canvas chart per selected channel.
//
// # Why not uPlot
//
// PLAN.org names uPlot, and for the app it describes (dozens of series, cursor
// sync, pan and zoom) that is the right call: it is canvas based, tens of
// kilobytes, and streams at display rate in a way a hand-rolled chart will not
// match without becoming uPlot.
//
// This grid is not that app. Every chart here shares one time axis fixed at the
// run length, has exactly one series, and never pans or zooms; what it has to
// do well is repaint a dozen small plots as chunks arrive without dropping
// frames. That is about a hundred and fifty lines against a third-party
// dependency, a vendored copy to keep current, and its licence to carry in the
// bundle. It also keeps the byte budget honest: the whole page is currently the
// bindings plus its own source and nothing else.
//
// The moment this page wants a shared cursor across charts, zoom, or forty
// series on one axis, vendor uPlot and delete this file. The `TrendGrid`
// surface is `setColumns`, `setSelection` and `draw`, which is small enough to
// reimplement over uPlot in an afternoon. That is the point of it being small.
//
// # Repaint policy
//
// Chunks arrive faster than a display refreshes. `draw` is therefore called
// from a `requestAnimationFrame` in `app.js` rather than from the message
// handler, so a hundred chunks in a frame cost one repaint, not a hundred.

import { decimateMinMax, formatValue, formatHours, niceRange } from "./format.js";

// Deliberately mid-saturation and mid-lightness so the same six read on both a
// white and a near-black background. TEP Studio follows the browser's
// `color-scheme` and cannot know which it will get.
const PALETTE = [
  "#4c78a8",
  "#e45756",
  "#54a24b",
  "#b279a2",
  "#eeca3b",
  "#72b7b2",
  "#ff9d4e",
  "#9d755d",
];

/** Colour for the nth selected channel. */
export function seriesColour(n) {
  return PALETTE[n % PALETTE.length];
}

export class TrendGrid {
  /**
   * @param {HTMLElement} container the element the figures are built into
   */
  constructor(container) {
    this.container = container;
    this.labels = [];
    this.units = [];
    // Zero-based offsets of `XMEAS(23..=41)`, which hold their value between
    // analyser reports. Drawing them with straight lines between samples would
    // invent a ramp the plant never had, so they are drawn as steps.
    this.stepped = new Set();
    this.selection = [];
    /** @type {{col: number, canvas: HTMLCanvasElement, value: HTMLElement}[]} */
    this.panels = [];
    this.totalHours = 1;

    // Canvases are sized in CSS pixels by the layout and in device pixels by
    // us. Without the observer a window resize leaves every chart blurry.
    this.observer =
      typeof ResizeObserver === "undefined"
        ? null
        : new ResizeObserver(() => this.resize());
  }

  /**
   * Column metadata, from the bindings. Labels and units come from
   * `tepsim_core::variables`, which the core checks against the `teprob.f`
   * header table, so nothing on an axis here was typed by hand.
   */
  setColumns(labels, units, sampledColumns) {
    this.labels = labels;
    this.units = units;
    this.stepped = new Set(sampledColumns);
  }

  /** The run length, which is the x axis. */
  setHours(hours) {
    this.totalHours = hours > 0 ? hours : 1;
  }

  /**
   * Rebuild the grid for a set of column indices into a packed row.
   *
   * @param {number[]} columns
   */
  setSelection(columns) {
    this.selection = [...columns];
    this.container.textContent = "";
    this.panels = [];
    if (this.observer) this.observer.disconnect();

    if (this.selection.length === 0) {
      const empty = document.createElement("p");
      empty.className = "muted";
      empty.textContent = "No channels selected. Pick some from the trend list.";
      this.container.append(empty);
      return;
    }

    this.selection.forEach((col, n) => {
      const figure = document.createElement("figure");
      figure.className = "trend";

      const caption = document.createElement("figcaption");
      const swatch = document.createElement("span");
      swatch.className = "swatch";
      swatch.style.background = seriesColour(n);
      const name = document.createElement("span");
      name.className = "trend-name";
      name.textContent = this.labels[col] ?? `column ${col}`;
      const value = document.createElement("span");
      value.className = "trend-value mono";
      value.textContent = "-";
      caption.append(swatch, name, value);

      const canvas = document.createElement("canvas");
      canvas.className = "trend-canvas";

      figure.append(caption, canvas);
      this.container.append(figure);
      this.panels.push({ col, canvas, value });
      if (this.observer) this.observer.observe(canvas);
    });

    this.resize();
  }

  /** Match every canvas's backing store to its CSS size and the pixel ratio. */
  resize() {
    const ratio = globalThis.devicePixelRatio || 1;
    for (const panel of this.panels) {
      const rect = panel.canvas.getBoundingClientRect();
      const w = Math.max(1, Math.round(rect.width * ratio));
      const h = Math.max(1, Math.round(rect.height * ratio));
      if (panel.canvas.width !== w || panel.canvas.height !== h) {
        panel.canvas.width = w;
        panel.canvas.height = h;
      }
    }
  }

  /**
   * Repaint every chart from the history.
   *
   * @param {import("./history.js").History} history
   */
  draw(history) {
    // `getComputedStyle` forces a style resolution, so it is read once per
    // frame rather than once per chart. With a dozen panels at display rate
    // that is the difference between a free call and a measurable one.
    const ink = getComputedStyle(this.container).color;
    this.panels.forEach((panel, n) => {
      this.drawPanel(panel, history, seriesColour(n), ink);
    });
  }

  drawPanel(panel, history, colour, ink) {
    const ctx = panel.canvas.getContext("2d");
    if (!ctx) return;
    const w = panel.canvas.width;
    const h = panel.canvas.height;
    ctx.clearRect(0, 0, w, h);

    if (history.count === 0) {
      panel.value.textContent = "-";
      return;
    }

    const extent = history.extent(panel.col);
    if (!extent) return;
    const { lo, hi } = niceRange(extent.lo, extent.hi);
    const span = hi - lo;

    const ratio = globalThis.devicePixelRatio || 1;
    const padTop = 4 * ratio;
    const padBottom = 12 * ratio;
    const padLeft = 44 * ratio;
    const padRight = 4 * ratio;
    const plotW = Math.max(1, w - padLeft - padRight);
    const plotH = Math.max(1, h - padTop - padBottom);

    const x = (hours) => padLeft + (hours / this.totalHours) * plotW;
    const y = (v) => padTop + (1 - (v - lo) / span) * plotH;

    // Frame and two gridlines. Enough to read a value off, not so much that a
    // 90-pixel-tall chart becomes a hatch pattern.
    ctx.strokeStyle = ink;
    ctx.globalAlpha = 0.18;
    ctx.lineWidth = Math.max(1, ratio * 0.75);
    ctx.beginPath();
    for (const frac of [0, 0.5, 1]) {
      const py = Math.round(padTop + frac * plotH) + 0.5;
      ctx.moveTo(padLeft, py);
      ctx.lineTo(padLeft + plotW, py);
    }
    ctx.stroke();
    ctx.globalAlpha = 1;

    // Axis labels.
    ctx.fillStyle = ink;
    ctx.globalAlpha = 0.65;
    ctx.font = `${Math.round(9 * ratio)}px ui-monospace, SFMono-Regular, Menlo, monospace`;
    ctx.textAlign = "right";
    ctx.textBaseline = "top";
    ctx.fillText(formatValue(hi), padLeft - 4 * ratio, padTop);
    ctx.textBaseline = "bottom";
    ctx.fillText(formatValue(lo), padLeft - 4 * ratio, padTop + plotH);
    ctx.textAlign = "left";
    ctx.fillText("0", padLeft, h);
    ctx.textAlign = "right";
    ctx.fillText(`${formatHours(this.totalHours)} h`, padLeft + plotW, h);
    ctx.globalAlpha = 1;

    // The trace. Two points per pixel column, which is all a pixel column can
    // show, and bounded by the canvas rather than by the run length.
    const points = decimateMinMax(history, panel.col, Math.round(plotW));
    const columns = points.length / 3;
    const stepped = this.stepped.has(panel.col);

    ctx.strokeStyle = colour;
    ctx.lineWidth = Math.max(1, ratio * 1.1);
    ctx.lineJoin = "round";
    ctx.beginPath();
    let previousY = null;
    for (let c = 0; c < columns; c += 1) {
      const px = x(points[c * 3]);
      const yLo = y(points[c * 3 + 1]);
      const yHi = y(points[c * 3 + 2]);
      if (previousY === null) {
        ctx.moveTo(px, yHi);
      } else if (stepped) {
        // Hold the previous reading until this one, which is what the analyser
        // actually did.
        ctx.lineTo(px, previousY);
      }
      // Within a pixel column the extremes are drawn as a vertical span, so a
      // spike survives decimation.
      ctx.lineTo(px, yHi);
      if (yLo !== yHi) ctx.lineTo(px, yLo);
      previousY = yLo;
    }
    ctx.stroke();

    const last = history.lastRow();
    const value = last ? last[panel.col] : Number.NaN;
    const unit = this.units[panel.col] ?? "";
    panel.value.textContent = `${formatValue(value)} ${unit}`.trim();
  }
}
