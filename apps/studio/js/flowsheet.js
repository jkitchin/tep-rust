// The process flow diagram, with live readings on it.
//
// Downs and Vogel's plant is five unit operations and thirteen streams, and the
// thing a trend plot cannot show is which of them is upstream of which. When
// IDV(1) steps the A/C feed composition, the interesting part is not that
// XMEAS(4) moved; it is that the reactor feed composition moved, then the
// pressure, then the purge, and the recycle loop carried it round again. A
// flowsheet is the only view where that reads as one event.
//
// # What is drawn is what the source does
//
// The topology below is the Fortran's, not the paper's figure. The two disagree
// on stream numbering: `FTM(1)` is the D feed, which the paper calls stream 2,
// and `FTM(3)` is the A feed, which the paper calls stream 1. Every
// reimplementation of TEP rediscovers this, and getting it wrong produces a
// plant that runs, looks plausible, and is wired up incorrectly. The mapping is
// in `book/src/process/plant.md`, established from `teprob.f` rather than from
// the figure, and the stream labels here use the paper's numbers because those
// are what the measurement descriptions use.
//
// # Where the numbers come from
//
// A readout names a Fortran index and nothing else. The short caption beside it
// is for the diagram; the authoritative description and unit come from
// `columnLabels()` and `columnUnits()` at run time, which the core checks
// against the `teprob.f` header table, and they are attached as the hover
// title. So a readout cannot claim to be something the bindings disagree with,
// and there is no second table of units here to fall out of date.

import { formatValue } from "./format.js";

const NS = "http://www.w3.org/2000/svg";

// One entry per number on the diagram. `meas` is a one-based `XMEAS` index,
// `mv` a one-based `XMV` index; a packed row is `[hours, XMEAS(1..41),
// XMV(1..12)]`, so the column is the index itself or 41 plus it.
//
// `x` and `y` are the anchor of the value text. `caption` is the short name
// drawn to its left; the full description is the hover title.
const READOUTS = [
  // The four feeds, along the left edge.
  { meas: 2, caption: "D feed", x: 150, y: 44 },
  { mv: 1, caption: "valve", x: 150, y: 58 },
  { meas: 3, caption: "E feed", x: 150, y: 90 },
  { mv: 2, caption: "valve", x: 150, y: 104 },
  { meas: 1, caption: "A feed", x: 150, y: 136 },
  { mv: 3, caption: "valve", x: 150, y: 150 },

  // Mixing zone outlet, the reactor's inlet.
  { meas: 6, caption: "feed rate", x: 296, y: 196 },

  // Reactor.
  { meas: 7, caption: "P", x: 400, y: 150 },
  { meas: 8, caption: "level", x: 400, y: 166 },
  { meas: 9, caption: "T", x: 400, y: 182 },
  { meas: 21, caption: "cw out", x: 400, y: 198 },
  { mv: 10, caption: "cw flow", x: 400, y: 214 },
  { mv: 12, caption: "agitator", x: 400, y: 230 },

  // Compressor and the recycle it drives.
  { meas: 20, caption: "work", x: 560, y: 42 },
  { meas: 5, caption: "recycle", x: 560, y: 58 },
  { mv: 5, caption: "recycle valve", x: 560, y: 74 },

  // Purge.
  { meas: 10, caption: "purge", x: 916, y: 44 },
  { mv: 6, caption: "valve", x: 916, y: 58 },
  { meas: 30, caption: "purge B", x: 916, y: 72 },

  // Separator.
  { meas: 13, caption: "P", x: 700, y: 150 },
  { meas: 12, caption: "level", x: 700, y: 166 },
  { meas: 11, caption: "T", x: 700, y: 182 },
  { meas: 22, caption: "cw out", x: 700, y: 198 },
  { mv: 11, caption: "cw flow", x: 700, y: 214 },
  { meas: 14, caption: "underflow", x: 700, y: 230 },
  { mv: 7, caption: "valve", x: 700, y: 246 },

  // Stripper.
  { meas: 16, caption: "P", x: 916, y: 258 },
  { meas: 15, caption: "level", x: 916, y: 274 },
  { meas: 18, caption: "T", x: 916, y: 290 },
  { meas: 19, caption: "steam", x: 916, y: 306 },
  { mv: 9, caption: "steam valve", x: 916, y: 322 },

  // The A and C feed enters the stripper, not the reactor.
  { meas: 4, caption: "A+C feed", x: 660, y: 388 },
  { mv: 4, caption: "valve", x: 660, y: 402 },

  // Product, out of the stripper base.
  { meas: 17, caption: "product", x: 916, y: 360 },
  { mv: 8, caption: "valve", x: 916, y: 374 },
  { meas: 40, caption: "G", x: 916, y: 390 },
  { meas: 41, caption: "H", x: 916, y: 404 },
];

/** Column of a readout in a packed row. */
function columnOf(readout) {
  return readout.meas !== undefined ? readout.meas : 41 + readout.mv;
}

/** A stable element id for a readout, so the update loop can find it. */
function idOf(readout) {
  return readout.meas !== undefined ? `pfd-m${readout.meas}` : `pfd-v${readout.mv}`;
}

// The static drawing. Vessels, pipes, arrowheads and titles: everything that
// does not change while a run is going.
//
// Laid out left to right in the order the fluid travels: feeds, mixing zone,
// reactor, condenser, separator, stripper, product. The recycle runs back along
// the top, which is where a loop belongs on a flowsheet, and the stripper
// overhead runs back along the bottom of the top band to the same mixing zone.
const SKELETON = `
  <defs>
    <marker id="pfd-arrow" viewBox="0 0 10 10" refX="9" refY="5"
            markerWidth="6" markerHeight="6" orient="auto-start-reverse">
      <path d="M 0 0 L 10 5 L 0 10 z" class="pfd-arrowhead" />
    </marker>
  </defs>

  <!-- Feed lines into the mixing zone header. -->
  <path class="pfd-pipe" d="M 40 50 H 232" marker-end="url(#pfd-arrow)" />
  <path class="pfd-pipe" d="M 40 96 H 232" marker-end="url(#pfd-arrow)" />
  <path class="pfd-pipe" d="M 40 142 H 232" marker-end="url(#pfd-arrow)" />
  <text class="pfd-stream" x="44" y="40">stream 2</text>
  <text class="pfd-stream" x="44" y="86">stream 3</text>
  <text class="pfd-stream" x="44" y="132">stream 1</text>

  <!-- Mixing zone. -->
  <rect class="pfd-vessel" x="232" y="30" width="46" height="130" rx="6" />
  <text class="pfd-unit" x="255" y="98" text-anchor="middle">mixing</text>
  <text class="pfd-unit" x="255" y="112" text-anchor="middle">zone</text>

  <!-- Mixing zone to reactor, stream 6. -->
  <path class="pfd-pipe" d="M 278 95 H 300 V 190 H 322" marker-end="url(#pfd-arrow)" />
  <text class="pfd-stream" x="302" y="186">stream 6</text>

  <!-- Reactor. -->
  <rect class="pfd-vessel" x="322" y="120" width="120" height="150" rx="10" />
  <text class="pfd-unit" x="382" y="112" text-anchor="middle">reactor</text>
  <line class="pfd-level" x1="322" y1="210" x2="442" y2="210" />
  <path class="pfd-cooling" d="M 332 282 H 432" marker-end="url(#pfd-arrow)" />
  <text class="pfd-stream" x="332" y="296">cooling water</text>

  <!-- Reactor to condenser to separator, stream 7. -->
  <path class="pfd-pipe" d="M 442 150 H 480" marker-end="url(#pfd-arrow)" />
  <rect class="pfd-vessel" x="480" y="132" width="56" height="36" rx="6" />
  <text class="pfd-unit" x="508" y="154" text-anchor="middle">cond</text>
  <path class="pfd-pipe" d="M 536 150 H 596" marker-end="url(#pfd-arrow)" />
  <text class="pfd-stream" x="446" y="144">stream 7</text>

  <!-- Separator. -->
  <rect class="pfd-vessel" x="596" y="120" width="110" height="150" rx="10" />
  <text class="pfd-unit" x="651" y="112" text-anchor="middle">separator</text>
  <line class="pfd-level" x1="596" y1="216" x2="706" y2="216" />
  <path class="pfd-cooling" d="M 606 282 H 696" marker-end="url(#pfd-arrow)" />
  <text class="pfd-stream" x="606" y="296">cooling water</text>

  <!-- Separator vapour to the compressor and back round as the recycle. -->
  <path class="pfd-pipe" d="M 651 120 V 60 H 522" marker-end="url(#pfd-arrow)" />
  <circle class="pfd-vessel" cx="500" cy="60" r="20" />
  <text class="pfd-unit" x="500" y="64" text-anchor="middle">C</text>
  <path class="pfd-pipe" d="M 480 60 H 255 V 30" marker-end="url(#pfd-arrow)" />
  <text class="pfd-stream" x="560" y="94">stream 8, recycle</text>

  <!-- Purge, off the separator vapour line. -->
  <path class="pfd-pipe" d="M 720 60 H 846" marker-end="url(#pfd-arrow)" />
  <path class="pfd-pipe" d="M 651 60 H 720" />
  <text class="pfd-stream" x="770" y="52">stream 9, purge</text>

  <!-- Separator underflow to the stripper, stream 10. -->
  <path class="pfd-pipe" d="M 706 240 H 786" marker-end="url(#pfd-arrow)" />
  <text class="pfd-stream" x="712" y="234">stream 10</text>

  <!-- Stripper. -->
  <rect class="pfd-vessel" x="786" y="200" width="100" height="150" rx="10" />
  <text class="pfd-unit" x="836" y="192" text-anchor="middle">stripper</text>
  <line class="pfd-level" x1="786" y1="290" x2="886" y2="290" />

  <!-- The A and C feed enters the stripper base as the stripping gas. -->
  <path class="pfd-pipe" d="M 700 394 H 836 V 350" marker-end="url(#pfd-arrow)" />
  <text class="pfd-stream" x="660" y="374">stream 4</text>

  <!-- Stripper overhead back to the mixing zone, stream 5. -->
  <path class="pfd-pipe" d="M 786 214 H 760 V 176 H 255 V 160" marker-end="url(#pfd-arrow)" />
  <text class="pfd-stream" x="560" y="172">stream 5, stripper overhead</text>

  <!-- Product out of the stripper base, stream 11. -->
  <path class="pfd-pipe" d="M 886 340 H 946" marker-end="url(#pfd-arrow)" />
  <text class="pfd-stream" x="890" y="334">stream 11</text>

  <!-- Steam to the stripper reboiler. -->
  <path class="pfd-cooling" d="M 786 320 H 700" marker-start="url(#pfd-arrow)" />
  <text class="pfd-stream" x="700" y="334">steam</text>
`;

/**
 * Build the diagram into an `<svg>` and return a handle that can update it.
 *
 * @param {SVGSVGElement} svg
 * @param {string[]} labels from `columnLabels()`
 * @param {string[]} units from `columnUnits()`
 */
export function buildFlowsheet(svg, labels, units) {
  svg.innerHTML = SKELETON;

  /** @type {{el: SVGTextElement, col: number, unit: string}[]} */
  const bound = [];

  for (const readout of READOUTS) {
    const col = columnOf(readout);
    const group = document.createElementNS(NS, "g");
    group.classList.add("pfd-readout");

    // The full description, from the bindings. A diagram short of space still
    // has room for a tooltip, and this is where the honest label lives.
    const title = document.createElementNS(NS, "title");
    title.textContent = labels[col] ?? `column ${col}`;
    group.append(title);

    const caption = document.createElementNS(NS, "text");
    caption.setAttribute("x", String(readout.x));
    caption.setAttribute("y", String(readout.y));
    caption.setAttribute("text-anchor", "end");
    caption.classList.add("pfd-caption");
    caption.textContent = readout.caption;

    const value = document.createElementNS(NS, "text");
    value.setAttribute("id", idOf(readout));
    value.setAttribute("x", String(readout.x + 6));
    value.setAttribute("y", String(readout.y));
    value.classList.add("pfd-value");
    value.textContent = "-";

    group.append(caption, value);
    svg.append(group);
    bound.push({ el: value, col, unit: units[col] ?? "" });
  }

  return {
    /**
     * Write the latest row onto the diagram.
     *
     * @param {Float64Array | null} row a packed row, or null to blank it
     */
    update(row) {
      for (const { el, col, unit } of bound) {
        el.textContent = row ? `${formatValue(row[col])} ${unit}`.trim() : "-";
      }
    },

    /** Every column the diagram reads, for tests and for sanity. */
    columns: bound.map((b) => b.col),
  };
}

// Exported for the Node tests, which check the binding table against the
// bindings' own column count rather than against a second copy of it here.
export { READOUTS, columnOf, idOf };
