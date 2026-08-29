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
// The one that catches everybody: **stream 4, the mixed A and C feed, does not
// enter the mixing zone.** `teprob.f:637-639` adds it to the stripper base,
// where it is the stripping gas. It is drawn entering the stripper from below
// for that reason, and the note beside it says so.
//
// # Two bands, and no pipe through a vessel
//
// The drawing is two horizontal bands. The top band is the vapour loop in the
// order the fluid travels: feeds, mixing zone, reactor, condenser, separator,
// with the compressor recycle and the purge running back along the header at
// the top. The bottom band is the stripper, which the separator underflow
// reaches by going round the right-hand side and which returns its overhead to
// the mixing zone along a lane between the two bands.
//
// That lane exists so the return has somewhere to run. A flowsheet whose pipes
// cross its vessels is not merely ugly: a line that appears to pass through the
// separator is a claim about the plant, and it is false. Every pipe here is
// orthogonal and every one of them stays in open space, which is checked by
// eye at each change and is the reason the layout is as tall as it is.
//
// # Where the numbers come from, and which numbers
//
// A readout names a Fortran index and nothing else. The short caption beside it
// is for the diagram; the authoritative description and unit come from
// `columnLabels()` and `columnUnits()` at run time, which the core checks
// against the `teprob.f` header table, and they are attached as the hover
// title. So a readout cannot claim to be something the bindings disagree with,
// and there is no second table of units here to fall out of date.
//
// The diagram carries all 22 continuous measurements and exactly three of the
// twelve valves. The rule is that a valve earns a place only where nothing else
// reports what it does: `XMV(5)` the compressor recycle, `XMV(10)` the reactor
// cooling water and `XMV(11)` the condenser cooling water move flows that no
// `XMEAS` measures. The other nine each sit directly upstream of a measurement
// that is already drawn (`XMV(1)` sets the D feed, and `XMEAS(2)` is the D
// feed), except `XMV(12)`, the agitator, which moves no material and holds at
// 50% in every published run. Space on a diagram is not free, and a number
// nobody reads costs the readable ones the room they needed. All twelve are on
// the trends, at full resolution, and all 53 channels are in the CSV.

import { formatReadout, shortUnit } from "./format.js";

const NS = "http://www.w3.org/2000/svg";

// One entry per number on the diagram. `meas` is a one-based `XMEAS` index,
// `mv` a one-based `XMV` index; a packed row is `[hours, XMEAS(1..41),
// XMV(1..12)]`, so the column is the index itself or 41 plus it.
//
// `x` and `y` anchor the pair: the caption ends at `x`, the value starts just
// after it. Readouts are grouped into blocks that sit in the open space beside
// the unit they belong to, at a 18-point pitch, and the block is placed so that
// the longest value it can ever show still clears whatever is next to it. The
// widest is a feed in `kscmh`, which is "-0.2505 kscmh" if a lost feed goes
// negative under noise: thirteen characters, about 94 units of width.
const READOUTS = [
  // The three feeds that reach the mixing zone, each on its own line in the
  // left margin. The fourth feed is at the stripper, below.
  { meas: 1, caption: "A feed", x: 70, y: 161 },
  { meas: 2, caption: "D feed", x: 70, y: 206 },
  { meas: 3, caption: "E feed", x: 70, y: 251 },

  // Compressor and the recycle it drives, under the top header.
  { meas: 20, caption: "work", x: 340, y: 82 },
  { meas: 5, caption: "recycle", x: 340, y: 100 },
  { mv: 5, caption: "recycle valve", x: 340, y: 118 },

  // Purge, at the end of the header. The B reading is the reason the purge
  // exists: it is the only exit for the inert.
  { meas: 10, caption: "purge", x: 1010, y: 62 },
  { meas: 30, caption: "purge B", x: 1010, y: 80 },

  // Reactor. The first row is the mixing zone outlet, which is the reactor's
  // inlet, and the last two are the cooling water that sets the temperature.
  { meas: 6, caption: "feed rate", x: 370, y: 352 },
  { meas: 7, caption: "P", x: 370, y: 370 },
  { meas: 8, caption: "level", x: 370, y: 388 },
  { meas: 9, caption: "T", x: 370, y: 406 },
  { meas: 21, caption: "cw out", x: 370, y: 424 },
  { mv: 10, caption: "cw valve", x: 370, y: 442 },

  // Separator, with the condenser cooling water that feeds it. XMEAS(22) is
  // "Separator Cooling Water Outlet Temp" and XMV(11) is "Condenser Cooling
  // Water Flow": one circuit, and the hover titles carry both names.
  { meas: 13, caption: "P", x: 700, y: 352 },
  { meas: 12, caption: "level", x: 700, y: 370 },
  { meas: 11, caption: "T", x: 700, y: 388 },
  { meas: 14, caption: "underflow", x: 700, y: 406 },
  { meas: 22, caption: "cw out", x: 700, y: 424 },
  { mv: 11, caption: "cw valve", x: 700, y: 442 },

  // Stripper.
  { meas: 16, caption: "P", x: 360, y: 520 },
  { meas: 15, caption: "level", x: 360, y: 538 },
  { meas: 18, caption: "T", x: 360, y: 556 },
  { meas: 19, caption: "steam", x: 360, y: 574 },

  // The A and C feed, entering the stripper base.
  { meas: 4, caption: "A+C feed", x: 300, y: 678 },

  // Product, out of the stripper base, with the two components it is specified
  // on.
  { meas: 17, caption: "product", x: 950, y: 618 },
  { meas: 40, caption: "G", x: 950, y: 640 },
  { meas: 41, caption: "H", x: 950, y: 662 },
];

/** Column of a readout in a packed row. */
function columnOf(readout) {
  return readout.meas !== undefined ? readout.meas : 41 + readout.mv;
}

/** A stable element id for a readout, so the update loop can find it. */
function idOf(readout) {
  return readout.meas !== undefined ? `pfd-m${readout.meas}` : `pfd-v${readout.mv}`;
}

// The static drawing: vessels, pipes, arrowheads and titles, everything that
// does not change while a run is going.
//
// The coordinate system is the 1200 by 720 viewBox in index.html. Vessels
// occupy y 140-330 in the top band and y 490-660 in the bottom one; the two
// lanes at y 46 and y 478 carry the recycle and the stripper overhead. Numbers
// live in the gaps that leaves, which is why they are where they are.
const SKELETON = `
  <defs>
    <marker id="pfd-arrow" viewBox="0 0 10 10" refX="9" refY="5"
            markerWidth="7" markerHeight="7" markerUnits="userSpaceOnUse"
            orient="auto-start-reverse">
      <path d="M 0 0 L 10 5 L 0 10 z" class="pfd-arrowhead" />
    </marker>
  </defs>

  <!-- The three feeds that enter the mixing zone. -->
  <path class="pfd-pipe" d="M 24 170 H 190" marker-end="url(#pfd-arrow)" />
  <path class="pfd-pipe" d="M 24 215 H 190" marker-end="url(#pfd-arrow)" />
  <path class="pfd-pipe" d="M 24 260 H 190" marker-end="url(#pfd-arrow)" />
  <text class="pfd-stream" x="24" y="185">stream 1</text>
  <text class="pfd-stream" x="24" y="230">stream 2</text>
  <text class="pfd-stream" x="24" y="275">stream 3</text>

  <!-- Mixing zone. Five things enter it: three feeds on the left, the recycle
       from the top, the stripper overhead from below left. -->
  <rect class="pfd-vessel" x="190" y="140" width="76" height="190" rx="6" />
  <text class="pfd-unit" x="228" y="292" text-anchor="middle">mixing</text>
  <text class="pfd-unit" x="228" y="308" text-anchor="middle">zone</text>

  <!-- Mixing zone to reactor, stream 6. -->
  <path class="pfd-pipe" d="M 266 235 H 330" marker-end="url(#pfd-arrow)" />
  <text class="pfd-stream" x="298" y="227" text-anchor="middle">stream 6</text>

  <!-- Reactor, with its level and its cooling coil. -->
  <rect class="pfd-vessel" x="330" y="140" width="130" height="190" rx="10" />
  <text class="pfd-unit" x="395" y="164" text-anchor="middle">reactor</text>
  <line class="pfd-level" x1="330" y1="250" x2="460" y2="250" />
  <path class="pfd-cooling" d="M 296 290 H 420 V 312 H 296" marker-end="url(#pfd-arrow)" />
  <text class="pfd-stream" x="288" y="305" text-anchor="end">cw</text>

  <!-- Reactor through the condenser to the separator, stream 7. -->
  <path class="pfd-pipe" d="M 460 180 H 500" marker-end="url(#pfd-arrow)" />
  <rect class="pfd-vessel" x="500" y="160" width="80" height="40" rx="6" />
  <text class="pfd-unit" x="540" y="185" text-anchor="middle">cond</text>
  <path class="pfd-pipe" d="M 580 180 H 620" marker-end="url(#pfd-arrow)" />
  <text class="pfd-stream" x="540" y="148" text-anchor="middle">stream 7</text>
  <path class="pfd-cooling" d="M 540 246 V 200" marker-end="url(#pfd-arrow)" />
  <text class="pfd-stream" x="548" y="230">cw</text>

  <!-- Separator. -->
  <rect class="pfd-vessel" x="620" y="140" width="120" height="190" rx="10" />
  <text class="pfd-unit" x="680" y="164" text-anchor="middle">separator</text>
  <line class="pfd-level" x1="620" y1="250" x2="740" y2="250" />

  <!-- The top header: separator vapour up, then west through the compressor
       and back to the mixing zone, and east to the purge. -->
  <path class="pfd-pipe" d="M 680 140 V 46" />
  <path class="pfd-pipe" d="M 680 46 H 452" marker-end="url(#pfd-arrow)" />
  <circle class="pfd-vessel" cx="430" cy="46" r="22" />
  <text class="pfd-unit" x="430" y="51" text-anchor="middle">K</text>
  <text class="pfd-unit" x="430" y="20" text-anchor="middle">compressor</text>
  <path class="pfd-pipe" d="M 408 46 H 228 V 140" marker-end="url(#pfd-arrow)" />
  <text class="pfd-stream" x="300" y="36" text-anchor="middle">stream 8, recycle</text>
  <path class="pfd-pipe" d="M 680 46 H 900" marker-end="url(#pfd-arrow)" />
  <text class="pfd-stream" x="790" y="36" text-anchor="middle">stream 9, purge</text>

  <!-- Separator underflow to the stripper, round the outside. -->
  <path class="pfd-pipe" d="M 740 280 H 1000 V 560 H 590" marker-end="url(#pfd-arrow)" />
  <text class="pfd-stream" x="1008" y="400">stream 10</text>

  <!-- Stripper. -->
  <rect class="pfd-vessel" x="470" y="490" width="120" height="170" rx="10" />
  <text class="pfd-unit" x="530" y="514" text-anchor="middle">stripper</text>
  <line class="pfd-level" x1="470" y1="590" x2="590" y2="590" />
  <path class="pfd-cooling" d="M 330 620 H 470" marker-end="url(#pfd-arrow)" />
  <text class="pfd-stream" x="322" y="624" text-anchor="end">steam</text>

  <!-- Stripper overhead back to the mixing zone, along the lane between the
       two bands. This is the pipe that must not cross anything. -->
  <path class="pfd-pipe" d="M 510 490 V 478 H 150 V 300 H 190"
        marker-end="url(#pfd-arrow)" />
  <text class="pfd-stream" x="330" y="470" text-anchor="middle">stream 5, stripper overhead</text>

  <!-- Product out of the stripper base, stream 11. -->
  <path class="pfd-pipe" d="M 590 630 H 840" marker-end="url(#pfd-arrow)" />
  <text class="pfd-stream" x="700" y="622" text-anchor="middle">stream 11, product</text>

  <!-- The A and C feed enters the stripper base as the stripping gas. -->
  <path class="pfd-pipe" d="M 200 690 H 490 V 660" marker-end="url(#pfd-arrow)" />
  <text class="pfd-stream" x="200" y="708">stream 4, A and C feed, to the stripper base</text>
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

    // The full description and the unit as the Fortran spells it. A diagram
    // short of space still has room for a tooltip, and this is where the
    // unabbreviated version lives.
    const title = document.createElementNS(NS, "title");
    const label = labels[col] ?? `column ${col}`;
    const unit = units[col] ?? "";
    title.textContent = unit ? `${label} (${unit})` : label;
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
    bound.push({ el: value, col, unit: shortUnit(unit) });
  }

  return {
    /**
     * Write the latest row onto the diagram.
     *
     * @param {Float64Array | null} row a packed row, or null to blank it
     */
    update(row) {
      for (const { el, col, unit } of bound) {
        el.textContent = row ? `${formatReadout(row[col])} ${unit}`.trim() : "-";
      }
    },

    /** Every column the diagram reads, for tests and for sanity. */
    columns: bound.map((b) => b.col),
  };
}

// Exported for the Node tests, which check the binding table against the
// bindings' own column count rather than against a second copy of it here.
export { READOUTS, columnOf, idOf };
