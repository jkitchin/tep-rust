# TEP Studio

The Tennessee Eastman Process, running in a browser tab. The plant is
`crates/tepsim-wasm` compiled to WebAssembly and driven on a Web Worker; this
directory is the interface to it.

## Run it

```sh
apps/studio/build.sh
python3 -m http.server 8000 --directory apps/studio/dist
open http://localhost:8000/
```

A module worker cannot be loaded from a `file://` URL, so the page has to be
served. Any static file server will do; `dist/` is the whole deployment.

Prerequisites are `rustup target add wasm32-unknown-unknown` and a
`wasm-bindgen-cli` whose version matches the `wasm-bindgen` dependency in
`Cargo.lock` exactly. `build.sh` reads the required version out of `Cargo.lock`
and tells you the command if it does not match. `wasm-opt` is used when it is on
`PATH` and skipped with a note when it is not.

## Test it

```sh
node --test --test-force-exit 'apps/studio/node/*.test.mjs'
```

Everything that can be checked without a browser is checked there: the worker
protocol driven under a shim of the Worker globals, the wasm module's
determinism digest against the value `crates/tepsim-wasm/tests/determinism.rs`
pins, throughput against the budget in `PLAN.org`, the link encoding, the
history store, the CSV export down to and including revoking the object URL,
and a static pass over the built page. Rendering is not tested and is not
pretended to be.

## What it does

A run streams from the worker as packed rows, `[hours, XMEAS(1..41),
XMV(1..12)]`, arriving as transferable `Float64Array`s. On the page they feed
three views: a process flow diagram with live readings on every unit and stream,
a grid of trend charts over a selectable subset of the 53 channels, and a panel
of the twenty `IDV` disturbances. The scenario lives in the URL fragment, so a
run travels as a link, and **Download CSV** saves what the run produced.

It travels as one token, `#q=tepsim.scenario.v1;seed=...`, which is what
`Scenario.text` produces and `Scenario.fromText` reads. Neither the page, the
worker nor `share.js` knows what is inside it. Before B-0054a all three
enumerated the fields of a scenario by hand, so a field added to `Scenario`
reached none of them and nothing failed: the link still opened and ran a
different scenario from the one it named. The format is versioned and its parser
is strict, so a link this build cannot honour is refused with a message rather
than quietly opening the baseline. See `crates/tepsim/src/text.rs`.

The first thing the page does is compare its own determinism digest against the
value a native run of the same commit produces. If they disagree, it says so at
the top and tells you not to trust anything below, because at that point the
invariant the whole validation ladder rests on has broken here.

## Decisions worth knowing about

**No `SharedArrayBuffer`.** It would let the worker and the page share memory
outright, and it requires COOP and COEP response headers. Neither GitHub Pages
nor a Hugging Face Static Space can set those, and free static hosting is the
reason this app exists. Chunks are transferred instead, which costs a pointer
move per chunk and nothing else.

**Toggling a disturbance restarts the run.** `tepsim::Simulation` takes its
scenario at construction and hands the disturbance vector to the driver there.
The bindings found no seam for changing it mid-run and declined to invent one,
because a run reachable only through the browser is a run nobody can reproduce.
`Sim::setFault` therefore records the request, `pendingRestart` goes true, and
the worker rebuilds from step zero. At the throughput measured below that is
fast enough to feel immediate, and `apps/studio/node/protocol.test.mjs` asserts
the rebuilt run is byte-identical to one started with the fault set.

**No third-party JavaScript, and no uPlot.** `PLAN.org` names uPlot, and for the
app it describes (dozens of series, cursor sync, pan and zoom) that is the right
call. This grid is not that app: every chart shares one time axis fixed at the
run length, has one series, and never pans or zooms. That came to about a
hundred and fifty lines in `js/chart.js`, against a vendored dependency to keep
current and its licence to carry into the bundle. The moment the page wants a
shared cursor, zoom, or forty series on one axis, vendor uPlot and delete that
file; the `TrendGrid` surface is three methods for exactly that reason.

**Trends are drawn from a decimating store.** `sampleEvery = 1` on a 48-hour run
is 172,800 samples, and the bindings will accept up to a million. `js/history.js`
keeps 20,000 rows and, when full, drops every other row and doubles its stride,
so the store always spans the whole run and loses only resolution. The page says
which it is doing in the "retained" readout.

**Pacing never touches the numbers.** The speed control gates when the worker
asks for the next chunk, never what is in it. A paced run and an unpaced one
emit identical bytes, which is asserted rather than asserted-in-a-comment.

**The download is the run, not the picture of it.** `js/csv.js` keeps a second
store beside the trend history, and that one keeps every sample. Decimating a
file would silently answer a different question from the one it was asked and
nothing downstream could tell, so instead the recorder has a 200,000 row ceiling
(86 MB retained, enough for 48 hours at `sampleEvery = 1`) and *stops* at it
rather than thinning what it has. A truncated file is a prefix of the run at
full resolution, and both the "recorded" readout and the file's own header say
so. It costs nothing to keep: the chunk arrives as a transferred
`Float64Array` the page already owns, `History` copies what it wants out of it,
and the recorder holds the same array.

**The scenario travels in the file, not in its name.** The first thing after
the title is `# scenario: tepsim.scenario.v1;...`, the token `Scenario.fromText`
reads, so a file in a downloads folder reproduces its own run rather than
describing it. Not the file name: the token contains `;` and `:`, which are
illegal in a Windows file name and which browsers rewrite on the download path,
so a name would arrive mangled and stop round-tripping without saying it had. A
comment also survives the file being renamed. The name carries what somebody
scanning a folder needs instead: `tep-idv6-5.8h-116rows-tripped-<checksum>.csv`.

**Ground truth is resolved per row, not per chunk.** The worker reports which
disturbances are on once per chunk, as of the chunk's last sample, and a chunk
can be hours of plant; handing that to every row would move an onset earlier by
up to a chunk, which is the exact quantity a detection-delay figure measures. So
the *age* of each active disturbance travels with the chunk and each row is
decided against its own clock. Ages rather than absolute onset times on purpose:
the driver forces `IDV(12)` on at eight hours, which is a sampling instant, and
converting between the two would put that row on the wrong side of the
comparison by one ulp.

## On Trunk

`PLAN.org` specifies a Leptos app built with Trunk. What was built instead is
vanilla JavaScript over the existing `tepsim-wasm` bindings, and the reason is
the byte budget: a Leptos front end is a second Rust wasm module on top of the
one that does the arithmetic, for a page whose interactive surface is a form,
a table and some canvases.

Trunk 0.21.14 was installed and does work on this crate. It builds the
`tepsim-wasm` cdylib and generates the glue. Two things stopped it being worth
adopting. Its bundled `wasm-opt` (binaryen 123) rejects the module outright,
because rustc 1.97 emits `memory.copy` and that binaryen build needs
`--enable-bulk-memory-opt` to accept it. And with the front end in JavaScript,
everything Trunk would do for us is `cargo build`, `wasm-bindgen` and
`wasm-opt`, which is what `build.sh` does in three commands, without a build
dependency that downloads its own toolchain at build time. If the app is ever
rewritten in Leptos, Trunk is the right tool and this note is the record of what
to expect from it.

## Measured

Built with `--profile release-wasm` (`opt-level = "s"`, thin LTO, `panic =
"abort"`, stripped), no `wasm-opt`, as of B-0073:

| | raw | gzip -9 |
|---|---|---|
| wasm module | 165,312 | 73,579 |
| wasm-bindgen glue | 52,791 | 10,110 |
| the app itself (HTML, CSS, 8 modules) | 102,460 | 37,927 |
| **total** | **320,563** | **121,616** |

`PLAN.org` budgets under 1.5 MB gzipped for the whole app. This is 7.9 percent
of it. `apps/studio/measure.sh` prints the table and fails if the budget is
exceeded. The CSV export (`js/csv.js` and its wiring) cost 7,481 gzipped bytes,
0.5 percent of the budget: 114,135 before, 121,616 after.

`wasm-opt -Oz` took the module from 87,298 to 79,759 bytes raw and from 35,359
to 35,284 gzipped when that was measured, on a smaller module than the one
above: worth having, worth almost nothing compressed, and not worth failing a
build over.

Throughput, on an Apple silicon laptop under Node 25, five consecutive 48-hour
closed-loop runs (172,800 one-second Euler steps each) after a warm-up:

```
48 h in 2210 ms  ->  78,193x real time
48 h in 2081 ms  ->  83,057x real time
48 h in 4250 ms  ->  40,662x real time
48 h in 2966 ms  ->  58,258x real time
48 h in 2754 ms  ->  62,748x real time
```

All five produced checksum `6ee4409dc3b1fe0f`. So a 48-hour run lands between
two and four and a half seconds, four hundred to eight hundred times the 100x
budget in `PLAN.org`; the spread is machine load rather than anything in the
simulator. `apps/studio/node/protocol.test.mjs` measures the same thing on a
24-hour run and fails below 100x.

## What has not been verified

Nothing visual, by the test suite. The build succeeds, the protocol is driven
under Node, the static structure of the built page is checked, and the numbers
are asserted, but no test renders this page. Layout, canvas rendering,
`ResizeObserver` behaviour, the clipboard fallback and the dark-mode palette are
all unexercised by anything that runs in CI.

The flowsheet is the exception, and only as a one-off: its layout was laid out
against a rendered page rather than against its own source. Headless Chrome
drew it in both themes, on the baseline plant and on a run that trips, and a
script walked the rendered SVG asking whether any label sat outside the
viewBox, touched another label, crossed a vessel border or was crossed by a
pipe, including a pass with every readout forced to the widest string the
formatter can produce. That is how the overlapping, clipped version was found:
it read perfectly well as source. Anyone moving a readout should do the same,
because the geometry that matters is the one the browser computes.

The download is checked as far as it can be off a browser: the file is built,
parsed back and compared against the run that produced it, and `saveCsv` is
driven against a stand-in for `document` and `URL` so that the blob, the file
name, the click and the revocation are all asserted. What no test here covers is
what a real browser does with the click, which is the one step that actually
writes the file.
