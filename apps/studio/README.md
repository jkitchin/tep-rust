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
history store, and a static pass over the built page. Rendering is not tested
and is not pretended to be.

## What it does

A run streams from the worker as packed rows, `[hours, XMEAS(1..41),
XMV(1..12)]`, arriving as transferable `Float64Array`s. On the page they feed
three views: a process flow diagram with live readings on every unit and stream,
a grid of trend charts over a selectable subset of the 53 channels, and a panel
of the twenty `IDV` disturbances. The scenario lives in the URL fragment, so a
run travels as a link.

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
"abort"`, stripped), no `wasm-opt`:

| | raw | gzip -9 |
|---|---|---|
| wasm module | 87,298 | 35,359 |
| wasm-bindgen glue | 50,309 | 9,627 |
| the app itself (HTML, CSS, 7 modules) | 79,121 | 29,024 |
| **total** | **216,728** | **74,010** |

`PLAN.org` budgets under 1.5 MB gzipped for the whole app. This is 4.8 percent
of it. `apps/studio/measure.sh` prints the table and fails if the budget is
exceeded.

`wasm-opt -Oz` takes the module from 87,298 to 79,759 bytes raw and from 35,359
to 35,284 gzipped: worth having, worth almost nothing compressed, and not worth
failing a build over.

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

Nothing visual. The build succeeds, the protocol is driven under Node, the
static structure of the built page is checked, and the numbers are asserted, but
no browser has rendered this page. Layout, the flowsheet's appearance at real
widths, canvas rendering, `ResizeObserver` behaviour, the clipboard fallback and
the dark-mode palette are all unexercised.
