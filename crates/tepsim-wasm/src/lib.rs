//! WebAssembly bindings for the Tennessee Eastman Process simulator.
//!
//! Compiles to a `cdylib` a browser loads, and exposes the [`tepsim`] facade to
//! JavaScript: build a [`Scenario`], construct a [`Sim`], step it in chunks,
//! and read the samples back as `Float64Array`s. Alongside that, the twenty
//! `IDV` disturbances and the 53 channel names, so a user interface can build
//! its whole layout before it has a single number.
//!
//! # The shape it is built for
//!
//! One `Sim` per Web Worker. The worker loops: call [`Sim::step_chunk`], post
//! the array with its `ArrayBuffer` in the transfer list, yield to the event
//! loop, repeat. Transferring rather than cloning makes a chunk cost a pointer
//! move instead of a copy, and yielding between chunks is what lets the worker
//! answer "stop", "toggle fault 6" or "change speed" while a run is in flight.
//! A worker that instead called `tepsim::Simulation::run` would be unavailable
//! for all 172,800 steps of a 48-hour run.
//!
//! `SharedArrayBuffer` would remove even the pointer move and is deliberately
//! not used: it requires COOP and COEP response headers, which neither GitHub
//! Pages nor a Hugging Face Static Space can set, and free static hosting is the
//! reason the browser app exists at all. See `PLAN.org`, "The browser
//! application".
//!
//! `crates/tepsim-wasm/www/` holds a self-contained harness, an HTML page and a
//! worker, that does exactly this. No build step and no bundler; the comment at
//! the top of `www/index.html` has the one command that generates the glue.
//!
//! # Row layout
//!
//! A packed row is `[hours, XMEAS(1..41), XMV(1..12)]`, [`ROW_WIDTH`] values
//! wide, and a chunk is those rows end to end. The 53 channels are
//! `tepsim::run::Sample::row`'s layout, which is what the correlation matrix,
//! the Tier 6 detectors and the published files all use. [`column_ids`] and
//! [`column_labels`] name them.
//!
//! # Determinism
//!
//! Cross-platform bit-identical output is a hard invariant (`PLAN.org`,
//! "Numerics and determinism"; Tier 9). These bindings are built not to weaken
//! it. They contain no `f32`, no `js_sys::Math`, no `Date`, no source of
//! randomness, and no floating-point arithmetic of their own except the one
//! division that reports progress. [`Sim::checksum`] is integer-only over IEEE
//! 754 bit patterns, so the digest cannot itself disagree between
//! architectures.
//!
//! [`self_check_digest`] runs a fixed one-hour baseline scenario and returns
//! its digest, and [`runner::tepsim_wasm_self_check_digest`] exports the same
//! number to any WebAssembly runtime without needing `wasm-bindgen` glue.
//! `tests/determinism.rs` pins the value; comparing a browser's against it is
//! the wasm half of Tier 9.
//!
//! # `unsafe`
//!
//! The crate cannot `forbid(unsafe_code)`, because `#[wasm_bindgen]` expands to
//! the `extern` declarations and exported symbols that are the whole point of
//! the boundary, and because the glue-free self-check export needs
//! `#[unsafe(no_mangle)]`. It writes no `unsafe` blocks of its own. The
//! `Float64Array`s handed to JavaScript are built with `Float64Array::from`,
//! which copies out of the wasm heap into a buffer JavaScript owns; the
//! zero-copy `Float64Array::view` is `unsafe` for a good reason, since the view
//! aliases the heap and dangles the moment an allocation grows it, and a
//! transferred buffer must be JavaScript-owned anyway.

// `#[wasm_bindgen]` generates the `unsafe extern` blocks that make up the
// JavaScript boundary, so `forbid(unsafe_code)` is not available here. This is
// the next-strongest thing, and it is what CLAUDE.md requires of the three FFI
// crates: no implicit unsafe operations inside an `unsafe fn` body.
#![deny(unsafe_op_in_unsafe_fn)]

pub mod channels;
pub mod digest;
pub mod runner;
pub mod sim;

pub use channels::{column_ids, column_labels, column_units};
pub use digest::Fnv1a64;
pub use runner::{ConfigError, Plan, ROW_WIDTH, Runner, self_check_digest};
pub use sim::{Fault, Scenario, Sim};
