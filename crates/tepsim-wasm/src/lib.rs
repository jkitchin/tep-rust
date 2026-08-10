//! WebAssembly bindings for the browser application.
//!
//! Empty on purpose: `wasm-bindgen` arrives in phase 8, alongside the Leptos
//! client in `apps/studio`. The simulation runs in a Web Worker and communicates
//! by `postMessage` with transferable buffers, deliberately not
//! `SharedArrayBuffer`, which would require COOP/COEP headers that static hosts
//! cannot set.
//!
//! # Status
//!
//! Skeleton. Lands in phase 8; see `BACKLOG.org`.

#![forbid(unsafe_code)]
