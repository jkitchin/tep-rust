//! Disturbances, events, and reproducible scenarios.
//!
//! Separates *where* a disturbance couples into the plant (`InjectionPoint`)
//! from *how* it evolves in time (`Profile`). The twenty canonical IDV faults
//! then become twenty entries in a table rather than twenty special cases, and
//! user-defined faults, continuous magnitudes, and composition all follow.
//!
//! A `Scenario` is the serializable unit of reproducible work: seed, duration,
//! integrator, control scheme, and a schedule of events, with a content hash so
//! a generated dataset is self-describing.
//!
//! # Status
//!
//! Skeleton. Lands in phases 3 and 6; see `BACKLOG.org`.

#![forbid(unsafe_code)]
