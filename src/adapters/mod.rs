//! Optional adapters for external math libraries.
//!
//! Enable feature flags (e.g. `nalgebra`) to add `Point` impls for
//! external vector types.

#[cfg(feature = "nalgebra")]
pub mod nalgebra;
