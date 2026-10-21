//! Optional adapters for external math libraries.
//!
//! The `nalgebra` feature implements `Point` for nalgebra vectors.

#[cfg(feature = "nalgebra")]
pub mod nalgebra;
