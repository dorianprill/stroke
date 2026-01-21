#![no_std]
#![forbid(unsafe_code)]
#![allow(incomplete_features)]
// this is needed to use expressions in const generics such as N-1 (see curve derivatives)
#![feature(generic_const_exprs)]

//! A const-generic, no-std spline library for Bezier and B-spline curves.
//!
//! The core abstraction is the [`Point`](crate::point::Point) trait, which models
//! a fixed-dimension vector space element. Optional capability traits such as
//! [`PointIndex`](crate::point::PointIndex), [`PointDot`](crate::point::PointDot),
//! and [`PointNorm`](crate::point::PointNorm) enable component access and
//! geometric helpers when required.
//!
//! # Examples
//! ```rust
//! use stroke::{Bezier, PointN};
//!
//! let curve = Bezier::<PointN<f32, 2>, 3>::new([
//!     PointN::new([0.0, 0.0]),
//!     PointN::new([1.0, 0.0]),
//!     PointN::new([1.0, 1.0]),
//! ]);
//!
//! let mid = curve.eval(0.5);
//! # let _ = mid;
//! ```
//!
//! # Feature flags
//! - `nalgebra`: implements `Point` for `nalgebra::SVector<T, D>` (add `nalgebra` as a dependency).

use tinyvec::ArrayVec;

// abstraction types
pub mod bezier_segment;
// specialized types
pub mod cubic_bezier;
pub mod line;
pub mod quadratic_bezier;
// generic types
pub mod bezier;
pub mod bspline;
pub mod point_generic;

// Traits
pub mod point;
pub mod path;
pub mod bspline_path;
pub mod find_root;

mod roots;

#[cfg(feature = "nalgebra")]
mod adapters;

// export common types at crate root
pub use bezier::Bezier;
pub use bspline::BSpline;
pub use cubic_bezier::CubicBezier;
pub use line::LineSegment;
pub use point::Point;
pub use point::{PointDot, PointIndex, PointNorm};
pub use point_generic::PointN;
pub use quadratic_bezier::QuadraticBezier;
pub use bspline_path::BSplinePath;
pub use find_root::FindRoot;
pub use path::BezierPath;

// Per-architecture epsilon used by tests.
#[cfg(test)]
#[cfg(target_pointer_width = "64")]
pub(crate) const EPSILON: f64 = f64::EPSILON;
#[cfg(test)]
#[cfg(not(target_pointer_width = "64"))]
pub(crate) const EPSILON: f32 = f32::EPSILON;
