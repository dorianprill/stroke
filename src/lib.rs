#![no_std]
#![forbid(unsafe_code)]
#![allow(incomplete_features)]
// this is needed to use expressions in const generics such as N-1 (see curve derivatives)
#![feature(generic_const_exprs)]

//! A const-generic, no-std spline library for Bezier and B-spline curves.
//!
//! The core abstraction is the [`Point`] trait, which models
//! a fixed-dimension vector space element. Optional capability traits such as
//! [`PointIndex`], [`PointDot`], and [`PointNorm`] enable component access and
//! geometric helpers when required.
//!
//! `Point::Scalar` must implement [`num_traits::Float`], so the library currently
//! targets floating-point numeric representations such as `f32` and `f64`.
//! Consumers using const-expression-based APIs, such as bounding boxes, should
//! enable `#![feature(generic_const_exprs)]` at their crate root on nightly.
//!
//! Curves and paths expose accuracy-controlled Gauss-Legendre arc length and
//! tight [`BoundingBox`] values. Their `control_box()` methods return the
//! conservative box of the control points. Arc-length operations return
//! [`ArcLengthError`] for invalid input, non-finite arithmetic or failure to converge.
//!
//! # Examples
//! ```rust
//! use stroke::{Bezier, ArrayPoint};
//!
//! let curve = Bezier::<ArrayPoint<f32, 2>, 3>::new([
//!     ArrayPoint::new([0.0, 0.0]),
//!     ArrayPoint::new([1.0, 0.0]),
//!     ArrayPoint::new([1.0, 1.0]),
//! ]).unwrap();
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
pub mod bounding_box;
// specialized types
pub mod cubic_bezier;
pub mod error;
pub mod line;
pub mod quadratic_bezier;
// generic types
pub mod array_point;
pub mod bezier;
pub mod bspline;

// Traits
pub mod bspline_path;
pub mod find_root;
pub mod path;
pub mod point;

mod arc_length;
mod bernstein;
mod conversions;
mod differential;
mod distance;
mod numeric;
mod path_parameter;
mod roots;

#[cfg(test)]
mod test_support;

#[cfg(test)]
mod geometry_tests;

#[cfg(test)]
mod geometry_regressions;

#[cfg(feature = "nalgebra")]
mod adapters;

// export common types at crate root
pub use arc_length::ArcLengthError;
pub use array_point::ArrayPoint;
pub use bezier::Bezier;
pub use bezier_segment::BezierSegment;
pub use bounding_box::BoundingBox;
pub use bspline::{BSpline, BSplineError};
pub use bspline_path::BSplinePath;
pub use cubic_bezier::CubicBezier;
pub use error::{BezierBasisError, CurveConstructionError, CurveParameterError};
pub use find_root::FindRoot;
pub use line::LineSegment;
pub use path::BezierPath;
pub use point::Point;
pub use point::{PointDot, PointIndex, PointNorm};
pub use quadratic_bezier::QuadraticBezier;
pub use roots::RootFindingError;
