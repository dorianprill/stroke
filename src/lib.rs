#![no_std]
#![forbid(unsafe_code)]
#![allow(incomplete_features)]
// Curve derivatives use array lengths such as N-1.
#![feature(generic_const_exprs)]

//! A const-generic, no-std spline library for Bezier and B-spline curves.
//!
//! Curves use the [`Point`] trait for addition, subtraction and scalar
//! multiplication in a fixed number of dimensions. Use [`ArrayPoint`] for
//! array-backed coordinates, or implement [`Point`] for your own type.
//! Coordinate indexing through [`PointIndex`] adds bounds queries.
//! [`PointDot`] and [`PointNorm`] add dot products and Euclidean lengths.
//!
//! `Point::Scalar` must implement [`num_traits::Float`], as `f32` and `f64` do.
//! Nightly Rust is required. For methods with const-expression return types,
//! such as bounding boxes, also enable `#![feature(generic_const_exprs)]`
//! in your crate root.
//!
//! Arc length uses adaptive Gauss-Legendre integration with a requested error
//! tolerance. Invalid inputs, non-finite arithmetic and exhausted subdivision
//! budgets return [`ArcLengthError`].
//!
//! `bounding_box()` bounds the curve along the coordinate axes. `control_box()`
//! bounds its control points and may be larger. Generic Beziers and B-splines
//! also have approximate oriented boxes through
//! [`Bezier::bounding_box_oriented_approx`] and [`BSpline::bounding_box_oriented_approx`].
//! Their rotation search has no guarantee of finding a minimum-volume box.
//!
//! # Examples
//!
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
//!
//! - `nalgebra`: implements `Point` for `nalgebra::SVector<T, D>`. Add `nalgebra`
//!   as a dependency in your crate to construct its vectors.

use tinyvec::ArrayVec;

pub mod bezier_segment;
pub mod bounding_box;
pub mod oriented_bounding_box;

pub mod cubic_bezier;
pub mod error;
pub mod line;
pub mod quadratic_bezier;

pub mod array_point;
pub mod bezier;
pub mod bspline;

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
pub use oriented_bounding_box::{
    OrientedBoundingBox, OrientedBoundingBoxDiagnostics, OrientedBoundingBoxError,
};
pub use path::BezierPath;
pub use point::Point;
pub use point::{PointDot, PointIndex, PointNorm};
pub use quadratic_bezier::QuadraticBezier;
pub use roots::RootFindingError;
