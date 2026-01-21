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

// this feature was needed for tinyvec < 2.0 to compile for const generic arrays like ArrayVec<[f32;N]>
//#![feature(min_const_generics)]

// NO LONGER NECESSARY (stabilized in 1.79)
// removes the need for generics with associated types to specify the
// associated type like P:Point instead of P: Point<Scalar=f64>
//#![feature(associated_type_bounds)]

// make splines usable as Fn in trait bounds
//#![feature(fn_traits)]

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

// Conditionally compiled newtype pattern used to determine which size float to use for internal constants
// so that the library can specialize internal types for the architecture for best performance
// TODO An FPU size definition is not available to the compiler, so

//   1. Either fix everything to 32-bit and accept performance loss (how it's done now)

//   2. Make this value a 'must-config' in cargo.toml (requires either external libraries
//      or must wait for additional float types)

// If it is a modern 64 bit architecture, it likely also has a 64 bit FPU
#[cfg(target_pointer_width = "64")]
type NativeFloat = f64;
#[cfg(target_pointer_width = "64")]
const EPSILON: f64 = f64::EPSILON;
// For now, we fix all non-64 bit architectures to 32bit floats
// as smaller-width architectures are more likely to have different int/float sizes if they have a fpu
#[cfg(not(target_pointer_width = "64"))]
type NativeFloat = f32;
#[cfg(not(target_pointer_width = "64"))]
const EPSILON: f32 = f32::EPSILON;
// This might change when/if Rust gets additional types like f16 or f24
