#![no_std]
#![forbid(unsafe_code)]
#![allow(incomplete_features)]
// let's get daring...
#![feature(const_generics)]
// this is needed to use expressions in const generics such as N-1 (see curve derivatives)
#![feature(const_evaluatable_checked)]
// this feature is was necessary to build an iterator over a
// const generic array (see PointN impl) before 1.51
//#![feature(array_value_iter)]
//  use of unstable library feature 'array_value_iter'
//  see issue #65798 <https://github.com/rust-lang/rust/issues/65798> for more information

// removes the need for generics with associated types to specify the
// associated type like P:Point instead of P: Point<Scalar=f64>
#![feature(associated_type_bounds)]

use core::ops::{Add, Mul, Sub};

extern crate num_traits;
use num_traits::float::Float;

extern crate tinyvec;
use tinyvec::ArrayVec;

pub mod bezier;
pub mod bezier_segment;
pub mod cubic_bezier;
pub mod line;
pub mod point;
pub mod point_generic;
pub mod quadratic_bezier;
//pub mod rational_bezier;
pub mod bspline;

// export common types at crate root
pub use bezier::Bezier;
pub use bspline::BSpline;
pub use cubic_bezier::CubicBezier;
pub use line::LineSegment;
pub use point::Point;
pub use point_generic::PointN;
pub use quadratic_bezier::QuadraticBezier;

// Conditionally compiled newtype pattern used to determine which size float to use
// so that the library can abstract over both 32 and 64 bit architectures (maybe even 16 bit)
// TODO there has to be a better architectural solution to this...
//   Option 1.: Is is possible to determine float size at build time?
//   Option 2.: Is it possible to make this value a 'must-config' in cargo.toml?
#[cfg(target_pointer_width = "64")]
type NativeFloat = f64;
#[cfg(target_pointer_width = "64")]
const EPSILON: f64 = f64::EPSILON;
// same thing for 32 bit
#[cfg(target_pointer_width = "32")]
type NativeFloat = f32;
#[cfg(target_pointer_width = "32")]
const EPSILON: f32 = f32::EPSILON;
