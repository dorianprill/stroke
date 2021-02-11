#![no_std]
#![forbid(unsafe_code)]
#![allow(incomplete_features)]
// let's get daring...
#![feature(const_generics)]
// this is needed to use expressions in const generics such as N-1 (see curve derivatives)
#![feature(const_evaluatable_checked)]
// this feature is currently necessary to build an iterator over a const generic array (see PointN impl)
#![feature(array_value_iter)]
//use of unstable library feature 'array_value_iter'
//see issue #65798 <https://github.com/rust-lang/rust/issues/65798> for more information

use core::ops::{Add, Sub, Mul};

extern crate num_traits;
use num_traits::{float::Float};

extern crate tinyvec;
use tinyvec::ArrayVec;

pub mod point;
pub mod point_generic;
pub mod line;
pub mod bezier;
pub mod quadratic_bezier;
pub mod cubic_bezier;
pub mod bezier_segment;
//pub mod rational_bezier;
pub mod bspline;


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