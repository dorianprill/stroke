#![no_std]
#![allow(incomplete_features)]
#![feature(const_generics)]

use core::ops::{Add, Sub, Mul};

extern crate num_traits;
use num_traits::{float::Float};

extern crate tinyvec;
use tinyvec::ArrayVec;

pub mod point;
pub mod point2;
pub mod point_generic;
pub mod line;
pub mod bezier;
pub mod quadratic_bezier;
pub mod cubic_bezier;
pub mod bezier_segment;
pub mod rational_bezier;
pub mod bspline;


// conditionally compiled newtype pattern used to determine which size float to use and for tests
// so that the library can abstract over both 32 and 64 bit architectures
// TODO there has to be a better architectural solution to this... if not, generate impls with macro
#[cfg(target_pointer_width = "64")]
type NativeFloat = f64;
// same thing for 32 bit
#[cfg(target_pointer_width = "32")]
type NativeFloat = f32;