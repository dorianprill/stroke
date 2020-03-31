#![no_std]
use core::ops::{Add, Sub, Mul};

extern crate num_traits;
use num_traits::{float::Float};

extern crate arrayvec;
use arrayvec::ArrayVec;

pub mod point2;
pub mod line;
pub mod quadratic_bezier;
pub mod cubic_bezier;
pub mod bezier_segment;


// conditionally compiled newtype pattern used to determine which size float to use in arclen() and for tests
// so that the library can abstract over both 32 and 64 bit architectures
// TODO there has to be a better architectural solution to this... if not, generate impls with macro
#[cfg(target_pointer_width = "64")]
type NativeFloat = f64;
#[cfg(target_pointer_width = "64")]
type NativeInt   = i64;
#[cfg(target_pointer_width = "64")]
type NativeUInt  = u64;

// same thing for 32 bit
#[cfg(target_pointer_width = "32")]
type NativeFloat = f32;
#[cfg(target_pointer_width = "32")]
type NativeInt   = i32;
#[cfg(target_pointer_width = "32")]
type NativeUInt  = u32;