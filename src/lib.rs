#![cfg_attr(not(test), no_std)]
#![forbid(unsafe_code)]
#![allow(incomplete_features)]
// this is needed to use expressions in const generics such as N-1 (see curve derivatives)
#![feature(generic_const_exprs)]
// this feature is needed for tinyvec < 2.0 to compile for const generic arrays like ArrayVec<[f32;N]>
//#![feature(min_const_generics)]
// removes the need for generics with associated types to specify the
// associated type like P:Point instead of P: Point<Scalar=f64>
#![feature(associated_type_bounds)]
// make splines usable as Fn in trait bounds
//#![feature(fn_traits)]

extern crate num_traits;
use num_traits::float::Float;

extern crate tinyvec;
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
pub use nalgebra;

// Traits
pub mod spline;

mod roots;

// export common types at crate root
pub use bezier::Bezier;
pub use bspline::BSpline;
pub use cubic_bezier::CubicBezier;
pub use line::LineSegment;

pub use quadratic_bezier::QuadraticBezier;
pub use spline::Spline;

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
const EPSILON: f64 = 0.000001;
// For now, we fix all non-64 bit architectures to 32bit floats
// as smaller-width architectures are more likely to have different int/float sizes if they have a fpu
#[cfg(not(target_pointer_width = "64"))]
type NativeFloat = f32;
#[cfg(not(target_pointer_width = "64"))]
const EPSILON: f32 = f32::EPSILON;
// This might change when/if Rust gets additional types like f16 or f24
