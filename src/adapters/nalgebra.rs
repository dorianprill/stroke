//! Nalgebra adapter implementations.
//!
//! Enable this adapter with the `nalgebra` feature to use `nalgebra::SVector<T, D>`
//! as a `Point`. Add `nalgebra` as a direct dependency to construct the vectors
//! in your own code.
//!
//! # Example
//! ```rust,no_run
//! use nalgebra::SVector;
//! use stroke::Bezier;
//!
//! let curve = Bezier::<SVector<f32, 2>, 3>::new([
//!     SVector::<f32, 2>::new(0.0, 0.0),
//!     SVector::<f32, 2>::new(1.0, 0.0),
//!     SVector::<f32, 2>::new(1.0, 1.0),
//! ]);
//!
//! let mid = curve.eval(0.5);
//! # let _ = mid;
//! ```
//!
//! The scalar type must satisfy `nalgebra::RealField` and `num_traits::Float`
//! (e.g. `f32` or `f64`).

use nalgebra::{RealField, SVector};
use num_traits::Float;

use crate::point::Point;

impl<T, const D: usize> Point for SVector<T, D>
where
    T: RealField + Float,
{
    type Scalar = T;
    const DIM: usize = D;
}
