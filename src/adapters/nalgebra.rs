//! Use nalgebra vectors as curve points.
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
//! ]).unwrap();
//!
//! let mid = curve.eval(0.5);
//! # let _ = mid;
//! ```
//!
//! The scalar type must implement `nalgebra::RealField` and `num_traits::Float`,
//! as `f32` and `f64` do.

use nalgebra::{RealField, SVector};
use num_traits::Float;

use crate::point::Point;

impl<T, const D: usize> Point for SVector<T, D>
where
    T: RealField + Float,
{
    type Scalar = T;
    const DIM: usize = D;

    fn is_finite(&self) -> bool {
        self.iter().all(|value| Float::is_finite(*value))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn vector_finiteness_checks_every_coordinate() {
        assert!(Point::is_finite(&SVector::<f64, 2>::new(0.0, -0.0)));
        assert!(!Point::is_finite(&SVector::<f64, 2>::new(0.0, f64::NAN)));
    }
}
