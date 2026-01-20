//! Point traits and capability helpers.

use core::ops::{Add, Index, Mul, Sub};
use num_traits::Float;

/// Core point abstraction used across the library.
///
/// `Point` models a fixed-dimension vector space element. Curve evaluation
/// only relies on these operations; geometry helpers are provided via
/// capability traits like [`PointIndex`] and [`PointNorm`].
pub trait Point:
    Add<Self, Output = Self> + Sub<Self, Output = Self> + Mul<Self::Scalar, Output = Self> + Copy
{
    /// Scalar type used for interpolation and arithmetic.
    type Scalar: Float;
    /// Compile-time dimension of the point.
    const DIM: usize;
}

/// Optional capability trait for component access via indexing.
///
/// This is required for methods that need per-axis values (e.g. bounding boxes).
pub trait PointIndex: Point + Index<usize, Output = Self::Scalar> {}

impl<T> PointIndex for T where T: Point + Index<usize, Output = Self::Scalar> {}

/// Optional capability trait for dot products.
///
/// Used by distance and projection helpers.
pub trait PointDot: PointIndex {
    fn dot(&self, other: &Self) -> Self::Scalar {
        let mut sum = <Self::Scalar as num_traits::NumCast>::from(0.0).unwrap();
        for i in 0..Self::DIM {
            sum = sum + self[i] * other[i];
        }
        sum
    }
}

impl<T> PointDot for T where T: PointIndex {}

/// Optional capability trait for norms.
///
/// Provides `squared_norm()`; it does not require allocation or dynamic sizing.
pub trait PointNorm: PointDot {
    fn squared_norm(&self) -> Self::Scalar {
        self.dot(self)
    }
}

impl<T> PointNorm for T where T: PointDot {}
