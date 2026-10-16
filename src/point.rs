//! Point traits and capability helpers.

use core::ops::{Add, Index, Mul, Sub};
use num_traits::{Float, Zero};

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

    /// Return whether every coordinate is finite.
    ///
    /// Numeric representations without NaN or infinity should return `true`.
    fn is_finite(&self) -> bool;
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
/// Provides stable norm calculations without allocation or dynamic sizing.
pub trait PointNorm: PointDot {
    fn squared_norm(&self) -> Self::Scalar {
        self.dot(self)
    }

    /// Return the Euclidean norm using scaled accumulation to avoid avoidable
    /// intermediate overflow and underflow.
    fn norm(&self) -> Self::Scalar {
        let mut scale = Self::Scalar::zero();
        for axis in 0..Self::DIM {
            if !self[axis].is_finite() {
                return self[axis].abs();
            }
            scale = scale.max(self[axis].abs());
        }
        if scale == Self::Scalar::zero() {
            return scale;
        }

        let mut sum = Self::Scalar::zero();
        for axis in 0..Self::DIM {
            let normalized = self[axis] / scale;
            sum = sum + normalized * normalized;
        }
        scale * sum.sqrt()
    }
}

impl<T> PointNorm for T where T: PointDot {}
