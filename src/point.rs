//! Point arithmetic, coordinate access and Euclidean lengths.

use core::ops::{Add, Index, Mul, Sub};
use num_traits::{Float, Zero};

/// A fixed-dimensional point with addition, subtraction and scalar multiplication.
///
/// These operations are enough to evaluate curves. Coordinate indexing through
/// [`PointIndex`] adds bounds queries. [`PointNorm`] adds Euclidean geometry.
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

/// Coordinate access through `point[axis]`, used by bounding-box queries.
///
/// Implemented automatically for any [`Point`] with the required indexing.
pub trait PointIndex: Point + Index<usize, Output = Self::Scalar> {}

impl<T> PointIndex for T where T: Point + Index<usize, Output = Self::Scalar> {}

/// Dot products for distance and projection calculations.
///
/// Implemented automatically for every [`PointIndex`].
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

/// Euclidean lengths for points and difference vectors.
///
/// Implemented automatically for every [`PointDot`].
pub trait PointNorm: PointDot {
    fn squared_norm(&self) -> Self::Scalar {
        self.dot(self)
    }

    /// Return the Euclidean length.
    /// Scale coordinates before squaring to limit intermediate overflow and underflow.
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
