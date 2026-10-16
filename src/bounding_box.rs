//! Axis-aligned bounding boxes.

use core::ops::Index;

use crate::PointIndex;
use num_traits::Float;

/// An axis-aligned bounding box in `D` dimensions.
///
/// Tightness relative to a curve depends on how the box was constructed.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BoundingBox<T, const D: usize> {
    axes: [(T, T); D],
}

impl<T: Copy, const D: usize> BoundingBox<T, D> {
    /// Construct a bounding box from ordered `(minimum, maximum)` axis bounds.
    pub(crate) fn from_axes_unchecked(axes: [(T, T); D]) -> Self {
        Self { axes }
    }

    /// Return the `(minimum, maximum)` bounds for every axis.
    pub fn axes(&self) -> &[(T, T); D] {
        &self.axes
    }

    /// Return the minimum and maximum coordinate on an axis.
    /// Panics if `axis >= D`.
    pub fn axis(&self, axis: usize) -> (T, T) {
        self.axes[axis]
    }

    /// Return whether the point is inside this box, including its boundary.
    pub fn contains<P>(&self, point: &P) -> bool
    where
        P: PointIndex<Scalar = T>,
        T: PartialOrd,
    {
        P::DIM == D
            && (0..D)
                .all(|axis| point[axis] >= self.axes[axis].0 && point[axis] <= self.axes[axis].1)
    }

    pub(crate) fn from_points<P>(points: &[P]) -> Self
    where
        P: PointIndex<Scalar = T>,
        T: Float,
    {
        let axes = core::array::from_fn(|axis| {
            let mut min = points[0][axis];
            let mut max = min;
            for point in &points[1..] {
                min = min.min(point[axis]);
                max = max.max(point[axis]);
            }
            (min, max)
        });
        Self { axes }
    }

    /// Convert this bounding box into its axis bounds.
    pub fn into_axes(self) -> [(T, T); D] {
        self.axes
    }

    /// Expand this box to include `other`.
    pub fn union(&self, other: &Self) -> Self
    where
        T: PartialOrd,
    {
        let axes = core::array::from_fn(|axis| {
            let min = if self.axes[axis].0 < other.axes[axis].0 {
                self.axes[axis].0
            } else {
                other.axes[axis].0
            };
            let max = if self.axes[axis].1 > other.axes[axis].1 {
                self.axes[axis].1
            } else {
                other.axes[axis].1
            };
            (min, max)
        });
        Self { axes }
    }
}

impl<T: Float, const D: usize> BoundingBox<T, D> {
    /// Construct a finite bounding box from minimum and maximum coordinates.
    /// Returns `None` if a coordinate is non-finite or a minimum exceeds its maximum.
    pub fn new(min: [T; D], max: [T; D]) -> Option<Self> {
        let mut valid = true;
        let axes = core::array::from_fn(|axis| {
            let lo = min[axis];
            let hi = max[axis];
            valid &= lo.is_finite() && hi.is_finite() && lo <= hi;
            (lo, hi)
        });
        valid.then_some(Self { axes })
    }
}

impl<T, const D: usize> Index<usize> for BoundingBox<T, D> {
    type Output = (T, T);

    fn index(&self, index: usize) -> &Self::Output {
        &self.axes[index]
    }
}
