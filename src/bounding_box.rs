//! Axis-aligned bounding boxes.

use core::ops::Index;

use crate::PointIndex;
use num_traits::Float;

/// An axis-aligned bounding box in `D` dimensions.
///
/// `bounding_box()` on a curve bounds its coordinate extrema. `control_box()`
/// uses its controls instead. Both methods return this type.
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

    /// Return the minimum coordinates.
    pub fn min(&self) -> [T; D] {
        self.axes.map(|(min, _)| min)
    }

    /// Return the maximum coordinates.
    pub fn max(&self) -> [T; D] {
        self.axes.map(|(_, max)| max)
    }

    /// Return whether the boxes overlap, including touching boundaries.
    pub fn intersects(&self, other: &Self) -> bool
    where
        T: PartialOrd,
    {
        (0..D).all(|i| self.axes[i].0 <= other.axes[i].1 && other.axes[i].0 <= self.axes[i].1)
    }

    /// Return the overlapping region, or `None` if the boxes are disjoint.
    /// Touching boundaries give a zero-width intersection.
    pub fn intersection(&self, other: &Self) -> Option<Self>
    where
        T: PartialOrd,
    {
        self.intersects(other).then(|| Self {
            axes: core::array::from_fn(|i| {
                let (a, b) = self.axes[i];
                let (c, d) = other.axes[i];
                (if a > c { a } else { c }, if b < d { b } else { d })
            }),
        })
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
    /// Return the center, avoiding intermediate overflow for finite bounds.
    pub fn center(&self) -> [T; D] {
        let half = T::one() / (T::one() + T::one());
        self.axes.map(|(min, max)| {
            if min.abs() <= T::max_value() * half && max.abs() <= T::max_value() * half {
                return (min + max) * half;
            }
            let width = max - min;
            if width.is_finite() {
                min + width * half
            } else {
                min * half + max * half
            }
        })
    }

    /// Return axis widths. A width can be infinite even when both bounds are finite.
    pub fn extent(&self) -> [T; D] {
        self.axes.map(|(min, max)| max - min)
    }

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
