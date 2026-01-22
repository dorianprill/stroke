//! Mixed Bezier path utilities.

use core::slice;

use num_traits::{Float, NumCast};
use crate::bezier_segment::BezierSegment;
use super::{ArrayVec, CubicBezier, LineSegment, Point, PointDot, PointIndex, PointNorm, QuadraticBezier};

const DEFAULT_LENGTH_STEPS: usize = 64;
const MAX_LENGTH_ITERS: usize = 32;

/// A path composed of mixed Bezier segments (line/quadratic/cubic).
///
/// Component-aware helpers (e.g. bounding boxes) require `PointIndex`.
///
/// # Examples
/// ```rust
/// use stroke::{BezierPath, LineSegment, PointN};
///
/// let mut path: BezierPath<PointN<f64, 2>, 2> = BezierPath::new();
/// path.push_line(LineSegment::new(
///     PointN::new([0.0, 0.0]),
///     PointN::new([1.0, 0.0]),
/// ));
/// path.push_line(LineSegment::new(
///     PointN::new([1.0, 0.0]),
///     PointN::new([2.0, 0.0]),
/// ));
///
/// assert_eq!(path.start().unwrap(), PointN::new([0.0, 0.0]));
/// assert_eq!(path.end().unwrap(), PointN::new([2.0, 0.0]));
///
/// let t = path.t_at_length(1.0).unwrap();
/// let p = path.point_at_length(1.0).unwrap();
/// let tangent = path.tangent(0.25).unwrap();
/// # let _ = (t, p, tangent);
/// ```
pub struct BezierPath<P, const N: usize>
where
    P: Point,
{
    segments: ArrayVec<[BezierSegment<P>; N]>,
}

impl<P, const N: usize> BezierPath<P, N>
where
    P: Point,
    [BezierSegment<P>; N]: tinyvec::Array<Item = BezierSegment<P>>,
{
    /// Create an empty path with capacity `N`.
    pub fn new() -> Self {
        BezierPath {
            segments: ArrayVec::new(),
        }
    }

    /// Return the number of segments currently stored.
    pub fn len(&self) -> usize {
        self.segments.len()
    }

    /// Return true if the path is empty.
    pub fn is_empty(&self) -> bool {
        self.segments.len() == 0
    }

    /// Iterate over the stored segments.
    pub fn segments(&self) -> slice::Iter<'_, BezierSegment<P>> {
        self.segments.iter()
    }

    /// Push a segment, returning false if the path is at capacity.
    pub fn push(&mut self, segment: BezierSegment<P>) -> bool {
        if self.segments.len() < self.segments.capacity() {
            self.segments.push(segment);
            true
        } else {
            false
        }
    }

    /// Push a line segment, returning false if the path is at capacity.
    pub fn push_line(&mut self, segment: LineSegment<P>) -> bool {
        self.push(segment.into())
    }

    /// Push a quadratic segment, returning false if the path is at capacity.
    pub fn push_quadratic(&mut self, segment: QuadraticBezier<P>) -> bool {
        self.push(segment.into())
    }

    /// Push a cubic segment, returning false if the path is at capacity.
    pub fn push_cubic(&mut self, segment: CubicBezier<P>) -> bool {
        self.push(segment.into())
    }

    /// Evaluate a point along the path for `t` in `[0, 1]`.
    /// Values outside the range are clamped. Returns None for empty paths.
    pub fn eval(&self, t: P::Scalar) -> Option<P> {
        let (index, local_t) = self.segment_parameter(t)?;
        Some(self.segments[index].eval(local_t))
    }

    /// Return the start point of the path.
    pub fn start(&self) -> Option<P> {
        let zero = <P::Scalar as NumCast>::from(0.0).unwrap();
        self.eval(zero)
    }

    /// Return the end point of the path.
    pub fn end(&self) -> Option<P> {
        let one = <P::Scalar as NumCast>::from(1.0).unwrap();
        self.eval(one)
    }

    /// Return a path with reversed direction.
    pub fn reverse(&self) -> Self {
        let mut reversed = BezierPath::new();
        for segment in self.segments.iter().rev() {
            let _ = reversed.push(segment.reverse());
        }
        reversed
    }

    /// Return the unit tangent direction at `t`.
    pub fn tangent(&self, t: P::Scalar) -> Option<P>
    where
        P: PointNorm,
    {
        let (index, local_t) = self.segment_parameter(t)?;
        Some(self.segments[index].tangent(local_t))
    }

    /// Return the curvature magnitude at `t`.
    pub fn curvature(&self, t: P::Scalar) -> Option<P::Scalar>
    where
        P: PointNorm + PointDot,
    {
        let (index, local_t) = self.segment_parameter(t)?;
        Some(self.segments[index].curvature(local_t))
    }

    /// Return the principal normal direction at `t`.
    ///
    /// Returns `None` if the path is empty or curvature is undefined.
    pub fn normal(&self, t: P::Scalar) -> Option<P>
    where
        P: PointNorm + PointDot,
    {
        let (index, local_t) = self.segment_parameter(t)?;
        self.segments[index].normal(local_t)
    }

    /// Approximate parameter `t` at arc length `s`.
    pub fn t_at_length_approx(&self, s: P::Scalar, nsteps: usize) -> Option<P::Scalar>
    where
        P: PointNorm,
    {
        let zero = <P::Scalar as NumCast>::from(0.0).unwrap();
        let one = <P::Scalar as NumCast>::from(1.0).unwrap();
        let half = <P::Scalar as NumCast>::from(0.5).unwrap();

        let total = self.arclen_partial(one, nsteps)?;
        if total <= P::Scalar::epsilon() {
            return Some(zero);
        }
        let target = s.clamp(zero, total);

        let mut lo = zero;
        let mut hi = one;
        for _ in 0..MAX_LENGTH_ITERS {
            let mid = (lo + hi) * half;
            let len = self.arclen_partial(mid, nsteps)?;
            if len < target {
                lo = mid;
            } else {
                hi = mid;
            }
        }

        Some((lo + hi) * half)
    }

    /// Approximate parameter `t` at arc length `s` using a default resolution.
    pub fn t_at_length(&self, s: P::Scalar) -> Option<P::Scalar>
    where
        P: PointNorm,
    {
        self.t_at_length_approx(s, DEFAULT_LENGTH_STEPS)
    }

    /// Evaluate the point at arc length `s`.
    pub fn point_at_length_approx(&self, s: P::Scalar, nsteps: usize) -> Option<P>
    where
        P: PointNorm,
    {
        let t = self.t_at_length_approx(s, nsteps)?;
        self.eval(t)
    }

    /// Evaluate the point at arc length `s` using a default resolution.
    pub fn point_at_length(&self, s: P::Scalar) -> Option<P>
    where
        P: PointNorm,
    {
        self.point_at_length_approx(s, DEFAULT_LENGTH_STEPS)
    }

    /// Return the bounding box across all segments. Returns None for empty paths.
    pub fn bounding_box(&self) -> Option<[(P::Scalar, P::Scalar); P::DIM]>
    where
        P: PointIndex,
        [P::Scalar; 1]: tinyvec::Array<Item = P::Scalar>,
        [P::Scalar; 2]: tinyvec::Array<Item = P::Scalar>,
        [P::Scalar; 3]: tinyvec::Array<Item = P::Scalar>,
        [P::Scalar; 4]: tinyvec::Array<Item = P::Scalar>,
    {
        let mut iter = self.segments.iter();
        let first = match iter.next() {
            Some(segment) => segment.bounding_box(),
            None => return None,
        };

        let mut bounds = first;
        for segment in iter {
            let segment_bounds = segment.bounding_box();
            for dim in 0..P::DIM {
                if segment_bounds[dim].0 < bounds[dim].0 {
                    bounds[dim].0 = segment_bounds[dim].0;
                }
                if segment_bounds[dim].1 > bounds[dim].1 {
                    bounds[dim].1 = segment_bounds[dim].1;
                }
            }
        }

        Some(bounds)
    }

    fn segment_parameter(&self, t: P::Scalar) -> Option<(usize, P::Scalar)> {
        let count = self.segments.len();
        if count == 0 {
            return None;
        }

        let zero = <P::Scalar as NumCast>::from(0.0).unwrap();
        let one = <P::Scalar as NumCast>::from(1.0).unwrap();
        let t = t.clamp(zero, one);

        let count_scalar = <P::Scalar as NumCast>::from(count as f64).unwrap();
        let scaled = t * count_scalar;
        if scaled >= count_scalar {
            return Some((count - 1, one));
        }

        let index = num_traits::cast::<P::Scalar, usize>(scaled.floor()).unwrap_or(0);
        let index_scalar = <P::Scalar as NumCast>::from(index as f64).unwrap();
        let local = scaled - index_scalar;
        Some((index, local))
    }

    fn arclen_partial(&self, t: P::Scalar, nsteps: usize) -> Option<P::Scalar>
    where
        P: PointNorm,
    {
        let zero = <P::Scalar as NumCast>::from(0.0).unwrap();
        let one = <P::Scalar as NumCast>::from(1.0).unwrap();
        let t = t.clamp(zero, one);
        if t <= zero {
            return Some(zero);
        }

        let nsteps = nsteps.max(1);
        let nsteps_scalar = <P::Scalar as NumCast>::from(nsteps as f64).unwrap();
        let dt = t / nsteps_scalar;
        let mut arclen = zero;
        let mut prev = self.eval(zero)?;

        for i in 1..=nsteps {
            let i_scalar = <P::Scalar as NumCast>::from(i as f64).unwrap();
            let ti = if i == nsteps { t } else { dt * i_scalar };
            let p = self.eval(ti)?;
            arclen = arclen + (p - prev).squared_norm().sqrt();
            prev = p;
        }

        Some(arclen)
    }
}

impl<P, const N: usize> Default for BezierPath<P, N>
where
    P: Point,
    [BezierSegment<P>; N]: tinyvec::Array<Item = BezierSegment<P>>,
{
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use crate::{PointN, PointNorm, EPSILON};
    use super::*;

    #[test]
    fn bezier_path_eval_segments() {
        let mut path: BezierPath<PointN<f64, 2>, 4> = BezierPath::new();
        path.push_line(LineSegment::new(
            PointN::new([0.0, 0.0]),
            PointN::new([1.0, 0.0]),
        ));
        path.push_line(LineSegment::new(
            PointN::new([1.0, 0.0]),
            PointN::new([1.0, 1.0]),
        ));

        let p0 = path.eval(0.0).unwrap();
        assert_eq!(p0, PointN::new([0.0, 0.0]));

        let p1 = path.eval(0.25).unwrap();
        assert!((p1 - PointN::new([0.5, 0.0])).squared_norm() < EPSILON);

        let p2 = path.eval(0.75).unwrap();
        assert!((p2 - PointN::new([1.0, 0.5])).squared_norm() < EPSILON);
    }

    #[test]
    fn bezier_path_bounds_union() {
        let mut path: BezierPath<PointN<f64, 2>, 4> = BezierPath::new();
        path.push_line(LineSegment::new(
            PointN::new([0.0, 0.0]),
            PointN::new([1.0, 0.0]),
        ));
        path.push_line(LineSegment::new(
            PointN::new([1.0, 0.0]),
            PointN::new([1.0, 1.0]),
        ));

        let bounds = path.bounding_box().unwrap();
        assert_eq!(bounds[0], (0.0, 1.0));
        assert_eq!(bounds[1], (0.0, 1.0));
    }

    #[test]
    fn bezier_path_clamps_out_of_range() {
        let mut path: BezierPath<PointN<f64, 2>, 2> = BezierPath::new();
        path.push_line(LineSegment::new(
            PointN::new([0.0, 0.0]),
            PointN::new([1.0, 0.0]),
        ));
        path.push_line(LineSegment::new(
            PointN::new([1.0, 0.0]),
            PointN::new([1.0, 1.0]),
        ));

        let start = path.eval(-1.0).unwrap();
        assert_eq!(start, PointN::new([0.0, 0.0]));

        let end = path.eval(2.0).unwrap();
        assert_eq!(end, PointN::new([1.0, 1.0]));
    }

    #[test]
    fn bezier_path_capacity() {
        let mut path: BezierPath<PointN<f64, 2>, 1> = BezierPath::new();
        let first = path.push_line(LineSegment::new(
            PointN::new([0.0, 0.0]),
            PointN::new([1.0, 0.0]),
        ));
        let second = path.push_line(LineSegment::new(
            PointN::new([1.0, 0.0]),
            PointN::new([2.0, 0.0]),
        ));

        assert!(first);
        assert!(!second);
        assert_eq!(path.len(), 1);
    }

    #[test]
    fn bezier_path_api_parity() {
        let mut path: BezierPath<PointN<f64, 2>, 2> = BezierPath::new();
        path.push_line(LineSegment::new(
            PointN::new([0.0, 0.0]),
            PointN::new([1.0, 0.0]),
        ));
        path.push_line(LineSegment::new(
            PointN::new([1.0, 0.0]),
            PointN::new([2.0, 0.0]),
        ));

        assert_eq!(path.start().unwrap(), PointN::new([0.0, 0.0]));
        assert_eq!(path.end().unwrap(), PointN::new([2.0, 0.0]));

        let reversed = path.reverse();
        let p0 = path.eval(0.25).unwrap();
        let p1 = reversed.eval(0.75).unwrap();
        assert!((p0 - p1).squared_norm() < EPSILON);

        let tangent = path.tangent(0.25).unwrap();
        assert!((tangent[0] - 1.0).abs() < EPSILON);
        assert!(tangent[1].abs() < EPSILON);

        let curvature = path.curvature(0.25).unwrap();
        assert!(curvature.abs() < EPSILON);
        assert!(path.normal(0.25).is_none());

        let t = path.t_at_length(1.0).unwrap();
        assert!((t - 0.5).abs() < 1e-6);

        let p = path.point_at_length(1.0).unwrap();
        assert!((p - PointN::new([1.0, 0.0])).squared_norm() < EPSILON);
    }

}
