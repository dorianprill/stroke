//! Mixed Bezier path utilities.

use crate::arc_length::validate_accuracy;
use crate::{ArcLengthError, BoundingBox};
use core::slice;

use super::{
    ArrayVec, CubicBezier, LineSegment, Point, PointDot, PointIndex, PointNorm, QuadraticBezier,
};
use crate::bezier_segment::BezierSegment;
use num_traits::{Float, NumCast};

/// A path composed of mixed Bezier segments (line/quadratic/cubic).
///
/// Component-aware helpers (e.g. bounding boxes) require `PointIndex`.
///
/// # Examples
/// ```rust
/// use stroke::{BezierPath, LineSegment, ArrayPoint};
///
/// let mut path: BezierPath<ArrayPoint<f64, 2>, 2> = BezierPath::new();
/// path.push_line(LineSegment::new(
///     ArrayPoint::new([0.0, 0.0]),
///     ArrayPoint::new([1.0, 0.0]),
/// ).unwrap());
/// path.push_line(LineSegment::new(
///     ArrayPoint::new([1.0, 0.0]),
///     ArrayPoint::new([2.0, 0.0]),
/// ).unwrap());
///
/// assert_eq!(path.start().unwrap(), ArrayPoint::new([0.0, 0.0]));
/// assert_eq!(path.end().unwrap(), ArrayPoint::new([2.0, 0.0]));
///
/// let t = path.t_at_length(1.0, 1e-12).unwrap().unwrap();
/// let p = path.point_at_length(1.0, 1e-12).unwrap().unwrap();
/// let tangent = path.tangent(0.25).unwrap();
/// # let _ = (t, p, tangent);
/// ```
pub struct BezierPath<P, const N: usize>
where
    P: Point,
{
    segments: ArrayVec<[BezierSegment<P>; N]>,
}

impl<P: Point, const N: usize> Clone for BezierPath<P, N>
where
    [BezierSegment<P>; N]: tinyvec::Array<Item = BezierSegment<P>>,
{
    fn clone(&self) -> Self {
        Self {
            segments: self.segments,
        }
    }
}

impl<P: Point + core::fmt::Debug, const N: usize> core::fmt::Debug for BezierPath<P, N>
where
    [BezierSegment<P>; N]: tinyvec::Array<Item = BezierSegment<P>>,
{
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("BezierPath")
            .field("segments", &self.segments.as_slice())
            .finish()
    }
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

    /// Return the fixed segment capacity.
    pub fn capacity(&self) -> usize {
        N
    }

    /// Return the number of additional segments that fit.
    pub fn remaining_capacity(&self) -> usize {
        N - self.len()
    }

    /// Remove all segments while retaining capacity.
    pub fn clear(&mut self) {
        self.segments.clear();
    }

    /// Remove and return the final segment, or `None` when empty.
    pub fn pop(&mut self) -> Option<BezierSegment<P>> {
        self.segments.pop()
    }

    /// Append a segment, returning it unchanged if the path is full.
    pub fn try_push(&mut self, segment: BezierSegment<P>) -> Result<(), BezierSegment<P>> {
        if self.len() == N {
            Err(segment)
        } else {
            self.segments.push(segment);
            Ok(())
        }
    }

    /// Iterate over the stored segments.
    pub fn segments(&self) -> slice::Iter<'_, BezierSegment<P>> {
        self.segments.iter()
    }

    /// Push a segment, returning false if the path is at capacity.
    pub fn push(&mut self, segment: BezierSegment<P>) -> bool {
        self.try_push(segment).is_ok()
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

    /// Evaluate the derivative with respect to the normalized path parameter.
    /// Values outside `[0, 1]` are clamped before taking the selected piece's
    /// derivative. Interior boundaries select the right-hand piece; the final
    /// endpoint selects the left. No derivative across a jump is implied.
    /// Returns `None` for an empty path or NaN parameter.
    pub fn derivative_at(&self, t: P::Scalar) -> Option<P> {
        let (index, local) = self.segment_parameter(t)?;
        let start = crate::path_parameter::boundary::<P::Scalar>(index, self.len());
        let end = crate::path_parameter::boundary::<P::Scalar>(index + 1, self.len());
        Some(crate::numeric::divide(
            self.segments[index].derivative_at(local),
            end - start,
        ))
    }

    /// Return the selected piece's speed per unit of normalized path parameter.
    pub fn speed(&self, t: P::Scalar) -> Option<P::Scalar>
    where
        P: PointNorm,
    {
        self.derivative_at(t).map(|value| value.norm())
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

    /// Return the sum of segment lengths using adaptive Gauss-Legendre quadrature.
    ///
    /// Empty paths have length zero. Gaps between segments do not add length.
    /// The positive finite `accuracy` is an absolute tolerance shared across segments.
    pub fn arclen(&self, accuracy: P::Scalar) -> Result<P::Scalar, ArcLengthError>
    where
        P: PointNorm,
    {
        let zero = <P::Scalar as NumCast>::from(0.0).unwrap();
        let one = <P::Scalar as NumCast>::from(1.0).unwrap();
        self.arclen_between(zero, one, accuracy)
    }

    /// Integrate length over `start..end` in the normalized path domain `[0, 1]`.
    /// Gaps do not contribute length; an empty path returns zero.
    pub fn arclen_between(
        &self,
        start: P::Scalar,
        end: P::Scalar,
        accuracy: P::Scalar,
    ) -> Result<P::Scalar, ArcLengthError>
    where
        P: PointNorm,
    {
        validate_accuracy(accuracy)?;
        let zero = <P::Scalar as NumCast>::from(0.0).unwrap();
        let one = <P::Scalar as NumCast>::from(1.0).unwrap();
        if !start.is_finite() || !end.is_finite() || start < zero || end > one || start > end {
            return Err(ArcLengthError::InvalidParameter);
        }
        if self.is_empty() || start == end {
            return Ok(zero);
        }
        let count = <P::Scalar as NumCast>::from(self.len()).unwrap();
        let tolerance =
            (accuracy / count).max(P::Scalar::min_positive_value() * P::Scalar::epsilon());
        let mut total = zero;
        for (i, segment) in self.segments.iter().enumerate() {
            if let Some((a, b)) = crate::path_parameter::local_range(start, end, i, self.len()) {
                total = total + segment.arclen_between(a, b, tolerance)?;
            }
        }
        if total.is_finite() {
            Ok(total)
        } else {
            Err(ArcLengthError::NonFiniteComputation)
        }
    }

    /// Approximate length using `segments` lines per curve (per knot span for B-splines).
    /// Disconnected pieces contribute only their own lengths.
    pub fn arclen_polyline(&self, segments: usize) -> P::Scalar
    where
        P: PointNorm,
    {
        self.segments.iter().fold(
            <P::Scalar as NumCast>::from(0.0).unwrap(),
            |sum, segment| sum + segment.arclen_polyline(segments),
        )
    }

    /// Return a path parameter at length `s`, clamped to the endpoints.
    /// Returns `Ok(None)` for an empty path. At a discontinuity the returned
    /// representable parameter retains the selected segment's ownership.
    pub fn t_at_length(
        &self,
        s: P::Scalar,
        accuracy: P::Scalar,
    ) -> Result<Option<P::Scalar>, ArcLengthError>
    where
        P: PointNorm,
    {
        match self.location_at_length(s, accuracy)? {
            Some((index, local)) => {
                crate::path_parameter::encode(index, local, self.len()).map(Some)
            }
            None => Ok(None),
        }
    }

    fn location_at_length(
        &self,
        s: P::Scalar,
        accuracy: P::Scalar,
    ) -> Result<Option<(usize, P::Scalar)>, ArcLengthError>
    where
        P: PointNorm,
    {
        validate_accuracy(accuracy)?;
        if !s.is_finite() {
            return Err(ArcLengthError::InvalidParameter);
        }
        if self.is_empty() {
            return Ok(None);
        }
        let zero = <P::Scalar as NumCast>::from(0.0).unwrap();
        let one = <P::Scalar as NumCast>::from(1.0).unwrap();
        if s <= zero {
            return Ok(Some((0, zero)));
        }
        let count = <P::Scalar as NumCast>::from(self.len()).unwrap();
        let tolerance =
            (accuracy / (count + one)).max(P::Scalar::min_positive_value() * P::Scalar::epsilon());
        let mut remaining = s;
        let mut any_length = false;
        for (index, segment) in self.segments.iter().enumerate() {
            let length = segment.arclen(tolerance)?;
            any_length |= length > zero;
            if length > zero && remaining < length {
                let local = segment.t_at_length(remaining, tolerance)?;
                return Ok(Some((index, local)));
            }
            remaining = remaining - length;
        }
        Ok(Some(if any_length {
            (self.len() - 1, one)
        } else {
            (0, zero)
        }))
    }

    /// Evaluate a point at length `s`. Returns `Ok(None)` for an empty path.
    pub fn point_at_length(
        &self,
        s: P::Scalar,
        accuracy: P::Scalar,
    ) -> Result<Option<P>, ArcLengthError>
    where
        P: PointNorm,
    {
        match self.location_at_length(s, accuracy)? {
            Some((index, local)) => Ok(Some(self.segments[index].eval(local))),
            None => Ok(None),
        }
    }

    /// Return the union of tight segment boxes, or `None` for an empty path.
    pub fn bounding_box(&self) -> Option<BoundingBox<P::Scalar, { P::DIM }>>
    where
        P: PointIndex,
        [P::Scalar; 1]: tinyvec::Array<Item = P::Scalar>,
        [P::Scalar; 2]: tinyvec::Array<Item = P::Scalar>,
        [P::Scalar; 3]: tinyvec::Array<Item = P::Scalar>,
        [P::Scalar; 4]: tinyvec::Array<Item = P::Scalar>,
    {
        let mut segments = self.segments.iter();
        let mut bounds = segments.next()?.bounding_box();
        for segment in segments {
            bounds = bounds.union(&segment.bounding_box());
        }
        Some(bounds)
    }

    /// Return a conservative box around all segment control points.
    /// Returns `None` for an empty path.
    pub fn control_box(&self) -> Option<BoundingBox<P::Scalar, { P::DIM }>>
    where
        P: PointIndex,
    {
        let mut segments = self.segments.iter();
        let mut bounds = segments.next()?.control_box();
        for segment in segments {
            bounds = bounds.union(&segment.control_box());
        }
        Some(bounds)
    }

    fn segment_parameter(&self, t: P::Scalar) -> Option<(usize, P::Scalar)> {
        crate::path_parameter::decode(t, self.segments.len())
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
    use super::*;
    use crate::ArrayPoint;
    use crate::test_support::{assert_approx_equal, assert_point_approx_equal};

    const ABSOLUTE_TOLERANCE: f64 = 1e-12;
    const RELATIVE_TOLERANCE: f64 = 1e-12;

    #[test]
    fn bezier_path_eval_segments() {
        let mut path: BezierPath<ArrayPoint<f64, 2>, 4> = BezierPath::new();
        path.push_line(
            LineSegment::new(ArrayPoint::new([0.0, 0.0]), ArrayPoint::new([1.0, 0.0])).unwrap(),
        );
        path.push_line(
            LineSegment::new(ArrayPoint::new([1.0, 0.0]), ArrayPoint::new([1.0, 1.0])).unwrap(),
        );

        let p0 = path.eval(0.0).unwrap();
        assert_eq!(p0, ArrayPoint::new([0.0, 0.0]));

        let p1 = path.eval(0.25).unwrap();
        assert_point_approx_equal(
            p1,
            ArrayPoint::new([0.5, 0.0]),
            ABSOLUTE_TOLERANCE,
            RELATIVE_TOLERANCE,
        );

        let p2 = path.eval(0.75).unwrap();
        assert_point_approx_equal(
            p2,
            ArrayPoint::new([1.0, 0.5]),
            ABSOLUTE_TOLERANCE,
            RELATIVE_TOLERANCE,
        );
    }

    #[test]
    fn bezier_path_bounds_union() {
        let mut path: BezierPath<ArrayPoint<f64, 2>, 4> = BezierPath::new();
        path.push_line(
            LineSegment::new(ArrayPoint::new([0.0, 0.0]), ArrayPoint::new([1.0, 0.0])).unwrap(),
        );
        path.push_line(
            LineSegment::new(ArrayPoint::new([1.0, 0.0]), ArrayPoint::new([1.0, 1.0])).unwrap(),
        );

        let bounds = path.bounding_box().unwrap();
        assert_eq!(bounds[0], (0.0, 1.0));
        assert_eq!(bounds[1], (0.0, 1.0));
    }

    #[test]
    fn bezier_path_clamps_out_of_range() {
        let mut path: BezierPath<ArrayPoint<f64, 2>, 2> = BezierPath::new();
        path.push_line(
            LineSegment::new(ArrayPoint::new([0.0, 0.0]), ArrayPoint::new([1.0, 0.0])).unwrap(),
        );
        path.push_line(
            LineSegment::new(ArrayPoint::new([1.0, 0.0]), ArrayPoint::new([1.0, 1.0])).unwrap(),
        );

        let start = path.eval(-1.0).unwrap();
        assert_eq!(start, ArrayPoint::new([0.0, 0.0]));

        let end = path.eval(2.0).unwrap();
        assert_eq!(end, ArrayPoint::new([1.0, 1.0]));
    }

    #[test]
    fn bezier_path_capacity() {
        let mut path: BezierPath<ArrayPoint<f64, 2>, 1> = BezierPath::new();
        let first = path.push_line(
            LineSegment::new(ArrayPoint::new([0.0, 0.0]), ArrayPoint::new([1.0, 0.0])).unwrap(),
        );
        let second = path.push_line(
            LineSegment::new(ArrayPoint::new([1.0, 0.0]), ArrayPoint::new([2.0, 0.0])).unwrap(),
        );

        assert!(first);
        assert!(!second);
        assert_eq!(path.len(), 1);
    }

    #[test]
    fn bezier_path_api_parity() {
        let mut path: BezierPath<ArrayPoint<f64, 2>, 2> = BezierPath::new();
        path.push_line(
            LineSegment::new(ArrayPoint::new([0.0, 0.0]), ArrayPoint::new([1.0, 0.0])).unwrap(),
        );
        path.push_line(
            LineSegment::new(ArrayPoint::new([1.0, 0.0]), ArrayPoint::new([2.0, 0.0])).unwrap(),
        );

        assert_eq!(path.start().unwrap(), ArrayPoint::new([0.0, 0.0]));
        assert_eq!(path.end().unwrap(), ArrayPoint::new([2.0, 0.0]));

        let reversed = path.reverse();
        let p0 = path.eval(0.25).unwrap();
        let p1 = reversed.eval(0.75).unwrap();
        assert_point_approx_equal(p0, p1, ABSOLUTE_TOLERANCE, RELATIVE_TOLERANCE);

        let tangent = path.tangent(0.25).unwrap();
        assert_approx_equal(tangent[0], 1.0, ABSOLUTE_TOLERANCE, RELATIVE_TOLERANCE);
        assert_approx_equal(tangent[1], 0.0, ABSOLUTE_TOLERANCE, RELATIVE_TOLERANCE);

        let curvature = path.curvature(0.25).unwrap();
        assert_approx_equal(curvature, 0.0, ABSOLUTE_TOLERANCE, RELATIVE_TOLERANCE);
        assert!(path.normal(0.25).is_none());

        let t = path.t_at_length(1.0, 1e-12).unwrap().unwrap();
        assert!((t - 0.5).abs() < 1e-6);

        let p = path.point_at_length(1.0, 1e-12).unwrap().unwrap();
        assert_point_approx_equal(p, ArrayPoint::new([1.0, 0.0]), 1e-9, 1e-9);
    }

    #[test]
    fn small_scale_path_remains_non_degenerate() {
        let scale = 1e-20_f64;
        let mut path: BezierPath<ArrayPoint<f64, 2>, 1> = BezierPath::new();
        assert!(path.push_line(
            LineSegment::new(ArrayPoint::new([0.0, 0.0]), ArrayPoint::new([scale, 0.0]),).unwrap()
        ));

        assert_approx_equal(
            path.t_at_length(scale * 0.5, 1e-12).unwrap().unwrap(),
            0.5,
            1e-9,
            1e-9,
        );
    }
}
