//! B-spline path utilities.

use crate::arc_length::validate_accuracy;
use crate::{ArcLengthError, BoundingBox};
use core::slice;

use super::{ArrayVec, BSpline, Point, PointDot, PointIndex, PointNorm};
use num_traits::{Float, NumCast};

/// A path composed of multiple B-spline segments of the same degree/knots layout.
///
/// Component-aware helpers (e.g. bounding boxes) require `PointIndex`.
///
/// # Examples
/// ```rust
/// use stroke::{BSpline, BSplinePath, PointN};
///
/// let knots: [f64; 4] = [0.0, 0.0, 1.0, 1.0];
/// let seg1: BSpline<PointN<f64, 2>, 4, 2, 1> = BSpline::new(
///     knots,
///     [
///         PointN::new([0.0, 0.0]),
///         PointN::new([1.0, 0.0]),
///     ],
/// )
/// .unwrap();
/// let seg2: BSpline<PointN<f64, 2>, 4, 2, 1> = BSpline::new(
///     knots,
///     [
///         PointN::new([1.0, 0.0]),
///         PointN::new([2.0, 0.0]),
///     ],
/// )
/// .unwrap();
///
/// let mut path: BSplinePath<PointN<f64, 2>, 4, 2, 1, 2> = BSplinePath::new();
/// path.push(seg1);
/// path.push(seg2);
///
/// assert_eq!(path.start().unwrap(), PointN::new([0.0, 0.0]));
/// assert_eq!(path.end().unwrap(), PointN::new([2.0, 0.0]));
///
/// let t = path.t_at_length(1.0, 1e-12).unwrap().unwrap();
/// let p = path.point_at_length(1.0, 1e-12).unwrap().unwrap();
/// let tangent = path.tangent(0.25).unwrap();
/// # let _ = (t, p, tangent);
/// ```
pub struct BSplinePath<P, const K: usize, const C: usize, const D: usize, const N: usize>
where
    P: Point,
{
    segments: ArrayVec<[BSpline<P, K, C, D>; N]>,
}

impl<P: Point, const K: usize, const C: usize, const D: usize, const N: usize> Clone
    for BSplinePath<P, K, C, D, N>
where
    [BSpline<P, K, C, D>; N]: tinyvec::Array<Item = BSpline<P, K, C, D>>,
{
    fn clone(&self) -> Self {
        Self {
            segments: self.segments,
        }
    }
}

impl<P, const K: usize, const C: usize, const D: usize, const N: usize> core::fmt::Debug
    for BSplinePath<P, K, C, D, N>
where
    P: Point + core::fmt::Debug,
    P::Scalar: core::fmt::Debug,
    [BSpline<P, K, C, D>; N]: tinyvec::Array<Item = BSpline<P, K, C, D>>,
{
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("BSplinePath")
            .field("segments", &self.segments.as_slice())
            .finish()
    }
}

impl<P, const K: usize, const C: usize, const D: usize, const N: usize> BSplinePath<P, K, C, D, N>
where
    P: Point,
    [BSpline<P, K, C, D>; N]: tinyvec::Array<Item = BSpline<P, K, C, D>>,
{
    /// Create an empty path with capacity `N`.
    pub fn new() -> Self {
        BSplinePath {
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
    pub fn pop(&mut self) -> Option<BSpline<P, K, C, D>> {
        self.segments.pop()
    }

    /// Append a segment, returning it unchanged if the path is full.
    pub fn try_push(&mut self, segment: BSpline<P, K, C, D>) -> Result<(), BSpline<P, K, C, D>> {
        if self.len() == N {
            Err(segment)
        } else {
            self.segments.push(segment);
            Ok(())
        }
    }

    /// Iterate over the stored segments.
    pub fn segments(&self) -> slice::Iter<'_, BSpline<P, K, C, D>> {
        self.segments.iter()
    }

    /// Push a segment, returning false if the path is at capacity.
    pub fn push(&mut self, segment: BSpline<P, K, C, D>) -> bool {
        self.try_push(segment).is_ok()
    }

    /// Evaluate a point along the path for `t` in `[0, 1]`.
    /// Values outside the range are clamped. Returns None for empty paths.
    pub fn eval(&self, t: P::Scalar) -> Option<P>
    where
        [(); D + 1]: Sized,
    {
        let (index, local_t) = self.segment_parameter(t)?;
        let segment = &self.segments[index];
        let (kmin, kmax) = segment.knot_domain();
        let one = <P::Scalar as NumCast>::from(1.0).unwrap();
        let t_mapped = kmin * (one - local_t) + kmax * local_t;
        segment.eval(t_mapped).ok()
    }

    /// Return the start point of the path.
    pub fn start(&self) -> Option<P>
    where
        [(); D + 1]: Sized,
    {
        let zero = <P::Scalar as NumCast>::from(0.0).unwrap();
        self.eval(zero)
    }

    /// Return the end point of the path.
    pub fn end(&self) -> Option<P>
    where
        [(); D + 1]: Sized,
    {
        let one = <P::Scalar as NumCast>::from(1.0).unwrap();
        self.eval(one)
    }

    /// Return a path with reversed direction.
    pub fn reverse(&self) -> Self {
        let mut reversed = BSplinePath::new();
        for segment in self.segments.iter().rev() {
            let _ = reversed.push(segment.reverse());
        }
        reversed
    }

    /// Evaluate the derivative with respect to the normalized path parameter.
    /// Includes both native knot-domain and path segment scaling. Values outside
    /// `[0, 1]` are clamped before taking the selected piece's derivative.
    /// Interior boundaries select the right-hand piece; the final endpoint uses
    /// the left. No derivative across a jump is implied.
    /// Returns `None` for an empty path, NaN, or failed differential evaluation.
    pub fn derivative_at(&self, t: P::Scalar) -> Option<P>
    where
        [(); D + 1]: Sized,
    {
        let (index, local) = self.segment_parameter(t)?;
        let segment = &self.segments[index];
        let (start, end) = segment.knot_domain();
        let one = <P::Scalar as NumCast>::from(1.0).unwrap();
        let parameter = start * (one - local) + end * local;
        let derivative = segment.derivative_at(parameter).ok()?;
        let width = end - start;
        let scaled = if width.is_finite() {
            derivative * width
        } else {
            let half = <P::Scalar as NumCast>::from(0.5).unwrap();
            (derivative * (end * half - start * half)) * (one + one)
        };
        let path_start = crate::path_parameter::boundary::<P::Scalar>(index, self.len());
        let path_end = crate::path_parameter::boundary::<P::Scalar>(index + 1, self.len());
        let result = crate::numeric::divide(scaled, path_end - path_start);
        result.is_finite().then_some(result)
    }

    /// Return the selected piece's speed per unit of normalized path parameter.
    pub fn speed(&self, t: P::Scalar) -> Option<P::Scalar>
    where
        P: PointNorm,
        [(); D + 1]: Sized,
    {
        let speed = self.derivative_at(t)?.norm();
        speed.is_finite().then_some(speed)
    }

    /// Return the unit tangent direction at `t`.
    pub fn tangent(&self, t: P::Scalar) -> Option<P>
    where
        P: PointNorm,
        [(); D - 1]: Sized,
        [(); C - 1]: Sized,
        [(); K - 2]: Sized,
        [(); D]: Sized,
        [(); D + 1]: Sized,
        [(); (D - 1) + 1]: Sized, // keep (D - 1) + 1; compiler doesn't normalize to D
    {
        let (index, local_t) = self.segment_parameter(t)?;
        let segment = &self.segments[index];
        let (kmin, kmax) = segment.knot_domain();
        let one = <P::Scalar as NumCast>::from(1.0).unwrap();
        let t_mapped = kmin * (one - local_t) + kmax * local_t;
        segment.tangent(t_mapped).ok()
    }

    /// Return the curvature magnitude at `t`.
    pub fn curvature(&self, t: P::Scalar) -> Option<P::Scalar>
    where
        P: PointNorm + PointDot,
        [(); D + 1]: Sized,
        [(); D - 2]: Sized,
    {
        let (index, local_t) = self.segment_parameter(t)?;
        let segment = &self.segments[index];
        let (kmin, kmax) = segment.knot_domain();
        let one = <P::Scalar as NumCast>::from(1.0).unwrap();
        let t_mapped = kmin * (one - local_t) + kmax * local_t;
        segment.curvature(t_mapped).ok()
    }

    /// Return the principal normal direction at `t`.
    ///
    /// Returns `None` if the path is empty or curvature is undefined.
    pub fn normal(&self, t: P::Scalar) -> Option<P>
    where
        P: PointNorm + PointDot,
        [(); D + 1]: Sized,
        [(); D - 2]: Sized,
    {
        let (index, local_t) = self.segment_parameter(t)?;
        let segment = &self.segments[index];
        let (kmin, kmax) = segment.knot_domain();
        let one = <P::Scalar as NumCast>::from(1.0).unwrap();
        let t_mapped = kmin * (one - local_t) + kmax * local_t;
        segment.normal(t_mapped).ok().flatten()
    }

    /// Return the sum of segment lengths using adaptive Gauss-Legendre quadrature.
    ///
    /// Empty paths have length zero. Gaps between segments do not add length.
    /// The positive finite `accuracy` is an absolute tolerance shared across segments.
    pub fn arclen(&self, accuracy: P::Scalar) -> Result<P::Scalar, ArcLengthError>
    where
        P: PointNorm,
        [(); D + 1]: Sized,
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
        [(); D + 1]: Sized,
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
                let (kmin, kmax) = segment.knot_domain();
                let a = kmin * (one - a) + kmax * a;
                let b = kmin * (one - b) + kmax * b;
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
        [(); D + 1]: Sized,
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
        [(); D + 1]: Sized,
    {
        match self.location_at_length(s, accuracy)? {
            Some((index, local)) => {
                let (start, end) = self.segments[index].knot_domain();
                let local = crate::numeric::unit_parameter(local, start, end);
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
        [(); D + 1]: Sized,
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
            return Ok(Some((0, self.segments[0].knot_domain().0)));
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
            (
                self.len() - 1,
                self.segments[self.len() - 1].knot_domain().1,
            )
        } else {
            (0, self.segments[0].knot_domain().0)
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
        [(); D + 1]: Sized,
    {
        match self.location_at_length(s, accuracy)? {
            Some((index, local)) => self.segments[index]
                .eval(local)
                .map(Some)
                .map_err(|_| ArcLengthError::CurveEvaluationFailed),
            None => Ok(None),
        }
    }

    /// Return the union of tight segment boxes, or `None` for an empty path.
    pub fn bounding_box(&self) -> Option<BoundingBox<P::Scalar, { P::DIM }>>
    where
        P: PointIndex,
        [(); D + 1]: Sized,
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

impl<P, const K: usize, const C: usize, const D: usize, const N: usize> Default
    for BSplinePath<P, K, C, D, N>
where
    P: Point,
    [BSpline<P, K, C, D>; N]: tinyvec::Array<Item = BSpline<P, K, C, D>>,
{
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::PointN;
    use crate::test_support::{assert_approx_equal, assert_point_approx_equal};

    const ABSOLUTE_TOLERANCE: f64 = 1e-12;
    const RELATIVE_TOLERANCE: f64 = 1e-12;

    #[test]
    fn bspline_path_eval_segments() {
        let knots: [f64; 4] = [0.0, 0.0, 1.0, 1.0];
        let s1: BSpline<PointN<f64, 2>, 4, 2, 1> =
            BSpline::new(knots, [PointN::new([0.0, 0.0]), PointN::new([1.0, 0.0])]).unwrap();
        let s2: BSpline<PointN<f64, 2>, 4, 2, 1> =
            BSpline::new(knots, [PointN::new([1.0, 0.0]), PointN::new([1.0, 1.0])]).unwrap();

        let mut path: BSplinePath<PointN<f64, 2>, 4, 2, 1, 4> = BSplinePath::new();
        path.push(s1);
        path.push(s2);

        let p1 = path.eval(0.25).unwrap();
        assert_point_approx_equal(
            p1,
            PointN::new([0.5, 0.0]),
            ABSOLUTE_TOLERANCE,
            RELATIVE_TOLERANCE,
        );

        let p2 = path.eval(0.75).unwrap();
        assert_point_approx_equal(
            p2,
            PointN::new([1.0, 0.5]),
            ABSOLUTE_TOLERANCE,
            RELATIVE_TOLERANCE,
        );
    }

    #[test]
    fn bspline_path_clamps_out_of_range() {
        let knots: [f64; 4] = [0.0, 0.0, 1.0, 1.0];
        let s1: BSpline<PointN<f64, 2>, 4, 2, 1> =
            BSpline::new(knots, [PointN::new([0.0, 0.0]), PointN::new([1.0, 0.0])]).unwrap();
        let s2: BSpline<PointN<f64, 2>, 4, 2, 1> =
            BSpline::new(knots, [PointN::new([1.0, 0.0]), PointN::new([1.0, 1.0])]).unwrap();

        let mut path: BSplinePath<PointN<f64, 2>, 4, 2, 1, 4> = BSplinePath::new();
        path.push(s1);
        path.push(s2);

        let start = path.eval(-1.0).unwrap();
        assert_eq!(start, PointN::new([0.0, 0.0]));

        let end = path.eval(2.0).unwrap();
        assert_eq!(end, PointN::new([1.0, 1.0]));
    }

    #[test]
    fn bspline_path_bounds_control_points() {
        let knots: [f64; 4] = [0.0, 0.0, 1.0, 1.0];
        let s1: BSpline<PointN<f64, 2>, 4, 2, 1> =
            BSpline::new(knots, [PointN::new([0.0, 2.0]), PointN::new([1.0, -1.0])]).unwrap();
        let s2: BSpline<PointN<f64, 2>, 4, 2, 1> =
            BSpline::new(knots, [PointN::new([-1.0, 0.5]), PointN::new([2.0, -2.0])]).unwrap();

        let mut path: BSplinePath<PointN<f64, 2>, 4, 2, 1, 4> = BSplinePath::new();
        path.push(s1);
        path.push(s2);

        let bounds = path.bounding_box().unwrap();
        assert_eq!(bounds[0], (-1.0, 2.0));
        assert_eq!(bounds[1], (-2.0, 2.0));
    }

    #[test]
    fn bspline_path_api_parity() {
        let knots: [f64; 4] = [0.0, 0.0, 1.0, 1.0];
        let s1: BSpline<PointN<f64, 2>, 4, 2, 1> =
            BSpline::new(knots, [PointN::new([0.0, 0.0]), PointN::new([1.0, 0.0])]).unwrap();
        let s2: BSpline<PointN<f64, 2>, 4, 2, 1> =
            BSpline::new(knots, [PointN::new([1.0, 0.0]), PointN::new([2.0, 0.0])]).unwrap();

        let mut path: BSplinePath<PointN<f64, 2>, 4, 2, 1, 2> = BSplinePath::new();
        path.push(s1);
        path.push(s2);

        assert_eq!(path.start().unwrap(), PointN::new([0.0, 0.0]));
        assert_eq!(path.end().unwrap(), PointN::new([2.0, 0.0]));

        let reversed = path.reverse();
        let p0 = path.eval(0.25).unwrap();
        let p1 = reversed.eval(0.75).unwrap();
        assert_point_approx_equal(p0, p1, ABSOLUTE_TOLERANCE, RELATIVE_TOLERANCE);

        let tangent = path.tangent(0.25).unwrap();
        assert_approx_equal(tangent[0], 1.0, ABSOLUTE_TOLERANCE, RELATIVE_TOLERANCE);
        assert_approx_equal(tangent[1], 0.0, ABSOLUTE_TOLERANCE, RELATIVE_TOLERANCE);

        let t = path.t_at_length(1.0, 1e-12).unwrap().unwrap();
        assert!((t - 0.5).abs() < 1e-6);

        let p = path.point_at_length(1.0, 1e-12).unwrap().unwrap();
        assert_point_approx_equal(p, PointN::new([1.0, 0.0]), 1e-9, 1e-9);
    }

    #[test]
    fn bspline_path_curvature_nonzero() {
        let knots: [f64; 6] = [0.0, 0.0, 0.0, 1.0, 1.0, 1.0];
        let seg: BSpline<PointN<f64, 2>, 6, 3, 2> = BSpline::new(
            knots,
            [
                PointN::new([0.0, 0.0]),
                PointN::new([1.0, 1.0]),
                PointN::new([2.0, 0.0]),
            ],
        )
        .unwrap();

        let mut path: BSplinePath<PointN<f64, 2>, 6, 3, 2, 1> = BSplinePath::new();
        path.push(seg);

        let curvature = path.curvature(0.5).unwrap();
        assert!(curvature > 0.0);

        let normal = path.normal(0.5).unwrap();
        assert!((normal.squared_norm() - 1.0).abs() < 1e-6);
    }

    #[test]
    fn small_scale_path_remains_non_degenerate() {
        let scale = 1e-20_f64;
        let knots = [0.0, 0.0, 1.0, 1.0];
        let segment = BSpline::<PointN<f64, 1>, 4, 2, 1>::new(
            knots,
            [PointN::new([0.0]), PointN::new([scale])],
        )
        .unwrap();
        let mut path: BSplinePath<PointN<f64, 1>, 4, 2, 1, 1> = BSplinePath::new();
        assert!(path.push(segment));

        assert_approx_equal(
            path.t_at_length(scale * 0.5, 1e-12).unwrap().unwrap(),
            0.5,
            1e-9,
            1e-9,
        );
    }
}
