//! Line segment curve type.

use super::{ArcLengthError, ArrayVec, BoundingBox, Point, PointIndex, PointNorm};
use crate::arc_length::validate_accuracy;
use crate::error::{CurveConstructionError, validate_control_points};
use crate::numeric::is_negligible;
use num_traits::{Float, NumCast};

/// A line segment between `start` and `end`.
///
/// `eval(t)` interpolates between the endpoints for `t` in `[0, 1]`.
/// This is also a degree-one Bezier curve.
///
/// Bounding boxes require `PointIndex`. Distance and differential geometry
/// methods require `PointNorm`.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct LineSegment<P> {
    pub(crate) start: P,
    pub(crate) end: P,
}

impl<P> LineSegment<P>
where
    P: Point,
{
    /// Create a new line segment from `start` to `end`.
    ///
    /// # Errors
    ///
    /// Returns [`CurveConstructionError::NonFiniteControlPoint`] with the index
    /// of the first control point containing a non-finite coordinate.
    pub fn new(start: P, end: P) -> Result<Self, CurveConstructionError> {
        validate_control_points(&[start, end])?;
        Ok(LineSegment { start, end })
    }

    /// Borrow the control points in curve order.
    pub fn control_points(&self) -> core::array::IntoIter<&P, 2> {
        [&self.start, &self.end].into_iter()
    }

    /// Return the control-point array, consuming the curve.
    pub fn into_control_points(self) -> [P; 2] {
        [self.start, self.end]
    }

    /// Replace all control points. Leave the curve unchanged if validation fails.
    pub fn set_control_points(&mut self, points: [P; 2]) -> Result<(), CurveConstructionError> {
        *self = Self::new(points[0], points[1])?;
        Ok(())
    }

    /// Return a curve with mapped control points. Reject non-finite results.
    /// Affine maps, such as translations, also map every point of the curve.
    pub fn try_map_control_points(
        &self,
        map: impl FnMut(P) -> P,
    ) -> Result<Self, CurveConstructionError> {
        let points = (*self).into_control_points().map(map);
        Self::new(points[0], points[1])
    }

    /// Evaluate the derivative with respect to `t` in `[0, 1]`.
    pub fn derivative_at(&self, t: P::Scalar) -> P {
        crate::Bezier::from(*self).derivative_at(t)
    }

    /// Return the speed, the magnitude of the derivative with respect to `t`.
    pub fn speed(&self, t: P::Scalar) -> P::Scalar
    where
        P: PointNorm,
    {
        self.derivative_at(t).norm()
    }

    /// Extract `start..end` and give the result a new parameter range `[0, 1]`.
    /// Requires finite `0 <= start <= end <= 1`. Equal endpoints give a constant curve.
    pub fn subcurve(
        &self,
        range: core::ops::Range<P::Scalar>,
    ) -> Result<Self, crate::CurveParameterError> {
        crate::Bezier::from(*self).subcurve(range).map(Self::from)
    }

    /// Split at a distance from the start, clamping finite distances to the endpoints.
    pub fn split_at_length(
        &self,
        length: P::Scalar,
        accuracy: P::Scalar,
    ) -> Result<(Self, Self), ArcLengthError>
    where
        P: PointNorm,
    {
        let (left, right) = crate::Bezier::from(*self).split_at_length(length, accuracy)?;
        Ok((Self::from(left), Self::from(right)))
    }

    /// Return the start point of the segment.
    pub fn start(&self) -> P {
        let zero = <P::Scalar as NumCast>::from(0.0).unwrap();
        self.eval(zero)
    }

    /// Return the end point of the segment.
    pub fn end(&self) -> P {
        let one = <P::Scalar as NumCast>::from(1.0).unwrap();
        self.eval(one)
    }

    /// Return a segment with reversed direction.
    pub fn reverse(&self) -> Self {
        LineSegment {
            start: self.end,
            end: self.start,
        }
    }

    /// Evaluate a point along the segment for `t` in `[0, 1]`.
    pub fn eval(&self, t: P::Scalar) -> P {
        crate::numeric::lerp(self.start, self.end, t)
    }

    /// Split the segment at `t` into two sub-segments.
    pub fn split(&self, t: P::Scalar) -> (Self, Self) {
        let ctrl_ab = self.eval(t);

        (
            LineSegment {
                start: self.start,
                end: ctrl_ab,
            },
            LineSegment {
                start: ctrl_ab,
                end: self.end,
            },
        )
    }

    /// Return the Euclidean distance to this segment using a scaled projection.
    /// Non-finite input points return NaN; an unrepresentable distance is infinite.
    pub fn distance_to_point(&self, p: P) -> P::Scalar
    where
        P: PointNorm,
    {
        let zero = <P::Scalar as NumCast>::from(0.0).unwrap();
        let one = <P::Scalar as NumCast>::from(1.0).unwrap();
        let half = <P::Scalar as NumCast>::from(0.5).unwrap();
        if !p.is_finite() || !self.start.is_finite() || !self.end.is_finite() {
            return P::Scalar::nan();
        }
        let start_distance = crate::distance::between(p, self.start);
        let end_distance = crate::distance::between(p, self.end);
        let best_endpoint = start_distance.min(end_distance);
        if best_endpoint == zero {
            return zero;
        }

        // Project displacements from the midpoint. On a huge segment, a small
        // offset can disappear if expressed as t=0.5 + offset. The high/low
        // midpoint pair preserves the rounding error in the center.
        let offset = |axis, scaled| {
            let (center, tail) = crate::numeric::midpoint_parts(self.start[axis], self.end[axis]);
            if scaled {
                (p[axis] * half - center * half) - tail * half
            } else {
                (p[axis] - center) - tail
            }
        };
        let direction_half = (0..P::DIM).any(|i| !(self.end[i] - self.start[i]).is_finite());
        let offset_half = (0..P::DIM).any(|i| !offset(i, false).is_finite());
        let direction = |axis| {
            if direction_half {
                self.end[axis] * half - self.start[axis] * half
            } else {
                self.end[axis] - self.start[axis]
            }
        };
        let mut direction_scale = zero;
        let mut offset_scale = zero;
        for axis in 0..P::DIM {
            direction_scale = direction_scale.max(direction(axis).abs());
            offset_scale = offset_scale.max(offset(axis, offset_half).abs());
        }
        // Compare with exact zero to preserve tiny nonzero segments.
        if direction_scale == zero {
            return start_distance;
        }
        if offset_scale == zero {
            return zero;
        }
        let mut dot = zero;
        let mut correction = zero;
        let mut squared_direction = zero;
        for axis in 0..P::DIM {
            let d = direction(axis) / direction_scale;
            let w = offset(axis, offset_half) / offset_scale;
            let product = d * w;
            let next = dot + product;
            correction = correction
                + if dot.abs() >= product.abs() {
                    (dot - next) + product
                } else {
                    (product - next) + dot
                };
            dot = next;
            squared_direction = squared_direction + d * d;
        }
        dot = dot + correction;
        // Compute twice the midpoint displacement in segment units.
        // Restore the scales after the dot products so their factors do not
        // overflow or underflow before cancellation.
        let shift = if offset_half { 1 } else { 0 } - if direction_half { 1 } else { 0 };
        let position = crate::numeric::scaled_product_ratio(
            &[offset_scale, dot.abs()],
            &[direction_scale, squared_direction],
            shift + 1,
        );
        if position >= one {
            return if dot < zero {
                start_distance
            } else {
                end_distance
            };
        }
        let mut distance = zero;
        for axis in 0..P::DIM {
            let d = direction(axis) / direction_scale;
            let native_offset = offset(axis, false);
            let residual = if native_offset.is_finite() {
                native_offset
                    - crate::numeric::scaled_product_ratio(
                        &[offset_scale, dot, d],
                        &[squared_direction],
                        if offset_half { 1 } else { 0 },
                    )
            } else {
                (offset(axis, true)
                    - crate::numeric::scaled_product_ratio(
                        &[offset_scale, dot, d],
                        &[squared_direction],
                        0,
                    ))
                    * (one + one)
            };
            distance = distance.hypot(residual);
        }
        if distance.is_nan() {
            distance
        } else {
            distance.min(best_endpoint)
        }
    }

    /// Return the derivative vector of the segment.
    pub fn derivative(&self) -> P {
        self.end - self.start
    }

    /// Return the unit tangent direction at `t`.
    pub fn tangent(&self, t: P::Scalar) -> P
    where
        P: PointNorm,
    {
        crate::Bezier::from(*self).tangent(t)
    }

    /// Return the curvature magnitude at `t`.
    ///
    /// Lines have zero curvature, so this always returns `0`.
    pub fn curvature(&self, _t: P::Scalar) -> P::Scalar
    where
        P: PointNorm,
    {
        <P::Scalar as NumCast>::from(0.0).unwrap()
    }

    /// Return the principal normal direction at `t`.
    ///
    /// Lines have zero curvature, so this always returns `None`.
    pub fn normal(&self, _t: P::Scalar) -> Option<P>
    where
        P: PointNorm,
    {
        None
    }

    /// Return the line's length. The unused subdivision count matches curved types.
    pub fn arclen_polyline(&self, _segments: usize) -> P::Scalar
    where
        P: PointNorm,
    {
        (self.end - self.start).norm()
    }

    /// Return the exact arc length of this segment.
    ///
    /// `accuracy` must be finite and positive, as for curved segments.
    /// A line needs no numerical integration.
    pub fn arclen(&self, accuracy: P::Scalar) -> Result<P::Scalar, ArcLengthError>
    where
        P: PointNorm,
    {
        let zero = <P::Scalar as NumCast>::from(0.0).unwrap();
        let one = <P::Scalar as NumCast>::from(1.0).unwrap();
        self.arclen_between(zero, one, accuracy)
    }

    /// Return the length over `start..end`, where `0 <= start <= end <= 1`.
    pub fn arclen_between(
        &self,
        start: P::Scalar,
        end: P::Scalar,
        accuracy: P::Scalar,
    ) -> Result<P::Scalar, ArcLengthError>
    where
        P: PointNorm,
    {
        let zero = <P::Scalar as NumCast>::from(0.0).unwrap();
        let one = <P::Scalar as NumCast>::from(1.0).unwrap();
        validate_accuracy(accuracy)?;
        if !start.is_finite() || !end.is_finite() || start < zero || end > one || start > end {
            return Err(ArcLengthError::InvalidParameter);
        }
        let fraction = end - start;
        let mut length = zero;
        for axis in 0..P::DIM {
            let a = self.start[axis];
            let b = self.end[axis];
            let difference = b - a;
            // Form the requested displacement before its norm. For opposite
            // signs the full subtraction may overflow even though its scaled
            // value is representable. Scaling the operands then cannot cause
            // cancellation, because they have opposite signs.
            let component = if difference.is_finite() {
                difference * fraction
            } else {
                b * fraction - a * fraction
            };
            length = length.hypot(component);
        }
        if length.is_finite() {
            Ok(length)
        } else {
            Err(ArcLengthError::NonFiniteComputation)
        }
    }

    /// Return the exact parameter `t` at arc length `s`.
    pub fn t_at_length(
        &self,
        s: P::Scalar,
        accuracy: P::Scalar,
    ) -> Result<P::Scalar, ArcLengthError>
    where
        P: PointNorm,
    {
        let length = self.arclen(accuracy)?;
        if !s.is_finite() {
            return Err(ArcLengthError::InvalidParameter);
        }
        let zero = <P::Scalar as NumCast>::from(0.0).unwrap();
        let one = <P::Scalar as NumCast>::from(1.0).unwrap();
        Ok(if length == zero {
            zero
        } else {
            (s / length).clamp(zero, one)
        })
    }

    /// Evaluate the point at arc length `s`.
    pub fn point_at_length(&self, s: P::Scalar, accuracy: P::Scalar) -> Result<P, ArcLengthError>
    where
        P: PointNorm,
    {
        Ok(self.eval(self.t_at_length(s, accuracy)?))
    }

    pub(crate) fn root(&self, a: P::Scalar, b: P::Scalar) -> ArrayVec<[P::Scalar; 1]>
    where
        [P::Scalar; 1]: tinyvec::Array<Item = P::Scalar>,
    {
        let mut r = ArrayVec::new();
        if is_negligible(a, a.abs().max(b.abs())) {
            return r;
        }
        r.push(-b / a);
        r
    }

    /// Return the tight axis-aligned bounding box of the line.
    pub fn bounding_box(&self) -> BoundingBox<P::Scalar, { P::DIM }>
    where
        P: PointIndex,
    {
        let mut bounds = [(
            <P::Scalar as NumCast>::from(0.0).unwrap(),
            <P::Scalar as NumCast>::from(0.0).unwrap(),
        ); P::DIM];

        for i in 0..P::DIM {
            if self.start[i] < self.end[i] {
                bounds[i] = (self.start[i], self.end[i]);
            } else {
                bounds[i] = (self.end[i], self.start[i]);
            }
        }

        BoundingBox::from_axes_unchecked(bounds)
    }

    /// Return the bounding box of the segment's control points.
    ///
    /// This is identical to [`Self::bounding_box`] for a line.
    pub fn control_box(&self) -> BoundingBox<P::Scalar, { P::DIM }>
    where
        P: PointIndex,
    {
        self.bounding_box()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ArrayPoint;
    use crate::CurveConstructionError;
    use crate::test_support::{assert_approx_equal, assert_point_approx_equal};

    const ABSOLUTE_TOLERANCE: f64 = 1e-12;
    const RELATIVE_TOLERANCE: f64 = 1e-12;

    #[test]
    fn rejects_non_finite_control_points() {
        let result = LineSegment::new(
            ArrayPoint::new([0.0_f64, 0.0]),
            ArrayPoint::new([f64::INFINITY, 0.0]),
        );

        assert_eq!(
            result,
            Err(CurveConstructionError::NonFiniteControlPoint { index: 1 })
        );
    }

    /// Check that the midpoint is equidistant from both endpoints.
    #[test]
    fn line_segment_interpolation() {
        let line = LineSegment {
            start: ArrayPoint::new([0f64, 1.77f64]),
            end: ArrayPoint::new([4.3f64, 3f64]),
        };

        let mid = line.eval(0.5);
        assert_approx_equal(
            (mid - line.start).squared_norm(),
            (mid - line.end).squared_norm(),
            ABSOLUTE_TOLERANCE,
            RELATIVE_TOLERANCE,
        );
    }

    /// Check distances using a 3-4-5 right triangle.
    #[test]
    fn line_segment_distance_to_point() {
        let line = LineSegment {
            start: ArrayPoint::new([0f64, 1f64, 0f64]),
            end: ArrayPoint::new([3f64, 1f64, 0f64]),
        };
        // dist to start should be 4; dist to end should be 5
        let p1 = ArrayPoint::new([0f64, 5f64, 0f64]);
        assert_approx_equal(
            line.distance_to_point(p1),
            4.0,
            ABSOLUTE_TOLERANCE,
            RELATIVE_TOLERANCE,
        );
        assert_approx_equal(
            (p1 - line.start).squared_norm().sqrt(),
            4.0,
            ABSOLUTE_TOLERANCE,
            RELATIVE_TOLERANCE,
        );
        assert_approx_equal(
            (p1 - line.end).squared_norm().sqrt(),
            5.0,
            ABSOLUTE_TOLERANCE,
            RELATIVE_TOLERANCE,
        );
        // dist to midpoint (t=0.5) should be 1
        let p2 = ArrayPoint::new([1.5f64, 2f64, 0f64]);
        assert_approx_equal(
            (p2 - line.eval(0.5)).squared_norm().sqrt(),
            1.0,
            ABSOLUTE_TOLERANCE,
            RELATIVE_TOLERANCE,
        );
    }

    #[test]
    fn line_segment_split_midpoint() {
        let line = LineSegment {
            start: ArrayPoint::new([0f64, 0f64]),
            end: ArrayPoint::new([4f64, 0f64]),
        };

        let (left, right) = line.split(0.5);
        assert_eq!(left.start, ArrayPoint::new([0.0, 0.0]));
        assert_eq!(left.end, ArrayPoint::new([2.0, 0.0]));
        assert_eq!(right.start, ArrayPoint::new([2.0, 0.0]));
        assert_eq!(right.end, ArrayPoint::new([4.0, 0.0]));
    }

    #[test]
    fn line_segment_distance_to_point_after_end() {
        let line = LineSegment {
            start: ArrayPoint::new([0f64, 0f64]),
            end: ArrayPoint::new([1f64, 0f64]),
        };
        let p = ArrayPoint::new([2f64, 0f64]);
        assert_approx_equal(
            line.distance_to_point(p),
            1.0,
            ABSOLUTE_TOLERANCE,
            RELATIVE_TOLERANCE,
        );
    }

    #[test]
    fn small_scale_line_remains_non_degenerate() {
        let scale = 1e-20_f64;
        let line =
            LineSegment::new(ArrayPoint::new([0.0, 0.0]), ArrayPoint::new([scale, 0.0])).unwrap();

        assert_point_approx_equal(
            line.tangent(0.5),
            ArrayPoint::new([1.0, 0.0]),
            ABSOLUTE_TOLERANCE,
            RELATIVE_TOLERANCE,
        );
        assert_approx_equal(
            line.t_at_length(scale * 0.5, 1e-12).unwrap(),
            0.5,
            ABSOLUTE_TOLERANCE,
            RELATIVE_TOLERANCE,
        );
        assert_approx_equal(
            line.distance_to_point(ArrayPoint::new([scale * 0.5, scale])),
            scale,
            1e-32,
            1e-12,
        );
    }

    #[test]
    fn line_root_is_scale_invariant() {
        let line = LineSegment::new(ArrayPoint::new([0.0_f64]), ArrayPoint::new([1.0])).unwrap();
        let roots = line.root(1e-30, -5e-31);

        assert_eq!(roots.len(), 1);
        assert_approx_equal(roots[0], 0.5, ABSOLUTE_TOLERANCE, RELATIVE_TOLERANCE);
    }

    #[test]
    fn line_segment_api_parity() {
        let line =
            LineSegment::new(ArrayPoint::new([0.0_f64, 0.0]), ArrayPoint::new([2.0, 0.0])).unwrap();

        assert_eq!(line.start(), line.eval(0.0));
        assert_eq!(line.end(), line.eval(1.0));

        let reversed = line.reverse();
        assert_eq!(reversed.start(), line.end());
        assert_eq!(reversed.end(), line.start());

        let tangent = line.tangent(0.3);
        assert_approx_equal(tangent[0], 1.0, ABSOLUTE_TOLERANCE, RELATIVE_TOLERANCE);
        assert_approx_equal(tangent[1], 0.0, ABSOLUTE_TOLERANCE, RELATIVE_TOLERANCE);

        let curvature = line.curvature(0.3);
        assert_approx_equal(curvature, 0.0, ABSOLUTE_TOLERANCE, RELATIVE_TOLERANCE);
        assert!(line.normal(0.3).is_none());

        let t = line.t_at_length(1.0, 1e-12).unwrap();
        assert_approx_equal(t, 0.5, ABSOLUTE_TOLERANCE, RELATIVE_TOLERANCE);

        let p = line.point_at_length(1.0, 1e-12).unwrap();
        assert_point_approx_equal(
            p,
            ArrayPoint::new([1.0, 0.0]),
            ABSOLUTE_TOLERANCE,
            RELATIVE_TOLERANCE,
        );
    }
}
