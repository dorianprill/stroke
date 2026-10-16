//! Sum type for specialized Bezier segments.

use super::{
    ArcLengthError, BoundingBox, CubicBezier, LineSegment, Point, PointDot, PointIndex, PointNorm,
    QuadraticBezier,
};
use num_traits::{Float, NumCast};

/// Sum type for line/quadratic/cubic Bezier segments.
///
/// Methods that need component access or norms add `PointIndex`/`PointNorm`
/// bounds as required.
#[derive(Copy, Clone, Debug, PartialEq)]
pub enum BezierSegment<P: Point> {
    Linear(LineSegment<P>),
    Quadratic(QuadraticBezier<P>),
    Cubic(CubicBezier<P>),
}

impl<P> BezierSegment<P>
where
    P: Point,
{
    /// Evaluate the segment at `t` in `[0, 1]`.
    pub fn eval(&self, t: P::Scalar) -> P {
        match self {
            BezierSegment::Linear(segment) => segment.eval(t),
            BezierSegment::Quadratic(segment) => segment.eval(t),
            BezierSegment::Cubic(segment) => segment.eval(t),
        }
    }

    /// Return the segment start point.
    pub fn start(&self) -> P {
        match self {
            BezierSegment::Linear(segment) => segment.start,
            BezierSegment::Quadratic(segment) => segment.start,
            BezierSegment::Cubic(segment) => segment.start,
        }
    }

    #[inline]
    /// Return the segment end point.
    pub fn end(&self) -> P {
        match self {
            BezierSegment::Linear(segment) => segment.end,
            BezierSegment::Quadratic(segment) => segment.end,
            BezierSegment::Cubic(segment) => segment.end,
        }
    }

    /// Return a segment with reversed direction.
    pub fn reverse(&self) -> Self {
        match self {
            BezierSegment::Linear(segment) => BezierSegment::Linear(segment.reverse()),
            BezierSegment::Quadratic(segment) => BezierSegment::Quadratic(segment.reverse()),
            BezierSegment::Cubic(segment) => BezierSegment::Cubic(segment.reverse()),
        }
    }

    #[inline]
    /// Return true if the segment is linear within `tolerance`.
    pub fn is_linear<F>(&self, tolerance: F) -> bool
    where
        P: PointNorm,
        F: Float + NumCast,
    {
        let tolerance = <P::Scalar as NumCast>::from(tolerance).unwrap();
        match self {
            BezierSegment::Linear(..) => true,
            BezierSegment::Quadratic(segment) => segment.is_linear(tolerance),
            BezierSegment::Cubic(segment) => segment.is_linear(tolerance),
        }
    }

    #[inline]
    /// Return the baseline line segment between start and end.
    pub fn baseline(&self) -> LineSegment<P> {
        match self {
            BezierSegment::Linear(segment) => *segment,
            BezierSegment::Quadratic(segment) => segment.baseline(),
            BezierSegment::Cubic(segment) => segment.baseline(),
        }
    }

    #[inline]
    /// Return the bounding box across all axes.
    pub fn bounding_box(&self) -> BoundingBox<P::Scalar, { P::DIM }>
    where
        P: PointIndex,
        [P::Scalar; 1]: tinyvec::Array<Item = P::Scalar>,
        [P::Scalar; 2]: tinyvec::Array<Item = P::Scalar>,
        [P::Scalar; 3]: tinyvec::Array<Item = P::Scalar>,
        [P::Scalar; 4]: tinyvec::Array<Item = P::Scalar>,
    {
        match self {
            BezierSegment::Linear(segment) => segment.bounding_box(),
            BezierSegment::Quadratic(segment) => segment.bounding_box(),
            BezierSegment::Cubic(segment) => segment.bounding_box(),
        }
    }

    /// Return a conservative box around the segment's control points.
    pub fn control_box(&self) -> BoundingBox<P::Scalar, { P::DIM }>
    where
        P: PointIndex,
    {
        match self {
            BezierSegment::Linear(segment) => segment.control_box(),
            BezierSegment::Quadratic(segment) => segment.control_box(),
            BezierSegment::Cubic(segment) => segment.control_box(),
        }
    }

    /// Return the arc length using adaptive Gauss-Legendre quadrature.
    pub fn arclen(&self, accuracy: P::Scalar) -> Result<P::Scalar, ArcLengthError>
    where
        P: PointNorm,
    {
        match self {
            BezierSegment::Linear(segment) => segment.arclen(accuracy),
            BezierSegment::Quadratic(segment) => segment.arclen(accuracy),
            BezierSegment::Cubic(segment) => segment.arclen(accuracy),
        }
    }

    /// Integrate length over `start..end`, where `0 <= start <= end <= 1`.
    pub fn arclen_between(
        &self,
        start: P::Scalar,
        end: P::Scalar,
        accuracy: P::Scalar,
    ) -> Result<P::Scalar, ArcLengthError>
    where
        P: PointNorm,
    {
        match self {
            BezierSegment::Linear(segment) => segment.arclen_between(start, end, accuracy),
            BezierSegment::Quadratic(segment) => segment.arclen_between(start, end, accuracy),
            BezierSegment::Cubic(segment) => segment.arclen_between(start, end, accuracy),
        }
    }

    /// Approximate the arc length using `segments` line segments.
    pub fn arclen_polyline(&self, segments: usize) -> P::Scalar
    where
        P: PointNorm,
    {
        match self {
            BezierSegment::Linear(segment) => segment.arclen_polyline(segments),
            BezierSegment::Quadratic(segment) => segment.arclen_polyline(segments),
            BezierSegment::Cubic(segment) => segment.arclen_polyline(segments),
        }
    }

    /// Split this segment into two sub-segments at `t`.
    pub fn split<F>(&self, t: F) -> (BezierSegment<P>, BezierSegment<P>)
    where
        F: Float + NumCast,
    {
        let t = <P::Scalar as NumCast>::from(t).unwrap();
        match self {
            BezierSegment::Linear(segment) => {
                let (a, b) = segment.split(t);
                (BezierSegment::Linear(a), BezierSegment::Linear(b))
            }
            BezierSegment::Quadratic(segment) => {
                let (a, b) = segment.split(t);
                (BezierSegment::Quadratic(a), BezierSegment::Quadratic(b))
            }
            BezierSegment::Cubic(segment) => {
                let (a, b) = segment.split(t);
                (BezierSegment::Cubic(a), BezierSegment::Cubic(b))
            }
        }
    }

    /// Evaluate the derivative with respect to the unit-domain parameter.
    pub fn derivative_at(&self, t: P::Scalar) -> P {
        match self {
            Self::Linear(curve) => curve.derivative_at(t),
            Self::Quadratic(curve) => curve.derivative_at(t),
            Self::Cubic(curve) => curve.derivative_at(t),
        }
    }

    /// Return derivative magnitude per unit of the unit-domain parameter.
    pub fn speed(&self, t: P::Scalar) -> P::Scalar
    where
        P: PointNorm,
    {
        self.derivative_at(t).norm()
    }

    /// Extract a finite ordered range in `[0, 1]`, reparameterized to `[0, 1]`.
    /// Equal endpoints yield a constant segment.
    pub fn subcurve(
        &self,
        range: core::ops::Range<P::Scalar>,
    ) -> Result<Self, crate::CurveParameterError> {
        match self {
            Self::Linear(curve) => curve.subcurve(range).map(Self::Linear),
            Self::Quadratic(curve) => curve.subcurve(range).map(Self::Quadratic),
            Self::Cubic(curve) => curve.subcurve(range).map(Self::Cubic),
        }
    }

    /// Split at a distance from the start, clamping finite distances to endpoints.
    pub fn split_at_length(
        &self,
        length: P::Scalar,
        accuracy: P::Scalar,
    ) -> Result<(Self, Self), ArcLengthError>
    where
        P: PointNorm,
    {
        match self {
            Self::Linear(curve) => curve
                .split_at_length(length, accuracy)
                .map(|(a, b)| (Self::Linear(a), Self::Linear(b))),
            Self::Quadratic(curve) => curve
                .split_at_length(length, accuracy)
                .map(|(a, b)| (Self::Quadratic(a), Self::Quadratic(b))),
            Self::Cubic(curve) => curve
                .split_at_length(length, accuracy)
                .map(|(a, b)| (Self::Cubic(a), Self::Cubic(b))),
        }
    }

    /// Map control points and validate the result without modifying this segment.
    /// Only affine maps geometrically transform the whole curve in the same way.
    pub fn try_map_control_points(
        &self,
        map: impl FnMut(P) -> P,
    ) -> Result<Self, crate::CurveConstructionError> {
        match self {
            Self::Linear(curve) => curve.try_map_control_points(map).map(Self::Linear),
            Self::Quadratic(curve) => curve.try_map_control_points(map).map(Self::Quadratic),
            Self::Cubic(curve) => curve.try_map_control_points(map).map(Self::Cubic),
        }
    }

    /// Return the unit tangent direction at `t`.
    pub fn tangent(&self, t: P::Scalar) -> P
    where
        P: PointNorm,
    {
        match self {
            BezierSegment::Linear(segment) => segment.tangent(t),
            BezierSegment::Quadratic(segment) => segment.tangent(t),
            BezierSegment::Cubic(segment) => segment.tangent(t),
        }
    }

    /// Return the curvature magnitude at `t`.
    pub fn curvature(&self, t: P::Scalar) -> P::Scalar
    where
        P: PointNorm + PointDot,
    {
        match self {
            BezierSegment::Linear(segment) => segment.curvature(t),
            BezierSegment::Quadratic(segment) => segment.curvature(t),
            BezierSegment::Cubic(segment) => segment.curvature(t),
        }
    }

    /// Return the principal normal direction at `t`.
    ///
    /// Returns `None` if the velocity is zero or curvature is undefined.
    pub fn normal(&self, t: P::Scalar) -> Option<P>
    where
        P: PointNorm + PointDot,
    {
        match self {
            BezierSegment::Linear(segment) => segment.normal(t),
            BezierSegment::Quadratic(segment) => segment.normal(t),
            BezierSegment::Cubic(segment) => segment.normal(t),
        }
    }

    /// Return parameter `t` at arc length `s` using adaptive quadrature.
    pub fn t_at_length(
        &self,
        s: P::Scalar,
        accuracy: P::Scalar,
    ) -> Result<P::Scalar, ArcLengthError>
    where
        P: PointNorm,
    {
        match self {
            BezierSegment::Linear(segment) => segment.t_at_length(s, accuracy),
            BezierSegment::Quadratic(segment) => segment.t_at_length(s, accuracy),
            BezierSegment::Cubic(segment) => segment.t_at_length(s, accuracy),
        }
    }

    /// Evaluate the point at arc length `s` using adaptive quadrature.
    pub fn point_at_length(&self, s: P::Scalar, accuracy: P::Scalar) -> Result<P, ArcLengthError>
    where
        P: PointNorm,
    {
        Ok(self.eval(self.t_at_length(s, accuracy)?))
    }
}

impl<P> From<LineSegment<P>> for BezierSegment<P>
where
    P: Point,
{
    fn from(s: LineSegment<P>) -> Self {
        BezierSegment::Linear(s)
    }
}

impl<P> From<QuadraticBezier<P>> for BezierSegment<P>
where
    P: Point,
{
    fn from(s: QuadraticBezier<P>) -> Self {
        BezierSegment::Quadratic(s)
    }
}

impl<P> From<CubicBezier<P>> for BezierSegment<P>
where
    P: Point,
{
    fn from(s: CubicBezier<P>) -> Self {
        BezierSegment::Cubic(s)
    }
}

impl<P> Default for BezierSegment<P>
where
    P: Point + Default,
{
    fn default() -> Self {
        BezierSegment::Linear(LineSegment {
            start: P::default(),
            end: P::default(),
        })
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
    fn bezier_segment_api_parity() {
        let line = LineSegment::new(PointN::new([0f64, 0f64]), PointN::new([2f64, 0f64])).unwrap();
        let segment = BezierSegment::from(line);

        assert_eq!(segment.start(), segment.eval(0.0));
        assert_eq!(segment.end(), segment.eval(1.0));

        let reversed = segment.reverse();
        let p0 = segment.eval(0.25);
        let p1 = reversed.eval(0.75);
        assert_point_approx_equal(p0, p1, ABSOLUTE_TOLERANCE, RELATIVE_TOLERANCE);

        let tangent = segment.tangent(0.3);
        assert_approx_equal(tangent[0], 1.0, ABSOLUTE_TOLERANCE, RELATIVE_TOLERANCE);
        assert_approx_equal(tangent[1], 0.0, ABSOLUTE_TOLERANCE, RELATIVE_TOLERANCE);

        let curvature = segment.curvature(0.3);
        assert_approx_equal(curvature, 0.0, ABSOLUTE_TOLERANCE, RELATIVE_TOLERANCE);
        assert!(segment.normal(0.3).is_none());

        let t = segment.t_at_length(1.0, 1e-12).unwrap();
        assert_approx_equal(t, 0.5, ABSOLUTE_TOLERANCE, RELATIVE_TOLERANCE);

        let p = segment.point_at_length(1.0, 1e-12).unwrap();
        assert_point_approx_equal(
            p,
            PointN::new([1.0, 0.0]),
            ABSOLUTE_TOLERANCE,
            RELATIVE_TOLERANCE,
        );
    }
}
