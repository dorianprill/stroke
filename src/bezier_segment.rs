//! Sum type for specialized Bezier segments.

use num_traits::{Float, NumCast};
use super::{CubicBezier, LineSegment, Point, PointDot, PointIndex, PointNorm, QuadraticBezier};

const DEFAULT_LENGTH_STEPS: usize = 64;

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
    pub fn bounding_box(&self) -> [(P::Scalar, P::Scalar); P::DIM]
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

    /// Approximate parameter `t` at arc length `s`.
    pub fn t_at_length_approx(&self, s: P::Scalar, nsteps: usize) -> P::Scalar
    where
        P: PointNorm,
    {
        match self {
            BezierSegment::Linear(segment) => segment.t_at_length_approx(s, nsteps),
            BezierSegment::Quadratic(segment) => segment.t_at_length_approx(s, nsteps),
            BezierSegment::Cubic(segment) => segment.t_at_length_approx(s, nsteps),
        }
    }

    /// Approximate parameter `t` at arc length `s` using a default resolution.
    pub fn t_at_length(&self, s: P::Scalar) -> P::Scalar
    where
        P: PointNorm,
    {
        self.t_at_length_approx(s, DEFAULT_LENGTH_STEPS)
    }

    /// Evaluate the point at arc length `s`.
    pub fn point_at_length_approx(&self, s: P::Scalar, nsteps: usize) -> P
    where
        P: PointNorm,
    {
        self.eval(self.t_at_length_approx(s, nsteps))
    }

    /// Evaluate the point at arc length `s` using a default resolution.
    pub fn point_at_length(&self, s: P::Scalar) -> P
    where
        P: PointNorm,
    {
        self.point_at_length_approx(s, DEFAULT_LENGTH_STEPS)
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
        BezierSegment::Linear(LineSegment::new(P::default(), P::default()))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{PointN, PointNorm, EPSILON};

    #[test]
    fn bezier_segment_api_parity() {
        let line = LineSegment::new(PointN::new([0f64, 0f64]), PointN::new([2f64, 0f64]));
        let segment = BezierSegment::from(line);

        assert_eq!(segment.start(), segment.eval(0.0));
        assert_eq!(segment.end(), segment.eval(1.0));

        let reversed = segment.reverse();
        let p0 = segment.eval(0.25);
        let p1 = reversed.eval(0.75);
        assert!((p0 - p1).squared_norm() < EPSILON);

        let tangent = segment.tangent(0.3);
        assert!((tangent[0] - 1.0).abs() < EPSILON);
        assert!(tangent[1].abs() < EPSILON);

        let curvature = segment.curvature(0.3);
        assert!(curvature.abs() < EPSILON);
        assert!(segment.normal(0.3).is_none());

        let t = segment.t_at_length(1.0);
        assert!((t - 0.5).abs() < EPSILON);

        let p = segment.point_at_length(1.0);
        assert!((p - PointN::new([1.0, 0.0])).squared_norm() < EPSILON);
    }
}
