//! Sum type for specialized Bezier segments.

use num_traits::{Float, NumCast};
use super::{CubicBezier, LineSegment, Point, PointIndex, PointNorm, QuadraticBezier};

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
