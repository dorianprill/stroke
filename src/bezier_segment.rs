use super::*;
use super::point::Point;
use super::line::LineSegment; 
use super::quadratic_bezier::QuadraticBezier; 
use super::cubic_bezier::CubicBezier;

#[derive(Copy, Clone, Debug, PartialEq)]
pub enum BezierSegment<P: Point<NativeFloat>>
{
    Linear(     LineSegment<P>),
    Quadratic(  QuadraticBezier<P>),
    Cubic(      CubicBezier<P>),
}


impl<P> BezierSegment<P>
where 
P: Point<NativeFloat>,
// Sub<P, Output = P>
//     + Add<P, Output = P>
//     + Mul<NativeFloat, Output = P>,
NativeFloat: Sub<NativeFloat, Output = NativeFloat> 
    + Mul<NativeFloat, Output = NativeFloat> {


    pub fn eval<F>(&self, t: F) -> P
    where 
    F: Float + Copy + Default + Into<NativeFloat> + From<NativeFloat>
    {
        match self {
            BezierSegment::Linear(segment) => segment.eval(t),
            BezierSegment::Quadratic(segment) => segment.eval(t),
            BezierSegment::Cubic(segment) => segment.eval(t),
        }
    }

    pub fn start(&self) -> P {
        match self {
            BezierSegment::Linear(segment) => segment.start,
            BezierSegment::Quadratic(segment) => segment.start,
            BezierSegment::Cubic(segment) => segment.start,
        }
    }

    #[inline]
    pub fn end(&self) -> P {
        match self {
            BezierSegment::Linear(segment) => segment.end,
            BezierSegment::Quadratic(segment) => segment.end,
            BezierSegment::Cubic(segment) => segment.end,
        }
    }

    #[inline]
    pub fn is_linear<F>(&self, tolerance: F) -> bool 
    where 
    F: Float + Into<NativeFloat>
    {
        match self {
            BezierSegment::Linear(..) => true,
            BezierSegment::Quadratic(segment) => segment.is_linear(tolerance.into()),
            BezierSegment::Cubic(segment) => segment.is_linear(tolerance.into()),
        }
    }

    #[inline]
    pub fn baseline(&self) -> LineSegment<P> {
        match self {
            BezierSegment::Linear(segment) => *segment,
            BezierSegment::Quadratic(segment) => segment.baseline(),
            BezierSegment::Cubic(segment) => segment.baseline(),
        }
    }

    /// Split this segment into two sub-segments.
    pub fn split<F>(&self, t: F) -> (BezierSegment<P>, BezierSegment<P>) 
    where 
    F: Float + Into<NativeFloat>
    {
        match self {
            BezierSegment::Linear(segment) => {
                let (a, b) = segment.split(t.into());
                (BezierSegment::Linear(a), BezierSegment::Linear(b))
            }
            BezierSegment::Quadratic(segment) => {
                let (a, b) = segment.split(t.into());
                (BezierSegment::Quadratic(a), BezierSegment::Quadratic(b))
            }
            BezierSegment::Cubic(segment) => {
                let (a, b) = segment.split(t.into());
                (BezierSegment::Cubic(a), BezierSegment::Cubic(b))
            }
        }
    }
}


impl<P> From<LineSegment<P>> for BezierSegment<P>
where
P: Point<NativeFloat> 
{
    fn from(s: LineSegment<P>) -> Self {
        BezierSegment::Linear(s)
    }
}

impl<P> From<QuadraticBezier<P>> for BezierSegment<P> 
where
P: Point<NativeFloat> 
{
    fn from(s: QuadraticBezier<P>) -> Self {
        BezierSegment::Quadratic(s)
    }
}

impl<P> From<CubicBezier<P>> for BezierSegment<P> 
where
P: Point<NativeFloat> 
{
    fn from(s: CubicBezier<P>) -> Self {
        BezierSegment::Cubic(s)
    }
}