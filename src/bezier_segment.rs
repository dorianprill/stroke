use super::*;
use super::point::Point;
use super::line::LineSegment; 
use super::quadratic_bezier::QuadraticBezier; 
use super::cubic_bezier::CubicBezier;

#[derive(Copy, Clone, Debug, PartialEq)]
pub enum BezierSegment<P: Point>
{
    Linear(     LineSegment<P>),
    Quadratic(  QuadraticBezier<P>),
    Cubic(      CubicBezier<P>),
}


impl<P> BezierSegment<P>
where 
P:  Sub<P, Output = P>
    + Add<P, Output = P>
    + Mul<NativeFloat, Output = P>
    + Point<Scalar = NativeFloat>,
NativeFloat: Sub<NativeFloat, Output = NativeFloat> 
+ Mul<NativeFloat, Output = NativeFloat> {
    pub fn eval<F>(&self, t: F) -> P
    where 
    F: Float 
    + Into<NativeFloat>,
    NativeFloat: Sub<F, Output = F> 
    + Add<F, Output = F>
    + Mul<F, Output = F> 
    + Into<F> {
        match self {
            BezierSegment::Linear(segment) => segment.eval(t.into()),
            BezierSegment::Quadratic(segment) => segment.eval(t.into()),
            BezierSegment::Cubic(segment) => segment.eval(t.into()),
        }
    }
}