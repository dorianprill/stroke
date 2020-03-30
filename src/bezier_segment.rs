use super::*;
use super::point2::{Point2, Coordinate, Distance};
use super::line::Line; 
use super::quadratic_bezier::QuadraticBezier; 
use super::cubic_bezier::CubicBezier;

#[derive(Copy, Clone, Debug, PartialEq)]
pub enum BezierSegment<P> {
    Linear(     Line<P>),
    Quadratic(  QuadraticBezier<P>),
    Cubic(      CubicBezier<P>),
}


impl<P> BezierSegment<P>
where 
P:  Sub<P, Output = P>
    + Add<P, Output = P>
    + Mul<NativeFloat, Output = P>
    + Distance<ScalarDist = NativeFloat>
    + Coordinate< Coordinate = NativeFloat>
    + Copy,
NativeFloat: Sub<NativeFloat, Output = NativeFloat> 
+ Mul<NativeFloat, Output = NativeFloat> {
    pub fn eval(&self, t: NativeFloat) -> P {
        match self {
            BezierSegment::Linear(segment) => segment.eval(t),
            BezierSegment::Quadratic(segment) => segment.eval(t),
            BezierSegment::Cubic(segment) => segment.eval(t),
        }
    }
}