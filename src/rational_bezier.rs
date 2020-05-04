use super::*;
use super::point::Point;
use super::line::LineSegment; 
use super::quadratic_bezier::QuadraticBezier; 
use super::cubic_bezier::CubicBezier;

use num_traits::{float::Float};

/// A Wrapper Class for Bezier curves that makes them rational, allowing 
/// representation of conics by the use of weights for the control points
/// A Rational Bezier is defined by:
///            SUM 0..n ( b(t, i, n) * p[i] * w[i] )
/// B(t,n) = ----------------------------------------
///              SUM 0..n ( b(t, i, n) * w[i] )
/// where b(t,i,n) is the bernstein polynomial of degree n, w[i] is the weights to apply for each control point p[i]
pub enum RationalBezier<P: Point>
{
    Linear(     LineSegment<P>,     [NativeFloat; 2]),
    Quadratic(  QuadraticBezier<P>, [NativeFloat; 3]),
    Cubic(      CubicBezier<P>,     [NativeFloat; 4]),
}


impl<P> RationalBezier<P>
where 
P:  Sub<P, Output = P>
    + Add<P, Output = P>
    + Mul<NativeFloat, Output = P>
    + Point<Scalar = NativeFloat>,
NativeFloat: Sub<NativeFloat, Output = NativeFloat> 
+ Mul<NativeFloat, Output = NativeFloat> {

    /// Creates a quadratic rational Bezier curve out of an existing cubic bezier curve
    pub fn from_linear<F>(&self, curve: LineSegment<P>, weights: [NativeFloat; 2]) -> RationalBezier<P>
    where
    F: Float,
    NativeFloat: Float + Into<F>
    {
        return RationalBezier::Linear(curve, weights)
    }

    /// Creates a quadratic rational Bezier curve out of an existing cubic bezier curve
    pub fn from_quadratic<F>(&self, curve: QuadraticBezier<P>, weights: [NativeFloat; 3]) -> RationalBezier<P>
    where
    F: Float,
    NativeFloat: Float + Into<F>
    {
        return RationalBezier::Quadratic(curve, weights)
    }

    /// Creates a cubic rational Bezier curve out of an existing cubic bezier curve
    pub fn from_cubic<F>(&self, curve: CubicBezier<P>, weights: [NativeFloat; 4]) -> RationalBezier<P>
    where
    F: Float,
    NativeFloat: Float + Into<F>
    {
        return RationalBezier::Cubic(curve, weights)
    }


    /// Evaluates the Bernstein Polynomial b(t,i,n) of the necessary rank (quad. or cubic) adjusted with weights
    /// b_w = sum_i=0_to_n( b(t,i,n) * w_i )
    /// and returns its inverse 1/b_w used as the normalization factor for the rational curve
    fn normalization_factor<F>(&self, t: F) -> F 
    where
    F: Float,
    NativeFloat: Float + Into<F>
    {
        match self {
            Self::Linear(_, weights) => (1.0.into() - t) * weights[0].into() + t * weights[1].into(),

            Self::Quadratic(_, weights) => (1.0.into() - t) * (1.0.into() - t) * weights[0].into()
                                                + 2.0.into() * t * (1.0.into() - t) * weights[1].into()
                                                + t * t * weights[2].into(),

            Self::Cubic(_, weights) => (1.0.into() - t)*(1.0.into() - t) * (1.0.into() - t) * weights[0].into() 
                                            + 3.0.into() * (1.0.into() - t) * (1.0.into() - t) * t * weights[1].into() 
                                            + 3.0.into() * t * t * (1.0.into() - t) * weights[2].into()
                                            + t * t * t * weights[3].into(),
        }
    }


        
    pub fn eval<F>(&self, t: F) -> P 
    where
    F: Float
        + Into<NativeFloat>,
    P: Point<Scalar = NativeFloat>,
    NativeFloat: Sub<F, Output = F> 
        + Add<F, Output = F>
        + Mul<F, Output = F> 
        + Into<F>
    {
        match self {
            Self::Linear(segment, _) => segment.eval(t.into()) * self.normalization_factor(t).into(),
            Self::Quadratic(segment, _) => segment.eval(t.into()) * self.normalization_factor(t).into(),
            Self::Cubic(segment, _) => segment.eval(t.into()) * self.normalization_factor(t).into(),
        }
    }


}