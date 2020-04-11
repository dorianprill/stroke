use super::*;
#[allow(unused_imports)]
use super::point2::{Point2, Coordinate, Distance};
#[allow(unused_imports)]
use super::line::{Line, LineSegment}; 
#[allow(unused_imports)]
use super::cubic_bezier::CubicBezier;

#[derive(Copy, Clone, Debug, PartialEq)]
pub struct QuadraticBezier<P>
{
    pub(crate) start:  P,
    pub(crate) ctrl:   P,
    pub(crate) end:    P,
}

impl<P> QuadraticBezier<P>
where
P: Add + Sub + Copy
    + Add<P, Output = P>
    + Sub<P, Output = P>
    + Mul<NativeFloat, Output = P>
    + Distance<ScalarDist = NativeFloat> 
    + Coordinate<Coordinate = NativeFloat>,
NativeFloat: Sub<NativeFloat, Output = NativeFloat> 
    + Mul<NativeFloat, Output = NativeFloat> {

    pub fn new(start: P, ctrl: P, end: P) -> Self {
        QuadraticBezier { 
            start, 
            ctrl, 
            end 
        }
    }

    pub fn eval(&self, t: NativeFloat) -> P {
        let t2:     NativeFloat = t * t;
        let one_t:  NativeFloat = 1.0 as NativeFloat - t;
        let one_t2: NativeFloat = one_t * one_t;

        self.start * one_t2
            + self.ctrl * 2.0 as NativeFloat * one_t * t
            + self.end * t2
    }

        /// Evaluate a CubicBezier curve at t using the numerically stable De Casteljau algorithm
    pub fn eval_casteljau(&self, t: NativeFloat) -> P 
    {
        // unrolled de casteljau algorithm
        // _1ab is the first iteration from first (a) to second (b) control point and so on
        let ctrl_1ab = self.start + (self.ctrl - self.start) * t;
        let ctrl_1bc   = self.ctrl + (self.end - self.ctrl) * t;
        // second iteration, final point on the curve
        let ctrl_2ab  = ctrl_1ab + (ctrl_1bc - ctrl_1ab) * t;
    
        return ctrl_2ab
    }

    /// Sample the x coordinate of the curve at t (expecting t between 0 and 1).
    pub fn x<F>(&self, t: F) -> F 
    where
    F : Float,
    P:  Sub<P, Output = P>
        + Add<P, Output = P>
        + Mul<F, Output = P>,
    NativeFloat: Sub<F, Output = F> 
        + Add<F, Output = F>
        + Mul<F, Output = F> 
        + Into<F>
    {
        let t2 = t * t;
        let one_t = 1 as NativeFloat - t;
        let one_t2 = one_t * one_t;

        self.start.x() * one_t2 + self.ctrl.x() * 2.0.into() * one_t * t + self.end.x() * t2
    }

    /// Sample the y coordinate of the curve at t (expecting t between 0 and 1).
    pub fn y<F>(&self, t: F) -> F 
    where
    F : Float,
    P:  Sub<P, Output = P>
        + Add<P, Output = P>
        + Mul<F, Output = P>,
    NativeFloat: Sub<F, Output = F> 
        + Add<F, Output = F>
        + Mul<F, Output = F> 
        + Into<F>
    {
        let t2 = t * t;
        let one_t = 1 as NativeFloat - t;
        let one_t2 = one_t * one_t;

        self.start.y() * one_t2 + self.ctrl.y() * 2.0.into() * one_t * t + self.end.y() * t2
    }

    /// Return the derivative function.
    /// The derivative is also a bezier curve but of degree n-1 - In the case of quadratic it is just a line.
    /// Since it returns the derivative function, eval() needs to be called separately
    pub fn derivative<F>(&self) -> LineSegment<P>
    where
    F: Float,
    P:  Sub<P, Output = P>
        + Add<P, Output = P>
        + Mul<F, Output = P>,
    NativeFloat: Sub<F, Output = F> 
        + Add<F, Output = F>
        + Mul<F, Output = F>
        + Into<F>
    {
        return LineSegment{
            start: (self.ctrl - self.start) * 2.0.into(),
            end:   (self.end - self.ctrl)   * 2.0.into()
        }
    }


    /// Sample the x coordinate of the curve's derivative at t (expecting t between 0 and 1).
    /// Convenience function for .derivative().eval(t).x()
    pub fn dx<F>(&self, t: F) -> F
    where
    F: Float,
    P:  Sub<P, Output = P>
        + Add<P, Output = P>
        + Mul<F, Output = P>,
    NativeFloat: Sub<F, Output = F> 
        + Mul<F, Output = F>
        + Float
        + Into<F>
    {
        let c0: F = t * 2.0.into() - 2.0.into();
        let c1: F =  2.0.into() - 4.0.into() * t;
        let c2: F = 2.0.into() * t;
        return self.start.x().into() * c0 + self.ctrl.x().into() * c1 + self.end.x().into() * c2;
    }


    /// Sample the y coordinate of the curve's derivative at t (for t between 0 and 1).
    pub fn dy<F>(&self, t: F) -> F
    where
    F: Float,
    P:  Sub<P, Output = P>
        + Add<P, Output = P>
        + Mul<F, Output = P>,
    NativeFloat: Sub<F, Output = F> 
        + Mul<F, Output = F>
        + Into<F>
    {
        let c0: F = t * 2.0.into() - 2.0.into();
        let c1: F =  2.0.into() - 4.0.into() * t;
        let c2: F = 2.0.into() * t;
        return self.start.y().into() * c0 + self.ctrl.y().into() * c1 + self.end.y().into() * c2;
    }

    /// Calculates the curvature of the curve at point t
    /// The curvature is the inverse of the radius of the tangential circle at t: k=1/r
    pub fn curvature<F>(&self, t: F) -> F
    where
    F: Float,
    P:  Sub<P, Output = P>
        + Add<P, Output = P>
        + Mul<F, Output = P>,
    NativeFloat: Sub<F, Output = F> 
        + Add<F, Output = F>
        + Mul<F, Output = F>
        + Into<F>
    {
        let d = self.derivative();
        let dd = d.derivative();
        let dx = d.x(t);
        let dy = d.y(t);
        let ddx = dd.x();
        let ddy = dd.y();
        let numerator = dx * ddy.into() - ddx * dy;
        let denominator = (dx*dx + dy*dy).powf(1.5.into());
        return numerator / denominator
    }

    /// Calculates the radius of the tangential circle at t
    /// It is the inverse of the curvature at t: r=1/k
    pub fn radius<F>(&self, t: F) -> F
    where
    F: Float,
    P:  Sub<P, Output = P>
        + Add<P, Output = P>
        + Mul<F, Output = P>,
    NativeFloat: Sub<F, Output = F> 
        + Add<F, Output = F>
        + Mul<F, Output = F>
        + Into<F>
    {
        return 1.0.into() / self.curvature(t)
    }


    /// Solve for the roots of the polynomial at^2 + bt + c
    /// Returns an ArrayVec of roots in the order
    /// needs to be called for x and y components separately
    pub(crate) fn real_roots<F>(&self, a: F, b: F, c: F) -> ArrayVec<[F; 2]>
    where
    F: Float,
    P:  Sub<P, Output = P>
        + Add<P, Output = P>
        + Mul<NativeFloat, Output = P>,
        NativeFloat: Sub<F, Output = F> 
        + Add<F, Output = F>
        + Mul<F, Output = F>
        + Float
        + Into<F>
    {

        let mut result = ArrayVec::new();
        let epsilon = 1e-5.into();

        // check if can be handled below cubic order
        if a.abs() < epsilon {
            if b.abs() < epsilon {
                // no solutions
                return result;
            }
            // is linear equation
            result.push(-c / b);
            return result;
        }
        // is quadratic equation
        let delta = b * b - 4.0.into() * a * c;
        if delta > 0.0.into() {
            let sqrt_delta = delta.sqrt();
            result.push((-b - sqrt_delta) / (2.0.into() * a));
            result.push((-b + sqrt_delta) / (2.0.into() * a));
        } else if delta.abs() < epsilon {
            result.push(-b / (2.0.into() * a));
        }
        return result;
    }

}