use super::*;
#[allow(unused_imports)]
use super::point2::{Point2, Coordinate, Distance};
use super::line::Line; 
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
    pub fn derivative<F>(&self) -> Line<P>
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
        return Line{
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


    /// Solve for the roots of the bezier curve
    /// Will return an array of roots in the order: [x1, y1, x2, y2] 
    /// All roots not positive will return NaN an need to be checked
    /// All roots not in [0,1] are also not meaningful in this context and can be discarded
    pub fn roots(&self) -> [NativeFloat; 4] 
    where
    P:  Sub<P, Output = P>
        + Add<P, Output = P>
        + Mul<NativeFloat, Output = P>,
    NativeFloat: Sub<NativeFloat, Output = NativeFloat> 
        + Mul<NativeFloat, Output = NativeFloat>
        {
        // by substituion we get the coefficients for the quadratic formula
        // r1,2 = ( -b +/- sqrt(b^2 - 4ac) ) / 2a 
        let a: P = self.start - (self.ctrl * 2. as NativeFloat) + self.end;
        let b: P = (self.ctrl - self.start) * 2. as NativeFloat;
        let c: P = self.start;
        let rx1: NativeFloat =  (b.x() * (-1.) + (b.x() * b.x() - a.x() * c.x() * 4.).sqrt() ) / (a.x() * 2.);
        let rx2: NativeFloat =  (b.x() * (-1.) - (b.x() * b.x() - a.x() * c.x() * 4.).sqrt() ) / (a.x() * 2.);

        let ry1: NativeFloat =  (b.y() * (-1.) + (b.y() * b.y() - a.y() * c.y() * 4.).sqrt() ) / (a.y() * 2.);
        let ry2: NativeFloat =  (b.y() * (-1.) - (b.y() * b.y() - a.y() * c.y() * 4.).sqrt() ) / (a.y() * 2.);

        return [rx1, ry1, rx2, ry2];
        
    }
}

