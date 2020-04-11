use super::*;
#[allow(unused_imports)]
use super::point2::{Point2, Coordinate, Distance};
#[allow(unused_imports)]
use super::line::{Line, LineSegment, LineEquation};
use super::quadratic_bezier::QuadraticBezier;

use num_traits::{float::Float};

/// Describes a cubic bezier curve through its control points.
/// To actually render the curve, the eval() function must be called at desired points t
/// t is an interpolation parameter between [0,1] (0% and 100%)
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct CubicBezier<P>
{
    start:  P,
    ctrl1:  P,
    ctrl2:  P,
    end:    P,
}

#[allow(dead_code)]
impl<P> CubicBezier<P> 
where 
    P: Add 
        + Sub 
        + Copy 
        + Mul<NativeFloat, Output = P>
        + Distance<ScalarDist = NativeFloat> 
        + Coordinate<Coordinate = NativeFloat>,
{

    pub fn new(start: P, ctrl1: P, ctrl2: P,  end: P) -> Self 
    {
        CubicBezier { 
            start, 
            ctrl1, 
            ctrl2, 
            end 
        }
    }

    /// Evaluate a CubicBezier curve at t by direct evaluation of the polynomial (not numerically stable)
    pub fn eval<F>(&self, t: F) -> P 
    where 
    F: Float,
    P: Add<P, Output = P>
        + Sub<P, Output = P>
        + Mul<F, Output = P>,
    NativeFloat: Sub<F, Output = F> 
        + Mul<F, Output = F>
    {
        return self.start * ((1.0-t) * (1.0-t) * (1.0-t))
                + self.ctrl1 * (3.0 * t * (1.0-t) * (1.0-t))
                + self.ctrl2 * (3.0 * t * t * (1.0-t))
                + self.end * (t * t * t);
    }

    /// Evaluate a CubicBezier curve at t using the numerically stable De Casteljau algorithm
    pub fn eval_casteljau<F>(&self, t: F) -> P 
    where 
    F: Float,
    P: Add<P, Output = P>
        + Sub<P, Output = P>
        + Mul<F, Output = P>,
    NativeFloat: Sub<F, Output = F> 
        + Mul<F, Output = F>
    {
        // unrolled de casteljau algorithm
        // _1ab is the first iteration from first (a) to second (b) control point and so on
        let ctrl_1ab = self.start + (self.ctrl1 - self.start) * t;
        let ctrl_1bc = self.ctrl1 + (self.ctrl2 - self.ctrl1) * t;
        let ctrl_1cd = self.ctrl2 + (self.end - self.ctrl2)   * t;
        // second iteration
        let ctrl_2ab  = ctrl_1ab + (ctrl_1bc - ctrl_1ab) * t;
        let ctrl_2bc  = ctrl_1bc + (ctrl_1cd - ctrl_1bc) * t;
        // third iteration, final point on the curve
        let ctrl_3ab = ctrl_2ab + (ctrl_2bc - ctrl_2ab) * t;

        return ctrl_3ab
    }


    pub fn x<F>(&self, t: F) -> F
    where
    F : Float,
    P: Mul<NativeFloat, Output = P>,
    NativeFloat: Sub<F, Output = F> 
        + Mul<F, Output = F> 
        + Into<F>
    {
        let t2 = t * t;
        let t3  = t2 * t;
        let one_t  = (1.0 as NativeFloat) - t;
        let one_t2 = one_t * one_t;
        let one_t3 = one_t2 * one_t;

        self.start.x().into() * one_t3
            + self.ctrl1.x().into() * 3.0.into() * one_t2 * t
            + self.ctrl2.x().into() * 3.0.into() * one_t * t2
            + self.end.x() * t3
    }

    pub fn y<F>(&self, t: F) -> F
    where
    F : Float,
    P: Mul<NativeFloat, Output = P>,
    NativeFloat: Sub<F, Output = F> 
        + Mul<F, Output = F> 
        + Into<F>
    {
        let t2 = t * t;
        let t3  = t2 * t;
        let one_t  = (1.0 as NativeFloat) - t;
        let one_t2 = one_t * one_t;
        let one_t3 = one_t2 * one_t;

        self.start.y().into() * one_t3
            + self.ctrl1.y().into() * 3.0.into() * one_t2 * t
            + self.ctrl2.y().into() * 3.0.into() * one_t * t2
            + self.end.y() * t3
    }

    /// Approximates the arc length of the curve by flattening it with straight line segments.
    /// This works quite well, at ~32 segments it should already provide an error < 0.5
    /// Remember arclen also works by linear approximation, not the integral, so we have to accept error!
    /// This approximation is unfeasable if desired accuracy is greater than 2 decimal places
    pub fn arclen<F>(&self, nsteps: usize) -> F
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
        let stepsize: NativeFloat = 1.0/(nsteps as NativeFloat);
        let mut arclen: NativeFloat = 0.0;
        for t in 1..nsteps {
            let t = t as NativeFloat * 1.0/(nsteps as NativeFloat);
            let p1 = self.eval_casteljau(t);
            let p2 = self.eval_casteljau(t+stepsize);

            arclen = arclen + p1.distance(p2);
        
        }
        return arclen.into()
    }


    pub fn split<F>(&self, t: F) -> (Self, Self)
    where
    F: Float,
    P:  Sub<P, Output = P>
        + Add<P, Output = P>
        + Mul<F, Output = P>,
    NativeFloat: Sub<F, Output = F> 
        + Mul<F, Output = F>,
    {
       // unrolled de casteljau algorithm
        // _1ab is the first iteration from first (a) to second (b) control point and so on
        let ctrl_1ab = self.start + (self.ctrl1 - self.start) * t;
        let ctrl_1bc   = self.ctrl1 + (self.ctrl2 - self.ctrl1) * t;
        let ctrl_1cd   = self.ctrl2 + (self.end - self.ctrl2) * t;
        // second iteration
        let ctrl_2ab  = ctrl_1ab + (ctrl_1bc - ctrl_1ab) * t;
        let ctrl_2bc  = ctrl_1bc + (ctrl_1cd - ctrl_1bc) * t;
        // third iteration, final point on the curve
        let ctrl_3ab = ctrl_2ab + (ctrl_2bc - ctrl_2ab) * t;

        return (
            CubicBezier {
                start: self.start,
                ctrl1: ctrl_1ab,
                ctrl2: ctrl_2ab,
                end: ctrl_3ab,
            },
            CubicBezier {
                start: ctrl_3ab,
                ctrl1: ctrl_2bc,
                ctrl2: ctrl_1cd,
                end: self.end,
            },
        );
    }

    /// Return the derivative function.
    /// The derivative is also a bezier curve but of degree n-1 (cubic->quadratic)
    /// ince it returns the derivative function, eval() needs to be called separately
    pub fn derivative<F>(&self) -> QuadraticBezier<P>
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
        return QuadraticBezier{
            start: (self.ctrl1 - self.start) * 3.0.into(),
            ctrl:  (self.ctrl2 - self.ctrl1) * 3.0.into(),
            end:   (self.end - self.ctrl2)   * 3.0.into()
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
        let t2 = t*t;
        let c0 = -3.0.into() * t + 6.0.into() * t - 3.0.into();
        let c1 = 9.0.into() * t2 - 12.0.into() * t + 3.0.into();
        let c2 = -9.0.into() * t2 + 6.0.into() * t;
        let c3 = 3.0.into() * t2;
        return self.start.x().into() * c0 
                + self.ctrl1.x().into() * c1 
                + self.ctrl2.x().into() * c2 
                + self.end.x().into() * c3
    }



    /// Sample the y coordinate of the curve's derivative at t (expecting t between 0 and 1).
    /// Convenience function for .derivative().eval(t).x()
    pub fn dy<F>(&self, t: F) -> F 
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
        let t2 = t*t;
        let c0 = -3.0.into() * t + 6.0.into() * t - 3.0.into();
        let c1 = 9.0.into() * t2 - 12.0.into() * t + 3.0.into();
        let c2 = -9.0.into() * t2 + 6.0.into() * t;
        let c3 = 3.0.into() * t2;
        return self.start.y().into() * c0 
                + self.ctrl1.y().into() * c1 
                + self.ctrl2.y().into() * c2 
                + self.end.y().into() * c3
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
        + Float
        + Into<F>
    {
        let d = self.derivative();
        let dd = d.derivative();
        let dx = d.x(t);
        let dy = d.y(t);
        let ddx = dd.x(t);
        let ddy = dd.y(t);
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
        + Float
        + Into<F>
    {
        return 1.0.into() / self.curvature(t)
    }


    pub fn baseline(&self) -> LineSegment<P> {
        LineSegment {
            start: self.start,
            end: self.end,
        }
    }


    fn non_point_is_linear<F>(&self, tolerance: F) -> bool 
    where
    F: Float,
    P:  Sub<P, Output = P>
        + Add<P, Output = P>
        + Mul<F, Output = P>,
    NativeFloat: Sub<F, Output = F> 
        + Add<F, Output = F>
        + Mul<F, Output = F>
        + Float
        + Into<F>
    {
        let line = self.baseline().to_line().equation();
        line.distance_to_point(self.ctrl1) <= tolerance
            && line.distance_to_point(self.ctrl2) <= tolerance
    }

    pub(crate) fn is_a_point<F>(&self, tolerance: F) -> bool 
    where
    F: Float,
    P:  Sub<P, Output = P>
        + Add<P, Output = P>
        + Mul<F, Output = P>,
    NativeFloat: Sub<F, Output = F> 
        + Add<F, Output = F>
        + Mul<F, Output = F>
        + Float
        + Into<F>
    {
        let tolerance_squared = tolerance * tolerance;
        // Use <= so that tolerance can be zero.
        self.start.distance(self.end).powi(2).into() <= tolerance_squared
            && self.start.distance(self.ctrl1).powi(2).into() <= tolerance_squared
            && self.end.distance(self.ctrl2).powi(2).into() <= tolerance_squared
    }

    /// Compute the real roots of the cubic bezier function
    /// of the form a*t^3 + b*t^2 + c*t + d
    /// using cardano's algorithm
    fn real_roots<F>(&self, a: F, b: F, c: F, d: F) -> ArrayVec<[F; 3]>
    where
    F: Float,
    P:  Sub<P, Output = P>
        + Add<P, Output = P>
        + Mul<F, Output = P>,
    NativeFloat: Sub<F, Output = F> 
        + Add<F, Output = F>
        + Mul<F, Output = F>
        + Float
        + Into<F>
    {
        let mut result = ArrayVec::new();

        let epsilon = 1e-5.into();
        let pi = 3.141592.into();

        // check if can be handled below cubic order
        if a.abs() < epsilon {
            if b.abs() < epsilon {
                if c.abs() < epsilon {
                    // no solutions
                    return result;
                }
                // is linear equation
                result.push(-d / c);
                return result;
            }
            // is quadratic equation
            let delta = c * c - 4.0.into() * b * d;
            if delta > 0.0.into() {
                let sqrt_delta = delta.sqrt();
                result.push((-c - sqrt_delta) / (2.0.into() * b));
                result.push((-c + sqrt_delta) / (2.0.into() * b));
            } else if delta.abs() < epsilon {
                result.push(-c / (2.0.into() * b));
            }
            return result;
        }

        // is cubic equation -> use cardano's algorithm
        let frac_1_3 = 1.0.into() / 3.0.into();

        let bn = b / a;
        let cn = c / a;
        let dn = d / a;
    
        let delta0 = (3.0.into() * cn - bn * bn) / 9.0.into();
        let delta1 = (9.0.into() * bn * cn - 27.0.into() * dn - 2.0.into() * bn * bn * bn) / 54.0.into();
        let delta_01 = delta0 * delta0 * delta0 + delta1 * delta1;
    
        if delta_01 >= 0.0.into() {
            let delta_p_sqrt = delta1 + delta_01.sqrt();
            let delta_m_sqrt = delta1 - delta_01.sqrt();
    
            let s = delta_p_sqrt.signum() * delta_p_sqrt.abs().powf(frac_1_3);
            let t = delta_m_sqrt.signum() * delta_m_sqrt.abs().powf(frac_1_3);
    
            result.push(-bn * frac_1_3 + (s + t));
    
            // Don't add the repeated root when s + t == 0.
            if (s - t).abs() < epsilon && (s + t).abs() >= epsilon {
                result.push(-bn * frac_1_3 - (s + t) / 2.0.into());
            }
        } else {
            let theta = (delta1 / (-delta0 * delta0 * delta0).sqrt()).acos();
            let two_sqrt_delta0 = 2.0.into() * (-delta0).sqrt();
            result.push(two_sqrt_delta0 * Float::cos(theta * frac_1_3) - bn * frac_1_3);
            result.push(
                two_sqrt_delta0 * Float::cos((theta + 2.0.into() * pi) * frac_1_3) - bn * frac_1_3,
            );
            result.push(
                two_sqrt_delta0 * Float::cos((theta + 4.0.into() * pi) * frac_1_3) - bn * frac_1_3,
            );
        }
    
        result
    }

    /// Return the parameter values corresponding to a given x coordinate.
    /// See also solve_t_for_x for monotonic curves.
    pub fn solve_t_for_x<F>(&self, x: F) -> ArrayVec<[F; 3]> 
    where
    F: Float,
    P:  Sub<P, Output = P>
        + Add<P, Output = P>
        + Mul<F, Output = P>,
    NativeFloat: Sub<F, Output = F> 
        + Add<F, Output = F>
        + Mul<F, Output = F>
        + Float
        + Into<F>
    {
        if self.is_a_point(0.0.into())
            || (self.non_point_is_linear(0.0.into()) && self.start.x() == self.end.x())
        {
            return ArrayVec::new();
        }

        self.solve_t_for_xy(x, self.start.x().into(), self.ctrl1.x().into(), self.ctrl2.x().into(), self.end.x().into())
    }

    /// Return the parameter values corresponding to a given y coordinate.
    /// See also solve_t_for_y for monotonic curves.
    pub fn solve_t_for_y<F>(&self, y: F) -> ArrayVec<[F; 3]> 
    where
    F: Float,
    P:  Sub<P, Output = P>
        + Add<P, Output = P>
        + Mul<F, Output = P>,
    NativeFloat: Sub<F, Output = F> 
        + Add<F, Output = F>
        + Mul<F, Output = F>
        + Float
        + Into<F> 
    {
        if self.is_a_point(0.0.into())
            || (self.non_point_is_linear(0.0.into()) && self.start.y() == self.end.y())
        {
            return ArrayVec::new();
        }

        self.solve_t_for_xy(y, self.start.y().into(), self.ctrl1.y().into(), self.ctrl2.y().into(), self.end.y().into())
    }

    /// Solves the cubic bezier function given the control points' x OR y values
    /// by solving the roots for x or y axis functions
    /// Returns those roots of the function that are in the interval [0.0, 1.0].
    /// This function is not exposed, it has wrappers solve_t_for_x() and solve_t_for_y()
    fn solve_t_for_xy<F>(
        &self,
        value: F,
        from: F,
        ctrl1: F,
        ctrl2: F,
        to: F,
    ) -> ArrayVec<[F; 3]> 
    where
    F: Float,
    P:  Sub<P, Output = P>
        + Add<P, Output = P>
        + Mul<F, Output = P>,
    NativeFloat: Sub<F, Output = F> 
        + Add<F, Output = F>
        + Mul<F, Output = F>
        + Float
        + Into<F>
    {
        let mut result = ArrayVec::new();

        let a = -from + 3.0.into() * ctrl1 - 3.0.into() * ctrl2 + to;
        let b = 3.0.into() * from - 6.0.into() * ctrl1 + 3.0.into() * ctrl2;
        let c = -3.0.into() * from + 3.0.into() * ctrl1;
        let d = from - value;

        let roots = self.real_roots(a, b, c, d);
        for root in roots {
            if root > 0.0.into() && root < 1.0.into() {
                result.push(root);
            }
        }

        result
    }

    /// Return the bounding box of the curve as two tuples ( (xmin, ymin), (xmax, ymax) )
    pub fn bounding_box<F>(&self) -> ((F,F), (F,F)) 
    where
    F: Float,
    P:  Sub<P, Output = P>
        + Add<P, Output = P>
        + Mul<F, Output = P>,
    NativeFloat: Sub<F, Output = F> 
        + Add<F, Output = F>
        + Mul<F, Output = F>
        + Float
        + Into<F>
    {
        // calculate coefficients for the derivative: at^2 + bt + c
        // from the expansion of the cubic bezier curve: sum_i=0_to_3( binomial(3, i) * t^i * (1-t)^(n-i) )
        // yields coeffcients
        // po: [1, -2,  1]
        // p1: [0,  2, -2]
        // p2: [0,  0,  1]
        //      c   b   a

        let derivative = self.derivative();
        // calculate coefficients for derivative
        let a = derivative.start + derivative.ctrl * -2.0.into() + self.end;
        let b = derivative.start * -2.0.into() + derivative.ctrl * 2.0.into();
        let c = self.start;

        // calculate roots for t over x axis and plug them into the bezier function
        //  to get x,y values (make vec 2 bigger for t=0,t=1 values)
        let mut xtremities: ArrayVec<[F; 4]> = ArrayVec::new();
        xtremities.extend(derivative.real_roots(
                                            a.x().into(), 
                                            b.x().into(), 
                                            c.x().into(), 
                                            ).into_iter()
                        );
        // only retain roots for which t is in [0..1] 
        xtremities.retain(|root| -> bool {root > &mut 0.0.into() && root < &mut 1.0.into()});
        // evaluates roots in original function
        for t in xtremities.iter_mut() {
            *t = self.eval_casteljau(*t).x().into();
        }
        // add y-values for start and end point as candidates
        xtremities.push(self.start.x().into()); 
        xtremities.push(self.end.x().into());
        // sort to get min and max values for bounding box
        xtremities.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap());

        // same for y...
        let mut ytremities: ArrayVec<[F; 4]> = ArrayVec::new();
        ytremities.extend(derivative.real_roots(
                                            a.y().into(), 
                                            b.y().into(), 
                                            c.y().into(), 
                                            ).into_iter()
                        );
        // only retain roots for which t is in [0..1] 
        ytremities.retain(|root| -> bool {root > &mut 0.0.into() && root < &mut 1.0.into()});
        // evaluates roots in original function
        for t in ytremities.iter_mut() {
            *t = self.eval_casteljau(*t).y().into();
        }
        // add y-values for start and end point as candidates
        ytremities.push(self.start.y().into()); 
        ytremities.push(self.end.y().into());
        // sort to get min and max values for bounding box
        ytremities.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap());
        // determine xmin, xmax, ymin, ymax, from the set {B(xroots), B(yroots), B(0), B(1)} 
        // (Intermediate control points can't form a boundary)
        let min = (xtremities[0], ytremities[0]);
        let max = (*xtremities.last().unwrap(), *ytremities.last().unwrap()); // can never be empty as we at least have the endpoints
        return (min, max)
    }

}


#[cfg(test)]
mod tests 
{
    use super::*;
    use crate::num_traits::{Pow};
    #[test]
    fn circle_approximation_error() 
    {
        // define closure for unit circle 
        let circle = |p: Point2<f64>| -> f64 { ( p.x.pow(2) as f64 
                                                + p.y.pow(2) as f64)
                                                .sqrt() - 1f64};

        // define control points for 4 bezier segments 
        // control points are chosen for minimum radial distance error
        // according to: http://spencermortensen.com/articles/bezier-circle/ 
        // TODO don't hardcode values
        let c               = 0.551915024494;
        let max_drift_perc  = 0.019608; // radial drift percent
        let max_error       = max_drift_perc * 0.01; // absolute max radial error

        let bezier_quadrant_1= CubicBezier{ start:  Point2{x:0f64,  y:1f64},
                                                ctrl1: Point2{x:c,     y:1f64},
                                                ctrl2: Point2{x:1f64,  y:c},
                                                end:   Point2{x:1f64,  y:0f64}};
        let bezier_quadrant_2 = CubicBezier{ start:  Point2{x:1f64,  y:0f64},
                                                ctrl1: Point2{x:1f64,     y:-c},
                                                ctrl2: Point2{x:c,  y:-1f64},
                                                end:   Point2{x:0f64,  y:-1f64}};
        let bezier_quadrant_3 = CubicBezier{ start:  Point2{x:0f64,  y:-1f64},
                                                ctrl1: Point2{x:-c,     y:-1f64},
                                                ctrl2: Point2{x:-1f64,  y:-c},
                                                end:   Point2{x:-1f64,  y:0f64}};
        let bezier_quadrant_4 = CubicBezier{ start:  Point2{x:-1f64,    y:0f64},
                                                ctrl1: Point2{x:-1f64,  y:c},
                                                ctrl2: Point2{x:-c,     y:1f64},
                                                end:   Point2{x:0f64,   y:1f64}};
        let nsteps =  1000;                                      
        for t in 0..nsteps {
            let t = t as f64 * 1f64/(nsteps as f64);

            let point = bezier_quadrant_1.eval(t);
            let contour = circle(point);
            assert!( contour.abs() <= max_error );

            let point = bezier_quadrant_2.eval(t);
            let contour = circle(point);
            assert!( contour.abs() <= max_error );

            let point = bezier_quadrant_3.eval(t);
            let contour = circle(point);
            assert!( contour.abs() <= max_error );

            let point = bezier_quadrant_4.eval(t);
            let contour = circle(point);
            assert!( contour.abs() <= max_error );
        }
    }


    #[test]
    fn circle_circumference_approximation() 
    {
        // define control points for 4 cubic bezier segments to best approximate a unit circle
        // control points are chosen for minimum radial distance error, see circle_approximation_error() in this file
        // given this, the circumference will also be close to 2*pi 
        // (remember arclen also works by linear approximation, not the true integral, so we have to accept error)!
        // This approximation is unfeasable if desired accuracy is greater than 2 decimal places (at 1000 steps)
        // TODO don't hardcode values, solve for them
        let c         = 0.551915024494;
        let max_error = 1e-2;
        let nsteps  = 1e3 as usize;
        let pi        = 3.14159265359;
        let tau       = 2.*pi;

        let bezier_quadrant_1= CubicBezier{ start:  Point2{x:0f64,  y:1f64},
                                                ctrl1: Point2{x:c,     y:1f64},
                                                ctrl2: Point2{x:1f64,  y:c},
                                                end:   Point2{x:1f64,  y:0f64}};
        let bezier_quadrant_2 = CubicBezier{ start:  Point2{x:1f64,  y:0f64},
                                                ctrl1: Point2{x:1f64,     y:-c},
                                                ctrl2: Point2{x:c,  y:-1f64},
                                                end:   Point2{x:0f64,  y:-1f64}};
        let bezier_quadrant_3 = CubicBezier{ start:  Point2{x:0f64,  y:-1f64},
                                                ctrl1: Point2{x:-c,     y:-1f64},
                                                ctrl2: Point2{x:-1f64,  y:-c},
                                                end:   Point2{x:-1f64,  y:0f64}};
        let bezier_quadrant_4 = CubicBezier{ start:  Point2{x:-1f64,    y:0f64},
                                                ctrl1: Point2{x:-1f64,  y:c},
                                                ctrl2: Point2{x:-c,     y:1f64},
                                                end:   Point2{x:0f64,   y:1f64}};
        let circumference = bezier_quadrant_1.arclen::<NativeFloat>(nsteps) +
                                bezier_quadrant_2.arclen::<NativeFloat>(nsteps) +
                                bezier_quadrant_3.arclen::<NativeFloat>(nsteps) +
                                bezier_quadrant_4.arclen::<NativeFloat>(nsteps);
        //dbg!(circumference);
        //dbg!(tau);
        assert!( ((tau + max_error) > circumference) && ((tau - max_error) < circumference) );
    }

    #[test]
    fn eval_equivalence() {
        // all eval methods should be approximately equivalent for well defined test cases
        // and not equivalent where numerical stability becomes an issue for normal eval
        let bezier = CubicBezier::new( 
            Point2::new(0f64,  1.77f64),
            Point2::new(1.1f64, -1f64),
            Point2::new(4.3f64,3f64),
            Point2::new(3.2f64, -4f64),
        );

        let max_err = 1e-14;
        let nsteps: usize =  1000;                                      
        for t in 0..nsteps {
            let t = t as f64 * 1f64/(nsteps as f64);
            let p1 = bezier.eval(t);
            let p2 = bezier.eval_casteljau(t);
            let err = p2-p1;
            //dbg!(p1);
            //dbg!(p2);
            assert!( (err.x.abs() < max_err) && (err.y.abs() < max_err) );
        }
    }

    #[test]
    fn split_equivalence() {
        // chose some arbitrary control points and construct a cubic bezier
        let bezier = CubicBezier{ 
                                start:  Point2{x:0f64,  y:1.77f64},
                                ctrl1: Point2{x:2.9f64, y:0f64},
                                ctrl2: Point2{x:4.3f64, y:3f64},
                                end:   Point2{x:3.2f64, y:-4f64}
                            };
        // split it at an arbitrary point
        let at = 0.5;
        let (left, right) = bezier.split(at);
        // compare left and right subcurves with parent curve
        // this is tricky as we have to map t->t/2 (for left) which will 
        // inevitably contain rounding errors from floating point ops.
        // instead, take the difference of the two points which must not exceed the absolute error
        // TODO update test to use norm() instead, once implemented for Point (maybe as trait?)
        let max_err = 1e-14;
        let nsteps: usize =  1000;                                      
        for t in 0..nsteps {
            let t = t as f64 * 1f64/(nsteps as f64);
            //dbg!(bezier.eval(t/2.0));
            //dbg!(left.eval(t));
            // left
            let mut err = bezier.eval(t/2.0) - left.eval(t);
            //dbg!(err);
            assert!( (err.x.abs() < max_err) && (err.y.abs() < max_err) );
            // right
            err = bezier.eval((t*0.5)+0.5) - right.eval(t);
            //dbg!(err);
            assert!( (err.x.abs() < max_err) && (err.y.abs() < max_err) );  
        }
    }


    #[test]
    fn bounding_box_contains() {
        // check if bounding box for a curve contains all points (with some approximation error)
        let bezier = CubicBezier{ start:  Point2{x:0f64,  y:1.77f64},
                                  ctrl1: Point2{x:2.9f64, y:0f64},
                                  ctrl2: Point2{x:4.3f64, y:-3f64},
                                  end:   Point2{x:3.2f64,  y:4f64}};

        let ((xmin, ymin), (xmax, ymax)) = bezier.bounding_box::<f64>();

        let max_err = 1e-2;

        let nsteps: usize =  100;                                      
        for t in 0..nsteps {
            let t = t as f64 * 1f64/(nsteps as f64);
            let p = bezier.eval_casteljau(t);
            dbg!(t);
            dbg!(p);
            dbg!(xmin-max_err, ymin-max_err, xmax+max_err, ymax+max_err);

            assert!( (p.x() >= (xmin-max_err) ) && (p.y() >= (ymin-max_err)) );
            assert!( (p.x() <= (xmax+max_err) ) && (p.y() <= (ymax+max_err)) );

        }
    }
}