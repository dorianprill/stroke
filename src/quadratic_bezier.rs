use core::default::Default;

use super::*;
use super::point::Point;
use super::line::LineSegment; 
//use super::cubic_bezier::CubicBezier;

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
    + Point<Scalar = NativeFloat>,
{

    /// Creates a new Instance of QuadraticBezier from the given control points
    pub fn new(start: P, ctrl: P, end: P) -> Self {
        QuadraticBezier { 
            start, 
            ctrl, 
            end 
        }
    }

    /// Evaluates the quadratic bezier curve at t using direct evaluation, which may not be numerically stable
    pub fn eval<F>(&self, t: F) -> P 
    where 
    F: Float,
    P: Add<P, Output = P>
        + Sub<P, Output = P>
        + Mul<F, Output = P>
        + Point,
    NativeFloat: Sub<F, Output = F> 
        + Mul<F, Output = F> 
        + Into<F>
    {
        let t2          = t * t;
        let one_t  = 1.0 - t;
        let one_t2      = one_t * one_t;

        self.start * one_t2
            + self.ctrl * 2.0 as NativeFloat * one_t * t
            + self.end * t2
    }

    /// Evaluates the cubic bezier curve at t using the numerically stable De Casteljau algorithm
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
        let ctrl_1ab = self.start + (self.ctrl - self.start) * t;
        let ctrl_1bc   = self.ctrl + (self.end - self.ctrl) * t;
        // second iteration, final point on the curve
        let ctrl_2ab  = ctrl_1ab + (ctrl_1bc - ctrl_1ab) * t;
    
        return ctrl_2ab
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
        let ctrl_1ab = self.start + (self.ctrl - self.start) * t;
        let ctrl_1bc   = self.ctrl + (self.end - self.ctrl) * t;
        // second iteration
        let ctrl_2ab = ctrl_1ab + (ctrl_1bc - ctrl_1ab) * t;

        return (
            QuadraticBezier {
                start: self.start,
                ctrl: ctrl_1ab,
                end: ctrl_2ab,
            },
            QuadraticBezier {
                start: ctrl_2ab,
                ctrl: ctrl_1bc,
                end: self.end,
            },
        );
    }

    /// Sample the x coordinate of the curve at t (expecting t between 0 and 1).
    pub fn axis<F>(&self, t: F, axis: usize) -> F 
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

        return self.start.axis(axis) * one_t2 
            + self.ctrl.axis(axis) * 2.0.into() * one_t * t 
            + self.end.axis(axis) * t2
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
        + Float
        + Into<F>
    {
        return LineSegment{
            start: (self.ctrl - self.start) * 2.0.into(),
            end:   (self.end - self.ctrl)   * 2.0.into()
        }
    }


    /// Sample the x coordinate of the curve's derivative at t (expecting t between 0 and 1).
    /// Convenience function for .derivative().eval(t).axis(dim)
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
        return self.start.axis(0).into() * c0 + self.ctrl.axis(0).into() * c1 + self.end.axis(0).into() * c2;
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
        return self.start.axis(1).into() * c0 + self.ctrl.axis(1).into() * c1 + self.end.axis(1).into() * c2;
    }

    // /// Calculates the curvature of the curve at point t
    // /// The curvature is the inverse of the radius of the tangential circle at t: k=1/r
    // pub fn curvature<F>(&self, t: F) -> F
    // where
    // F: Float,
    // P:  Sub<P, Output = P>
    //     + Add<P, Output = P>
    //     + Mul<F, Output = P>,
    // NativeFloat: Sub<F, Output = F> 
    //     + Add<F, Output = F>
    //     + Mul<F, Output = F>
    //     + Into<F>
    // {
    //     let d = self.derivative();
    //     let dd = d.derivative();
    //     let dx = d.x(t);
    //     let dy = d.y(t);
    //     let (ddx, ddy) = dd;
    //     let numerator = dx * ddy.into() - ddx * dy;
    //     let denominator = (dx*dx + dy*dy).powf(1.5.into());
    //     return numerator / denominator
    // }

    // /// Calculates the radius of the tangential circle at t
    // /// It is the inverse of the curvature at t: r=1/k
    // pub fn radius<F>(&self, t: F) -> F
    // where
    // F: Float,
    // P:  Sub<P, Output = P>
    //     + Add<P, Output = P>
    //     + Mul<F, Output = P>,
    // NativeFloat: Sub<F, Output = F> 
    //     + Add<F, Output = F>
    //     + Mul<F, Output = F>
    //     + Into<F>
    // {
    //     return 1.0.into() / self.curvature(t)
    // }


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
            let t = t as NativeFloat * 1.0.into()/(nsteps as NativeFloat).into();
            let p1 = self.eval_casteljau(t);
            let p2 = self.eval_casteljau(t+stepsize.into());

            arclen = arclen + (p1-p2).squared_length().sqrt();
        
        }
        return arclen.into()
    }


    /// Solve for the roots of the polynomial at^2 + bt + c
    /// Returns an ArrayVec of roots in the order
    /// needs to be called for x and y components separately
    pub(crate) fn real_roots<F>(&self, a: F, b: F, c: F) -> ArrayVec<[F; 2]>
    where
    F:  Float
        + Default,
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

        // check if can be handled below quadratic order
        if a.abs() < EPSILON.into() {
            if b.abs() < EPSILON.into() {
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
        } else if delta.abs() < EPSILON.into() {
            result.push(-b / (2.0.into() * a));
        }
        return result;
    }

    pub fn baseline(&self) -> LineSegment<P> {
        LineSegment {
            start: self.start,
            end: self.end,
        }
    }



    pub fn is_linear<F>(&self, tolerance: F) -> bool 
    where
    F: Float
        + Default
        + Into<NativeFloat>,
    P: Mul<F, Output = P>,
    NativeFloat: Sub<F, Output = F> 
        + Add<F, Output = F>
        + Mul<F, Output = F>
        + Float
        + Into<F> 
    {
        // if start and end are (nearly) the same
        // TODO using squred length vs machine epsilon OK?
        if (self.start-self.end).squared_length() < EPSILON {
            return false;
        } 
        // else check if ctrl points lie on baseline i.e. all points are colinear
        self.are_points_colinear(tolerance)
    }



    fn are_points_colinear<F>(&self, tolerance: F) -> bool 
    where
    F: Float
        + Default
        + Into<NativeFloat>,
    P: Mul<F, Output = P>,
    NativeFloat: Sub<F, Output = F> 
        + Add<F, Output = F>
        + Mul<F, Output = F>
        + Float
        + Into<F>
    {
        let lineeq = self.baseline();//.to_line().equation();
        lineeq.distance_to_point(self.ctrl) <= tolerance
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
        (self.start-self.end).squared_length().into() <= tolerance_squared
            && (self.start-self.ctrl).squared_length().into() <= tolerance_squared
    }



    /// Solves the cubic bezier function given the control points' x OR y values
    /// by solving the roots for the function over any of the axis
    /// Returns those roots of the function that are in the interval [0.0, 1.0].
    fn solve_t_for_axis<F>(&self, value: F, axis: usize) -> ArrayVec<[F; 3]> 
    where
    F:  Float
        + Default
        + Into<NativeFloat>,
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
        if self.is_a_point(EPSILON.into())
            || (self.are_points_colinear(0.0.into()))// && self.start.y() == self.end.y())
        {
            return result
        }
        // these are just the x or y components of the points
        let a = self.start.axis(axis) 
                    + self.ctrl.axis(axis) * -2.0.into() 
                    + self.end.axis(axis).into();
        let b = self.start.axis(axis) * -2.0.into() 
                    + self.ctrl.axis(axis) * 2.0.into();
        let c = self.start.axis(axis) - value;

        let roots = self.real_roots(a, b, c);
        for root in roots {
            if root > 0.0.into() && root < 1.0.into() {
                result.push(root);
            }
        }

        result
    }

    /// Return the bounding box of the curve as an array of (min, max) tuples for each dimension (its index)
    pub fn bounding_box<F>(&self) -> [(F,F); P::DIM]
    where
    F:  Float
        + Default,
    P:  Sub<P, Output = P>
        + Add<P, Output = P>
        + Mul<F, Output = P>,
    NativeFloat: Sub<F, Output = F> 
        + Add<F, Output = F>
        + Mul<F, Output = F>
        + Float
        + Into<F>
    {
        let mut bounds = [(0.0.into(), 0.0.into()); P::DIM];
        let derivative = self.derivative();
        // calculate coefficients for the derivative as a function of t: at + b
        // po: [1, -1]
        // p1: [0,  1]
        //      b   a
        let a = derivative.start * -1.0.into() + derivative.end;
        let b = derivative.start;

        for (dim, _) in a.into_iter().enumerate() {
            // calculate roots for t over x axis and plug them into the bezier function
            //  to get x,y values (make vec 2 bigger for t=0,t=1 values)
            let mut extrema: ArrayVec<[F; 3]> = ArrayVec::new();
            extrema.extend(derivative.root(a.axis(dim).into(),
                                                    b.axis(dim).into())
                                                        .into_iter());
            // only retain roots for which t is in [0..1] 
            extrema.retain(|root| -> bool {root > &mut 0.0.into() && root < &mut 1.0.into()});
            // evaluates roots in original function
            for t in extrema.iter_mut() {
                *t = self.eval_casteljau(*t).axis(dim).into();
            }
            // add y-values for start and end point as candidates
            extrema.push(self.start.axis(dim).into()); 
            extrema.push(self.end.axis(dim).into());
            // sort to get min and max values for bounding box
            extrema.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap());

            // determine xmin, xmax, ymin, ymax, from the set {B(xroots), B(yroots), B(0), B(1)} 
            // (Intermediate control points can't form a boundary)
            // unwrap() is ok as it always at least contains the endpoints
            bounds[dim] = (extrema[0], *extrema.last().unwrap()); 
        }
        return bounds
    }

}

#[cfg(test)]
mod tests 
{
    use super::*;
    //use crate::num_traits::{Pow};
    use super::point_generic::PointN;
    //TODO test needs to be adapted for 8 segments of quadratic order
    // #[test]
    // fn circle_approximation_error() 
    // {
    //     // define closure for unit circle 
    //     let circle = |p: Point2<f64>| -> f64 { ( p.x.pow(2) as f64 
    //                                             + p.y.pow(2) as f64)
    //                                             .sqrt() - 1f64};

    //     // define control points for 4 bezier segments 
    //     // control points are chosen for minimum radial distance error
    //     // according to: http://spencermortensen.com/articles/bezier-circle/ 
    //     // TODO don't hardcode values
    //     let c               = 0.551915024494;
    //     let max_drift_perc  = 0.019608; // radial drift percent
    //     let max_error       = max_drift_perc * 0.01; // absolute max radial error

    //     let bezier_quadrant_1= QuadraticBezier{ start:  Point2{x:0f64,  y:1f64},
    //                                             ctrl: Point2{x:1f64,  y:c},
    //                                             end:   Point2{x:1f64,  y:0f64}};
    //     let bezier_quadrant_2 = QuadraticBezier{ start:  Point2{x:1f64,  y:0f64},
    //                                             ctrl: Point2{x:c,  y:-1f64},
    //                                             end:   Point2{x:0f64,  y:-1f64}};
    //     let bezier_quadrant_3 = QuadraticBezier{ start:  Point2{x:0f64,  y:-1f64},
    //                                             ctrl: Point2{x:-1f64,  y:-c},
    //                                             end:   Point2{x:-1f64,  y:0f64}};
    //     let bezier_quadrant_4 = QuadraticBezier{ start:  Point2{x:-1f64,    y:0f64},
    //                                             ctrl: Point2{x:-c,     y:1f64},
    //                                             end:   Point2{x:0f64,   y:1f64}};
    //     let nsteps =  1000;                                      
    //     for t in 0..nsteps {
    //         let t = t as f64 * 1f64/(nsteps as f64);

    //         let point = bezier_quadrant_1.eval(t);
    //         let contour = circle(point);
    //         assert!( contour.abs() <= max_error );

    //         let point = bezier_quadrant_2.eval(t);
    //         let contour = circle(point);
    //         assert!( contour.abs() <= max_error );

    //         let point = bezier_quadrant_3.eval(t);
    //         let contour = circle(point);
    //         assert!( contour.abs() <= max_error );

    //         let point = bezier_quadrant_4.eval(t);
    //         let contour = circle(point);
    //         assert!( contour.abs() <= max_error );
    //     }
    // }


    
    //TODO test needs to be adapted for 8 segments of quadratic order
    // #[test]
    // fn circle_circumference_approximation() 
    // {
    //     // define control points for 8 quadratic bezier segments to best approximate a unit circle
    //     // control points are chosen for minimum radial distance error, see circle_approximation_error() in this file
    //     // given this, the circumference will also be close to 2*pi 
    //     // (remember arclen also works by linear approximation, not the true integral, so we have to accept error)!
    //     // This approximation is unfeasable if desired accuracy is greater than 2 decimal places (at 1000 steps)
    //     // TODO don't hardcode values, solve for them
    //     let c         = 0.551915024494;
    //     let max_error = 1e-2;
    //     let nsteps  = 1e3 as usize;
    //     let pi        = 3.14159265359;
    //     let tau       = 2.*pi;

    //     let bezier_quadrant_1= QuadraticBezier{ start:  Point2{x:0f64,  y:1f64},
    //                                             ctrl: Point2{x:1f64,  y:c},
    //                                             end:   Point2{x:1f64,  y:0f64}};
    //     let bezier_quadrant_2 = QuadraticBezier{ start:  Point2{x:1f64,  y:0f64},
    //                                             ctrl: Point2{x:c,  y:-1f64},
    //                                             end:   Point2{x:0f64,  y:-1f64}};
    //     let bezier_quadrant_3 = QuadraticBezier{ start:  Point2{x:0f64,  y:-1f64},
    //                                             ctrl: Point2{x:-1f64,  y:-c},
    //                                             end:   Point2{x:-1f64,  y:0f64}};
    //     let bezier_quadrant_4 = QuadraticBezier{ start:  Point2{x:-1f64,    y:0f64},
    //                                             ctrl: Point2{x:-c,     y:1f64},
    //                                             end:   Point2{x:0f64,   y:1f64}};
    //     let circumference = bezier_quadrant_1.arclen::<NativeFloat>(nsteps) +
    //                             bezier_quadrant_2.arclen::<NativeFloat>(nsteps) +
    //                             bezier_quadrant_3.arclen::<NativeFloat>(nsteps) +
    //                             bezier_quadrant_4.arclen::<NativeFloat>(nsteps);
    //     //dbg!(circumference);
    //     //dbg!(tau);
    //     assert!( ((tau + max_error) > circumference) && ((tau - max_error) < circumference) );
    // }

    #[test]
    fn eval_equivalence() {
        // all eval methods should be approximately equivalent for well defined test cases
        // and not equivalent where numerical stability becomes an issue for normal eval
        let bezier = QuadraticBezier::new( 
            PointN::new([0f64,  1.77f64]),
            PointN::new([4.3f64, 3f64]),
            PointN::new([3.2f64, -4f64]),
        );

        let max_err = 1e-14;
        let nsteps: usize =  1000;                                      
        for t in 0..=nsteps {
            let t = t as f64 * 1f64/(nsteps as f64);
            let p1 = bezier.eval(t);
            let p2 = bezier.eval_casteljau(t);
            let err = p2-p1;
            //dbg!(p1);
            //dbg!(p2);
            for axis in err {
                assert!(axis.abs() < max_err);
            }
        }
    }

    #[test]
    fn split_equivalence() {
        // chose some arbitrary control points and construct a cubic bezier
        let bezier = QuadraticBezier{ 
                        start:  PointN::new([0f64, 1.77f64]),
                        ctrl:   PointN::new([4.3f64, 3f64]),
                        end:    PointN::new([3.2f64, -4f64])
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
        let nsteps: usize = 1000;                                      
        for t in 0..=nsteps {
            let t = t as f64 * 1f64/(nsteps as f64);
            //dbg!(bezier.eval(t/2.0));
            //dbg!(left.eval(t));
            // left
            let mut err = bezier.eval(t/2.0) - left.eval(t);
            for axis in err {
                assert!(axis.abs() < max_err);
            }
            // right
            err = bezier.eval((t*0.5)+0.5) - right.eval(t);
            //dbg!(err);
            for axis in err {
                assert!(axis.abs() < max_err);
            }
        }
    }


    #[test]
    fn bounding_box_contains() {
        // check if bounding box for a curve contains all points (with some approximation error)
        let bezier = QuadraticBezier{ 
                    start:  PointN::new([0f64, 1.77f64]),
                    ctrl:   PointN::new([4.3f64, -3f64]),
                    end:    PointN::new([3.2f64, 4f64])
        };

        let bounds = bezier.bounding_box::<f64>();

        let max_err = 1e-2;

        let nsteps: usize =  100;                                      
        for t in 0..=nsteps {
            let t = t as f64 * 1f64/(nsteps as f64);
            let p = bezier.eval_casteljau(t);
            //dbg!(t);
            //dbg!(p);
            //dbg!(xmin-max_err, ymin-max_err, xmax+max_err, ymax+max_err);

            for (idx, axis) in p.into_iter().enumerate() {
                assert!( (axis >= (bounds[idx].0 - max_err)) && (axis <= (bounds[idx].1 + max_err)) )
            }

        }
    }
}
