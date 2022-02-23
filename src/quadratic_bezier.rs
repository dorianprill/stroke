use super::LineSegment;
use super::Point;
use super::*;

#[derive(Copy, Clone, Debug, PartialEq)]
pub struct QuadraticBezier<P: Point> {
    pub(crate) start: P,
    pub(crate) ctrl: P,
    pub(crate) end: P,
}

impl<P: Point> QuadraticBezier<P>
where
    P: Point,
{
    /// Creates a new instance of QuadraticBezier from the given control points
    pub fn new(start: P, ctrl: P, end: P) -> Self {
        QuadraticBezier { start, ctrl, end }
    }

    /// Evaluates the quadratic bezier curve at 't' using direct evaluation, which may not be numerically stable
    pub fn eval(&self, t: P::Scalar) -> P {
        let t2 = t * t;
        let one_t = -t + 1.0;
        let one_t2 = one_t * one_t;

        self.start * one_t2 + self.ctrl * 2.0 * one_t * t + self.end * t2
    }

    /// Evaluates the cubic bezier curve at t using the numerically stable De Casteljau algorithm
    pub fn eval_casteljau(&self, t: P::Scalar) -> P {
        // unrolled de casteljau algorithm
        // _1ab is the first iteration from first (a) to second (b) control point and so on
        let ctrl_1ab = self.start + (self.ctrl - self.start) * t;
        let ctrl_1bc = self.ctrl + (self.end - self.ctrl) * t;
        // second iteration, return final point on the curve ctrl_2ab
        ctrl_1ab + (ctrl_1bc - ctrl_1ab) * t
    }

    pub fn control_points(&self) -> [P; 3] {
        [self.start, self.ctrl, self.end]
    }

    pub fn split(&self, t: P::Scalar) -> (Self, Self) {
        // unrolled de casteljau algorithm
        // _1ab is the first iteration from first (a) to second (b) control point and so on
        let ctrl_1ab = self.start + (self.ctrl - self.start) * t;
        let ctrl_1bc = self.ctrl + (self.end - self.ctrl) * t;
        // second iteration
        let ctrl_2ab = ctrl_1ab + (ctrl_1bc - ctrl_1ab) * t;

        (
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
        )
    }

    /// Sample the a particular coordinate axis of the curve at t (expecting t between 0 and 1).
    /// Shortcut for curve.eval(t).axis(k)
    /// This function can panic! 
    /// TODO may add something like const_assert for Point's const DIM
    pub fn axis(&self, t: P::Scalar, axis: usize) -> P::Scalar {
        let t2 = t * t;
        let one_t = -t + 1.0;
        let one_t2 = one_t * one_t;

        self.start.axis(axis) * one_t2
            + self.ctrl.axis(axis) * 2.0 * one_t * t
            + self.end.axis(axis) * t2
    }

    /// Return the derivative curve.
    /// The derivative is also a bezier curve but of degree n-1.
    /// In the case of a quadratic derivative it is just a line segment
    /// which also implementes eval(), as it is just a linear bezier curve.
    pub fn derivative(&self) -> LineSegment<P> {
        LineSegment {
            start: (self.ctrl - self.start) * 2.0,
            end: (self.end - self.ctrl) * 2.0,
        }
    }

    /// Direct Derivative - Sample the axis coordinate at 'axis' of the curve's derivative at t
    /// without creating a new curve. This is a convenience function for .derivative().eval(t).axis(n)  
    /// Parameters:
    ///   t: the sampling parameter on the curve interval [0..1]
    ///   axis: the index of the coordinate axis [0..N]
    /// Returns:
    ///   Scalar value of the points own type type F  
    /// May be deprecated in the future.  
    /// This function can cause out of bounds panic when axis is larger than dimension of P
    pub fn dd(&self, t: P::Scalar, axis: usize) -> P::Scalar {
        let t = t.into();
        let c0 = t * 2.0 - 2.0;
        let c1 = 2.0 - 4.0 * t;
        let c2 = 2.0 * t;

        self.start.axis(axis) * c0 + self.ctrl.axis(axis) * c1 + self.end.axis(axis) * c2
    }

    // /// Calculates the curvature of the curve at point t
    // /// The curvature is the inverse of the radius of the tangential circle at t: k=1/r
    // pub fn curvature(&self, t: P::Scalar) -> F
    // where
    // F: P::Scalarloat,
    // P::Scalar: Sub<F, Output = F>
    //     + Add<F, Output = F>
    //     + Mul<F, Output = F>
    //     + Into
    //     + From
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
    // pub fn radius(&self, t: P::Scalar) -> F
    // where
    // F: P::Scalarloat,
    // P::Scalar: Sub<F, Output = F>
    //     + Add<F, Output = F>
    //     + Mul<F, Output = F>
    //     + Into
    //     + From
    // {
    //     return 1.0.into() / self.curvature(t)
    // }

    /// Approximates the arc length of the curve by flattening it with straight line segments.
    /// This works quite well, at ~32 segments it should already provide an error < 0.5
    /// Remember arclen also works by linear approximation, not the integral, so we have to accept error!
    /// This approximation is unfeasable if desired accuracy is greater than 2 decimal places
    pub fn arclen(&self, nsteps: usize) -> P::Scalar {
        let stepsize = P::Scalar::from(1.0 / (nsteps as NativeFloat));
        let mut arclen: P::Scalar = 0.0.into();
        for t in 1..nsteps {
            let t = P::Scalar::from(t as NativeFloat * 1.0 / (nsteps as NativeFloat));
            let p1 = self.eval_casteljau(t);
            let p2 = self.eval_casteljau(t + stepsize);

            arclen = arclen + (p1 - p2).squared_length().sqrt();
        }
        arclen
    }

    /// Solve for the roots of the polynomial at^2 + bt + c
    /// Returns an ArrayVec of roots in the order
    /// needs to be called for x and y components separately
    pub(crate) fn real_roots(
        &self,
        a: P::Scalar,
        b: P::Scalar,
        c: P::Scalar,
    ) -> ArrayVec<[P::Scalar; 2]> {
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
        let delta = b * b - a * c * 4.0;
        if delta > 0.0.into() {
            let sqrt_delta = delta.sqrt();
            result.push((-b - sqrt_delta) / (a * 2.0));
            result.push((-b + sqrt_delta) / (a * 2.0));
        } else if delta.abs() < EPSILON.into() {
            result.push(-b / (a * 2.0));
        }
        result
    }

    /// Returns the line segment formed by the curve's start and endpoint
    pub fn baseline(&self) -> LineSegment<P> {
        LineSegment {
            start: self.start,
            end: self.end,
        }
    }

    /// Checks if, given some tolerance, the curve can be considered equal to a line segment
    pub fn is_linear(&self, tolerance: P::Scalar) -> bool {
        // if start and end are (nearly) the same
        // TODO using squared length vs machine epsilon OK?
        if (self.start - self.end).squared_length() < EPSILON.into() {
            return false;
        }
        // else check if ctrl points lie on baseline i.e. all points are colinear
        self.are_points_colinear(tolerance)
    }

    /// Determines if, given some tolerance, all of the control points are colinear
    /// This private function is wrapped publically by is_linear()
    fn are_points_colinear(&self, tolerance: P::Scalar) -> bool {
        let line = self.baseline();
        line.distance_to_point(self.ctrl) <= tolerance
    }

    /// Determines if, given some tolerance, the control points of the curve can be considered equal.
    /// If true, the curve is just a singular point
    pub fn is_a_point(&self, tolerance: P::Scalar) -> bool {
        let tolerance_squared = tolerance * tolerance;
        // Use <= so that tolerance can be zero.
        (self.start - self.end).squared_length() <= tolerance_squared
            && (self.start - self.ctrl).squared_length() <= tolerance_squared
    }

    /// Solves the quadratic bezier function given a particular coordinate axis value
    /// by solving the roots for the axis functions
    /// Parameters:
    /// value: the coordinate value on the particular axis
    /// axis: the index of the axis
    /// Returns those roots of the function that are in the interval [0.0, 1.0].
    fn solve_t_for_axis(&self, value: P::Scalar, axis: usize) -> ArrayVec<[P::Scalar; 3]> {
        let mut result = ArrayVec::new();
        if self.is_a_point(EPSILON.into())
            || (self.are_points_colinear(0.0.into())
                && (self.start - self.end).squared_length() < EPSILON.into())
        {
            return result;
        }
        // these are just the x or y components of the points
        let a = self.start.axis(axis) + self.ctrl.axis(axis) * -2.0 + self.end.axis(axis);
        let b = self.start.axis(axis) * -2.0 + self.ctrl.axis(axis) * 2.0;
        let c = self.start.axis(axis) - value.into();

        let roots = self.real_roots(a, b, c);
        for root in roots {
            if root > 0.0.into() && root < 1.0.into() {
                result.push(root);
            }
        }

        result
    }

    /// Return the bounding box of the curve as an array of (min, max) tuples for each dimension (its index)
    pub fn bounding_box(&self) -> [(P::Scalar, P::Scalar); P::DIM] {
        let mut bounds = [(0.0.into(), 0.0.into()); P::DIM];
        let derivative = self.derivative();
        // calculate coefficients for the derivative as a function of t: at + b
        // po: [1, -1]
        // p1: [0,  1]
        //      b   a
        let a = derivative.start * -1.0 + derivative.end;
        let b = derivative.start;

        for (dim, _) in a.into_iter().enumerate() {
            // calculate roots for t over x axis and plug them into the bezier function
            //  to get x,y values (make vec 2 bigger for t=0,t=1 values)
            let mut extrema: ArrayVec<[P::Scalar; 3]> = ArrayVec::new();
            extrema.extend(derivative.root(a.axis(dim), b.axis(dim)).into_iter());
            // only retain roots for which t is in [0..1]
            extrema.retain(|root| -> bool { root > &mut 0.0.into() && root < &mut 1.0.into() });
            // evaluates roots in original function
            for t in extrema.iter_mut() {
                *t = self.eval_casteljau(*t).axis(dim);
            }
            // add y-values for start and end point as candidates
            extrema.push(self.start.axis(dim));
            extrema.push(self.end.axis(dim));
            // sort to get min and max values for bounding box
            extrema.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap());

            // determine xmin, xmax, ymin, ymax, from the set {B(xroots), B(yroots), B(0), B(1)}
            // (Intermediate control points can't form a boundary)
            // unwrap() is ok as it always at least contains the endpoints
            bounds[dim] = (extrema[0], *extrema.last().unwrap());
        }
        bounds
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    //use crate::num_traits::{Pow};
    use super::PointN;
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
    //     let circumference = bezier_quadrant_1.arclen::<P::Scalar>(nsteps) +
    //                             bezier_quadrant_2.arclen::<P::Scalar>(nsteps) +
    //                             bezier_quadrant_3.arclen::<P::Scalar>(nsteps) +
    //                             bezier_quadrant_4.arclen::<P::Scalar>(nsteps);
    //     //dbg!(circumference);
    //     //dbg!(tau);
    //     assert!( ((tau + max_error) > circumference) && ((tau - max_error) < circumference) );
    // }

    #[test]
    fn eval_equivalence() {
        // all eval methods should be approximately equivalent for well defined test cases
        // and not equivalent where numerical stability becomes an issue for normal eval
        let bezier = QuadraticBezier::new(
            PointN::new([0f64, 1.77f64]),
            PointN::new([4.3f64, 3f64]),
            PointN::new([3.2f64, -4f64]),
        );

        let nsteps: usize = 1000;
        for t in 0..=nsteps {
            let t = t as f64 * 1f64 / (nsteps as f64);
            let p1 = bezier.eval(t);
            let p2 = bezier.eval_casteljau(t);
            let err = p2 - p1;
            assert!(err.squared_length() < EPSILON);
        }
    }

    #[test]
    fn split_equivalence() {
        // chose some arbitrary control points and construct a cubic bezier
        let bezier = QuadraticBezier {
            start: PointN::new([0f64, 1.77f64]),
            ctrl: PointN::new([4.3f64, 3f64]),
            end: PointN::new([3.2f64, -4f64]),
        };
        // split it at an arbitrary point
        let at = 0.5;
        let (left, right) = bezier.split(at);
        // compare left and right subcurves with parent curve
        // take the difference of the two points which must not exceed the absolute error
        let nsteps: usize = 1000;
        for t in 0..=nsteps {
            let t = t as f64 * 1f64 / (nsteps as f64);
            // left
            let mut err = bezier.eval(t / 2.0) - left.eval(t);
            assert!(err.squared_length() < EPSILON);
            // right
            err = bezier.eval((t * 0.5) + 0.5) - right.eval(t);
            assert!(err.squared_length() < EPSILON);
        }
    }

    #[test]
    fn bounding_box_contains() {
        // check if bounding box for a curve contains all points (with some approximation error)
        let bezier = QuadraticBezier {
            start: PointN::new([0f64, 1.77f64]),
            ctrl: PointN::new([4.3f64, -3f64]),
            end: PointN::new([3.2f64, 4f64]),
        };

        let bounds = bezier.bounding_box();

        let max_err = 1e-2;

        let nsteps: usize = 100;
        for t in 0..=nsteps {
            let t = t as f64 * 1f64 / (nsteps as f64);
            let p = bezier.eval_casteljau(t);

            for (idx, axis) in p.into_iter().enumerate() {
                assert!((axis >= (bounds[idx].0 - max_err)) && (axis <= (bounds[idx].1 + max_err)))
            }
        }
    }
}
