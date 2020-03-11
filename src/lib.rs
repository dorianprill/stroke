#![no_std]
use core::ops::{Add, Sub, Mul};

extern crate num_traits;
use num_traits::{float::Float};

#[derive(Debug, Copy, Clone)]
#[allow(dead_code)]
pub struct Point2<T>
{
    x: T,
    y: T,
}

impl<T> Point2<T> 
where T: Add + Add<T,Output=T> + Sub + Mul + Mul<T,Output=T> + Clone {
    /// Creates a new Point2<T>, which requires that 
    /// T implements Add, Sub, Mul, and Clone
    pub fn new(x: T, y: T) -> Self {
        Point2 {
            x: x,
            y,
        }
    }
}

pub trait Distance {
    type ScalarDist;
    fn distance(&self, other: Self) -> Self::ScalarDist;
}

// conditionally compiled newtype pattern used to determine which size float to use in arclen() and for tests
// so that the library can abstract over both 32 and 64 bit architectures
// TODO there has to be a better architectural solution to this... if not, generate impls with macro
#[cfg(target_pointer_width = "64")]
type NativeFloat = f64;
// same thing for 32 bit
#[cfg(target_pointer_width = "32")]
type NativeFloat = f32;



impl Distance for Point2<f64> {
    type ScalarDist = NativeFloat;
    fn distance(&self, other: Self) -> NativeFloat {
        ( ((self.x - other.x) * (self.x - other.x))
            + ((self.y - other.y) * (self.y - other.y)) ) .sqrt()
    }
}

impl Distance for Point2<i64> {
    type ScalarDist = NativeFloat;
    fn distance(&self, other: Self) -> NativeFloat {
        ( ((self.x - other.x) * (self.x - other.x)) as NativeFloat
            + ((self.y - other.y) * (self.y - other.y)) as NativeFloat) .sqrt()
    }
}

impl Distance for Point2<u64> {
    type ScalarDist = NativeFloat;
    fn distance(&self, other: Self) -> NativeFloat {
        ( ((self.x - other.x) * (self.x - other.x)) as NativeFloat
            + ((self.y - other.y) * (self.y - other.y)) as NativeFloat ) .sqrt()
    }
}





impl<T> PartialEq for Point2<T> 
where T: PartialOrd {
    fn eq(&self, other: &Self) -> bool {
        (self.x == other.x) && (self.y == other.y)
    }
}


impl<T> Add for Point2<T>
where
    T: Add<Output=T>,
{
    type Output = Self;

    fn add(self, other: Point2<T>) -> Point2<T> {
        Point2 {
            x: self.x + other.x,
            y: self.y + other.y
        }
    }
}

impl<T> Sub for Point2<T>
where 
    T: Sub<Output=T>
{
    type Output = Self;

    fn sub(self, other: Self) -> Self::Output {
        Point2 {
            x: self.x - other.x,
            y: self.y - other.y,
        }
    }
}

impl<T,U> Mul<U> for Point2<T>
where
    // How you have the mulitplication done is mulitpling T * U => T, this
    // trait bounds for T will specify this requirement as the mul operator is
    // translated to using the first operand as self and the second as rhs. 
    T: Mul<U,Output=T> + Copy,
    U: Clone,
{
    type Output = Point2<T>;

    fn mul(self, _rhs: U) -> Point2<T> {
        return Point2{x: self.x * _rhs.clone(), y: self.y * _rhs}
    }
}


#[allow(dead_code)]
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
    P: Add + Sub + Copy + Distance<ScalarDist = NativeFloat>,
{

    pub fn new(start: P, ctrl1: P, ctrl2: P,  end: P) -> Self {
        CubicBezier { 
            start, 
            ctrl1, 
            ctrl2, 
            end 
        }
    }


    /// Evaluate a CubicBezier curve at t by direct evaluation of the polynomial (not numerically stable)
    fn eval<F>(&self, t: F) -> P 
    where 
    F: Float,
    P: Add<P, Output = P>
        + Sub<P, Output = P>
        + Mul<F, Output = P>,
    NativeFloat: Sub<F, Output = F> + Mul<F, Output = F>
    {
        return self.start * ((1.0-t) * (1.0-t) * (1.0-t))
                + self.ctrl1 * (3.0 * t * (1.0-t) * (1.0-t))
                + self.ctrl2 * (3.0 * t * t * (1.0-t))
                + self.end * (t * t * t);
    }

    /// Evaluate a CubicBezier curve at t using the numerically stable De Casteljau algorithm
    fn eval_casteljau<F>(&self, t: F) -> P 
    where 
    F: Float,
    P: Add<P, Output = P>
        + Sub<P, Output = P>
        + Mul<F, Output = P>,
    NativeFloat: Sub<F, Output = F> + Mul<F, Output = F>
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

        return ctrl_3ab
    }

    /// Approximates the arc length of the curve by flattening it with straight line segments.
    /// This works quite well, at ~32 segments it should already provide an error < 0.5
    /// Remember arclen also works by linear approximation, not the integral, so we have to accept error!
    /// This approximation is unfeasable if desired accuracy is greater than 2 decimal places
    fn arclen<F>(&self, nsteps: usize) -> NativeFloat
    where 
    F: Float,
    P: Add<P, Output = P>
        + Sub<P, Output = P>
        + Mul<F, Output = P>,
    NativeFloat: Sub<F, Output = F> + Mul<F, Output = F> + Into<F>
    {
        let stepsize: NativeFloat = 1.0/(nsteps as NativeFloat);
        let mut arclen: NativeFloat = 0.0;
        for t in 1..nsteps {
            let t = t as NativeFloat * 1.0.into()/(nsteps as NativeFloat).into();
            let p1 = self.eval_casteljau(t);
            let p2 = self.eval_casteljau(t+stepsize.into());

            arclen = arclen + p1.distance(p2);
        
        }
        return arclen
    }


    fn split<F>(&self, t: F) -> (Self, Self)
    where
    F: Float,
    P:  Sub<P, Output = P>
        + Add<P, Output = P>
        + Mul<F, Output = P>,
    NativeFloat: Sub<F, Output = F> + Mul<F, Output = F>,
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

}


#[cfg(test)]
mod tests 
{
    use super::*;
    use crate::num_traits::{Pow};

    #[test]
    fn eval_equivalence() {
        // all eval methods should be approximately equivalent for well defined test cases
        // and not equivalent where numerical stability becomes an issue for normal eval
        let bezier = CubicBezier{ start:  Point2{x:0f64,  y:1.77f64},
                                  ctrl1: Point2{x:2.9f64, y:0f64},
                                  ctrl2: Point2{x:4.3f64, y:-3f64},
                                  end:   Point2{x:3.2f64,  y:4f64}};

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
        let bezier = CubicBezier{ start:  Point2{x:0f64,  y:1.77f64},
                                  ctrl1: Point2{x:2.9f64, y:0f64},
                                  ctrl2: Point2{x:4.3f64, y:-3f64},
                                  end:   Point2{x:3.2f64,  y:4f64}};
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

}
