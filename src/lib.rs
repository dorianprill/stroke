use core::ops::{Add, Sub, Mul};

extern crate num_traits;
use num_traits::{float::Float, Pow};

#[derive(Debug, Copy, Clone)]
#[allow(dead_code)]
pub struct Point2<T>
{
    x: T,
    y: T,
}

impl<T> Point2<T> 
where T: Add + Sub + Mul + Clone + Mul<T,Output=T> {
    // Let this be the only way users can create Point2, which requires that 
    // T implement Add, Sub, Mul, and Clone
    pub fn new(x: T, y: T) -> Self {
        Point2 {
            x: x,
            y,
        }
    }

    // TODO fix norm() trait bound issue
    // pub fn norm(&self) -> T {
    //     ((self.x * self.x)
    //         + (self.y * self.y)).sqrt()
    // }
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
pub struct CubicBezier<T>
{
    start:  T,
    ctrl1:  T,
    ctrl2:  T,
    end:    T,
}


#[allow(dead_code)]
impl<P> CubicBezier<P> 
where 
    P: Add + Sub + Copy,
{

    pub fn new(start: P, ctrl1: P, ctrl2: P,  end: P) -> Self {
        CubicBezier { 
            start, 
            ctrl1, 
            ctrl2, 
            end 
        }
    }

    fn eval<F>(&self, t: F) -> P 
    where 
    F: Float,
    P: Add<P, Output = P>
        + Sub<P, Output = P>
        + Mul<F, Output = P>,
    f64: Sub<F, Output = F> + Mul<F, Output = F>
    {
        return self.start * ((1.0-t) * (1.0-t) * (1.0-t))
                + self.ctrl1 * (3.0 * t * (1.0-t) * (1.0-t))
                + self.ctrl2 * (3.0 * t * t * (1.0-t))
                + self.end * (t * t * t);
    }


    fn split<F>(&mut self, t: F) -> (Self, Self)
    where
        F: Float,
        P:  Sub<P, Output = P>
            + Add<P, Output = P>
            + Mul<F, Output = P>,
        f64: Sub<F, Output = F> + Mul<F, Output = F>,
    {
        let ctrl1a   = self.start + (self.ctrl1 - self.start) * t;
        let ctrl2a   = self.ctrl1 + (self.ctrl2 - self.ctrl1) * t;
        let ctrl1aa  = ctrl1a + (ctrl2a - ctrl1a) * t;
        let ctrl3a   = self.ctrl2 + (self.end - self.ctrl2) * t;
        let ctrl2aa  = ctrl2a + (ctrl3a - ctrl2a) * t;
        let ctrl1aaa = ctrl1aa + (ctrl2aa - ctrl1aa) * t;

        return (
            CubicBezier {
                start: self.start,
                ctrl1: ctrl1a,
                ctrl2: ctrl1aa,
                end: ctrl1aaa,
            },
            CubicBezier {
                start: ctrl1aaa,
                ctrl1: ctrl2aa,
                ctrl2: ctrl3a,
                end: self.end,
            },
        );
    }

}


#[cfg(test)]
mod tests 
{
    use super::*;
    use std::f64;
    use crate::num_traits::Pow;

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
            let mut point = bezier_quadrant_1.eval(t);
            let mut contour = circle(point);
            assert!( contour.abs() <= max_error );
            point = bezier_quadrant_2.eval(t);
            contour = circle(point);
            assert!( contour.abs() <= max_error );
            point = bezier_quadrant_3.eval(t);
            contour = circle(point);
            assert!( contour.abs() <= max_error );
            point = bezier_quadrant_4.eval(t);
            contour = circle(point);
            assert!( contour.abs() <= max_error );
        }
    }

    #[test]
    fn split_equivalence() {
        // chose some arbitrary control points and construct a cubic bezier
        let mut bezier = CubicBezier{ start:  Point2{x:0f64,  y:1.77f64},
                                  ctrl1: Point2{x:2.9f64, y:0f64},
                                  ctrl2: Point2{x:4.3f64, y:-3f64},
                                  end:   Point2{x:3.2f64,  y:4f64}};
        // split it at an arbitrary point
        let at = 0.5;
        let (left, right) = bezier.split(at);
        // compare left and right subcurves with parent curve
        // this is tricky as we have to map t->t/2 (for left) which will 
        // inevitably contain rounding errors from floating point ops
        // instead, take the difference of the two points does not exceed an absolute error
        // TODO update test to use norm() instead, once implemented for Point (maybe as trait?)
        let max_err = 1e-14;
        let nsteps =  1000;                                      
        for t in 0..nsteps {
            let t = t as f64 * 1f64/(nsteps as f64);
            dbg!(bezier.eval(t/2.0));
            dbg!(left.eval(t));
            // left
            let mut err = bezier.eval(t/2.0) - left.eval(t);
            dbg!(err);
            assert!( (err.x.abs() < max_err) && (err.y.abs() < max_err) );
            // right
            err = bezier.eval((t*0.5)+0.5) - right.eval(t);
            dbg!(err);
            assert!( (err.x.abs() < max_err) && (err.y.abs() < max_err) );  
        }
    }
}
