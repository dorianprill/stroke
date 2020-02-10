use core::ops::{Add, Sub, Mul};

extern crate num_traits;
use num_traits::{float::Float};

#[derive(Debug, Copy, Clone)]
#[allow(dead_code)]
pub struct Point2D<T>
{
    x: T,
    y: T,
}

impl<T> Point2D<T> 
where T: Add + Sub + Mul + Clone {
    // Let this be the only way users can create Point2D, which requires that 
    // T implement Add, Sub, Mul, and Clone
    pub fn new(x: T, y: T) -> Self {
        Point2D {
            x: x,
            y, // Shortcut (I think this is the preferred) way.
        }
    }
}


impl<T> Add for Point2D<T>
where
    T: Add<Output=T>,
{
    type Output = Self;

    fn add(self, other: Point2D<T>) -> Point2D<T> {
        Point2D {
            x: self.x + other.x,
            y: self.y + other.y
        }
    }
}

impl<T> Sub for Point2D<T>
where 
    T: Sub<Output=T>
{
    type Output = Self;

    fn sub(self, other: Self) -> Self::Output {
        Point2D {
            x: self.x - other.x,
            y: self.y - other.y,
        }
    }
}

impl<T,U> Mul<U> for Point2D<T>
where
    // How you have the mulitplication done is mulitpling T * U => T, this
    // trait bounds for T will specify this requirement as the mul operator is
    // translated to using the first operand as self and the second as rhs. 
    T: Mul<U,Output=T>,
    U: Clone,
{
    type Output = Point2D<T>;

    fn mul(self, _rhs: U) -> Point2D<T> {
        return Point2D{x: self.x * _rhs.clone(), y: self.y * _rhs}
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
impl<T> CubicBezier<T> 
where 
    T: Add + Sub + Mul + Copy,
{

    pub fn new(start: T, ctrl1: T, ctrl2: T,  end: T) -> Self {
        CubicBezier { 
            start, 
            ctrl1, 
            ctrl2, 
            end 
        }
    }

    fn eval<F>(&mut self, t: F) -> T 
    where 
        F: Float,
        T: Add<T, Output=T> + Add<F, Output=T> + Sub<F, Output=T> + Mul<F, Output=T>,
        f64: Sub<F, Output=F> + Mul<F, Output=F>, // this is the primitive type f64
    {
        return self.start * ((1.0-t) * (1.0-t) * (1.0-t))
                + self.ctrl1 * (3.0 * t * (1.0-t) * (1.0-t))
                + self.ctrl2 * (3.0 * t * t * (1.0-t))
                + self.end * (t * t * t);
    }


    fn split<F>(&mut self, t: F) -> (Self, Self)
    where
        F: Float,
        T: Add<F, Output = T>
            + Sub<T, Output = T>
            + Add<T, Output = T>
            + Sub<F, Output = T>
            + Mul<F, Output = T>,
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
        let circle = |p: Point2D<f64>| -> f64 { ( p.x.pow(2) as f64 
                                                + p.y.pow(2) as f64)
                                                .sqrt() - 1f64};

        // define control points for 4 bezier segments 
        // control points are chosen for minimum radial distance error
        // according to: http://spencermortensen.com/articles/bezier-circle/ 
        // TODO don't hardcode values
        let c               = 0.551915024494;
        let max_drift_perc  = 0.019608; // radial drift percent
        let max_error       = max_drift_perc * 0.01; // absolute max radial error

        let cubic_bezier_circle = CubicBezier{ start:  Point2D{x:0f64,  y:1f64},
                                                ctrl1: Point2D{x:c,     y:1f64},
                                                ctrl2: Point2D{x:1f64,  y:c},
                                                end:   Point2D{x:1f64,  y:0f64}};
        let nsteps =  1000;                                      
        for t in -nsteps..nsteps {
            let t = t as f64 * 1f64/(nsteps as f64);
            let point = cubic_bezier_circle.eval(t);
            let contour = circle(point);
            assert!( contour.abs() <= max_error );
        }
    }
}
