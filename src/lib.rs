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


pub trait Coordinate {
    type Coordinate;
    fn x(&self) -> Self::Coordinate;
    fn y(&self) -> Self::Coordinate;
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
#[cfg(target_pointer_width = "64")]
type NativeInt   = i64;
#[cfg(target_pointer_width = "64")]
type NativeUInt  = u64;

// same thing for 32 bit
#[cfg(target_pointer_width = "32")]
type NativeFloat = f32;
#[cfg(target_pointer_width = "32")]
type NativeInt   = i32;
#[cfg(target_pointer_width = "32")]
type NativeUInt  = u32;


impl Distance for Point2<NativeFloat> {
    type ScalarDist = NativeFloat;
    fn distance(&self, other: Self) -> NativeFloat {
        ( ((self.x - other.x) * (self.x - other.x))
            + ((self.y - other.y) * (self.y - other.y)) ) .sqrt()
    }
}

impl Distance for Point2<NativeInt> {
    type ScalarDist = NativeFloat;
    fn distance(&self, other: Self) -> NativeFloat {
        ( ((self.x - other.x) * (self.x - other.x)) as NativeFloat
            + ((self.y - other.y) * (self.y - other.y)) as NativeFloat) .sqrt()
    }
}

impl Distance for Point2<NativeUInt> {
    type ScalarDist = NativeFloat;
    fn distance(&self, other: Self) -> NativeFloat {
        ( ((self.x - other.x) * (self.x - other.x)) as NativeFloat
            + ((self.y - other.y) * (self.y - other.y)) as NativeFloat ) .sqrt()
    }
}


impl Coordinate for Point2<NativeFloat> {
    type Coordinate = NativeFloat;
    fn x(&self) -> NativeFloat {
        self.x
    }
    fn y(&self) -> NativeFloat {
        self.y
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



#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Line<P>
{
    start:  P,
    end:    P,
}

impl<P> Line<P> 
where
P: Add + Sub + Copy
    + Add<P, Output = P>
    + Sub<P, Output = P>
    + Mul<NativeFloat, Output = P>
    + Distance<ScalarDist = NativeFloat> 
    + Coordinate<Coordinate = NativeFloat>,
NativeFloat: Sub<NativeFloat, Output = NativeFloat> 
    + Mul<NativeFloat, Output = NativeFloat> {

    pub fn new(start: P, end: P) -> Self {
        Line { 
            start, 
            end 
        }
    }

    pub fn eval(&self, t: NativeFloat) -> P {
        return self.start + (self.end - self.start) * t
    }

    /// Sample the x coordinate of the segment at t (expecting t between 0 and 1).
    pub fn x<F>(&self, t: F) -> F 
    where
    F : Float,
    P: Mul<NativeFloat, Output = P>,
    NativeFloat: Sub<F, Output = F> 
        + Add<F, Output = F>
        + Mul<F, Output = F> 
        + Into<F>
    {
        self.start.x() + (self.end.x() - self.start.x().into()) * t
    }

    /// Sample the y coordinate of the segment at t (expecting t between 0 and 1).
    pub fn y<F>(&self, t: F) -> F 
    where
    F : Float,
    P: Mul<NativeFloat, Output = P>,
    NativeFloat: Sub<F, Output = F> 
        + Add<F, Output = F>
        + Mul<F, Output = F> 
        + Into<F>
    {
        self.start.y() + (self.end.y() - self.start.y().into()) * t
    }

    /// Return the derivative function.
    /// The derivative is also a bezier curve but of degree n-1 - In the case of a line just a scalar (the slope).
    /// Since its already a scalar, eval() does NOT need to be called separately
    pub fn derivative(&self) -> Point2<NativeFloat>
    where
    P:  Sub<P, Output = P>
        + Add<P, Output = P>
        + Mul<NativeFloat, Output = P>,
    NativeFloat: Sub<NativeFloat, Output = NativeFloat> 
        + Mul<NativeFloat, Output = NativeFloat>
        + Float
    {
        return Point2{
            x: self.end.x() - self.start.x(),
            y: self.end.y() - self.start.y()
        }
    }

    pub fn roots(&self) -> NativeFloat
    where
    P:  Sub<P, Output = P>
        + Add<P, Output = P>
        + Mul<NativeFloat, Output = P>,
    NativeFloat: Sub<NativeFloat, Output = NativeFloat> 
        + Mul<NativeFloat, Output = NativeFloat>
    {
        let slope = (self.end.y() - self.start.y()) / (self.end.x() - self.end.y());
        let intercept = self.start.y() - slope * self.start.x();
        return -intercept/slope
    }
}



#[derive(Copy, Clone, Debug, PartialEq)]
pub struct QuadraticBezier<P>
{
    start:  P,
    ctrl:   P,
    end:    P,
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
        let ctrl_1bc   = self.ctrl1 + (self.ctrl2 - self.ctrl1) * t;
        let ctrl_1cd   = self.ctrl2 + (self.end - self.ctrl2) * t;
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
    fn arclen<F>(&self, nsteps: usize) -> F
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


    fn split<F>(&self, t: F) -> (Self, Self)
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


}




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
