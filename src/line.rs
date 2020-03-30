use super::*;
use super::point2::*;
#[allow(unused_imports)]
use super::cubic_bezier::CubicBezier;

#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Line<P>
{
    pub(crate) start:  P,
    pub(crate) end:    P,
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
    + Mul<NativeFloat, Output = NativeFloat> 
{

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


#[cfg(test)]
mod tests 
{

   
}