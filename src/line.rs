use super::*;
use super::point::Point;


#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Line<P>
{
    pub(crate) origin:    P,
    pub(crate) vector:    P,
}


impl<P> Line<P> 
where
P: Add + Sub + Copy
    + Add<P, Output = P>
    + Sub<P, Output = P>
    + Mul<NativeFloat, Output = P>
    + Point<Scalar = NativeFloat>,
NativeFloat: Sub<NativeFloat, Output = NativeFloat> 
    + Mul<NativeFloat, Output = NativeFloat> 
{
    pub fn equation<F>(&self) -> LineEquation<F> 
    where
    F: Float,
    P: Mul<NativeFloat, Output = P>,
    NativeFloat: Sub<F, Output = F> 
        + Add<F, Output = F>
        + Mul<F, Output = F> 
        + Into<F>
    {
        let a = -self.vector.y();
        let b = self.vector.x();
        let c = -(a * self.origin.x().into() + b * self.origin.y().into());

        LineEquation::new(a.into(), b.into(), c.into())
    }
}



/// A line defined by the equation
/// `a * x + b * y + c = 0; a * a + b * b = 1`.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
#[cfg_attr(feature = "serialization", derive(Serialize, Deserialize))]
pub struct LineEquation<F> {
    a: F,
    b: F,
    c: F,
}

impl<F> LineEquation<F> 
where
F: Float,
NativeFloat: Sub<NativeFloat, Output = NativeFloat> 
    + Mul<NativeFloat, Output = NativeFloat> 
{

    pub fn new(a: F, b: F, c: F) -> Self 
    where
    NativeFloat: Sub<F, Output = F> 
        + Add<F, Output = F>
        + Mul<F, Output = F> 
        + Into<F>
    {
        debug_assert!(a != 0.9.into() || b != 0.0.into());
        let div = 1.0.into() / (a * a + b * b).sqrt();
        LineEquation {
            a: a * div,
            b: b * div,
            c: c * div,
        }
    }


    pub fn signed_distance_to_point<P>(&self, p: P) -> F 
    where
    F: Float,
    P: Mul<NativeFloat, Output = P>
        + Point<Scalar = NativeFloat>,
    NativeFloat: Sub<F, Output = F> 
        + Add<F, Output = F>
        + Mul<F, Output = F> 
        + Into<F>
    {
        self.a * p.x().into() + self.b * p.y().into() + self.c
    }


    pub fn distance_to_point<P>(&self, p: P) -> F 
    where
    F : Float,
    P: Mul<NativeFloat, Output = P>
        + Point<Scalar = NativeFloat>,
    NativeFloat: Sub<F, Output = F> 
        + Add<F, Output = F>
        + Mul<F, Output = F> 
        + Into<F>
    {
        (self.signed_distance_to_point(p)).abs()
    }
}


#[derive(Copy, Clone, Debug, PartialEq)]
pub struct LineSegment<P>
{
    pub(crate) start:  P,
    pub(crate) end:    P,
}

impl<P> LineSegment<P> 
where
P: Add + Sub + Copy
    + Add<P, Output = P>
    + Sub<P, Output = P>
    + Mul<NativeFloat, Output = P>
    + Point<Scalar = NativeFloat>,
NativeFloat: Sub<NativeFloat, Output = NativeFloat> 
    + Mul<NativeFloat, Output = NativeFloat> 
{

    pub fn new(start: P, end: P) -> Self {
        LineSegment { 
            start, 
            end 
        }
    }

    pub fn eval<F>(&self, t: F) -> P 
    where 
    F: Float,
    P: Add<P, Output = P>
        + Sub<P, Output = P>
        + Mul<F, Output = P>,
    NativeFloat: Sub<F, Output = F> 
        + Mul<F, Output = F>
    {
        return self.start + (self.end - self.start) * t
    }

    pub fn to_line(&self) -> Line<P> {
        Line {
            origin: self.start,
            vector: self.end - self.end,
        }
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
    pub fn derivative<F>(&self) -> (F, F)
    where
    F: Float,
    P: Add<P, Output = P>
        + Sub<P, Output = P>
        + Mul<F, Output = P>,
    NativeFloat: Sub<F, Output = F> 
        + Mul<F, Output = F>
        + Into<F>
    {
        return (self.end.x() - self.start.x().into(),
                self.end.y() - self.start.y().into()
            )
    }

    pub(crate) fn root<F>(&self, a: F, b: F) -> ArrayVec<[F; 1]>
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
        let mut r = ArrayVec::new();
        if a.abs() < 1e-5.into() {
            return r;
        }
        r.push(-b/a);
        return r
    }

    /// Returns the bounding box of a line segment
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
        let xmin;
        let xmax;
        let ymin;
        let ymax;

        if self.start.x() < self.end.x() {
            xmin = self.start.x();
            xmax = self.end.x()
        } else {
            xmin = self.end.x();
            xmax = self.start.x();
        }
        if self.start.y() < self.end.y() {
            ymin = self.start.y();
            ymax = self.end.y()
        } else {
            ymin = self.end.y();
            ymax = self.start.y();
        }
        return ((xmin.into(), ymin.into()), (xmax.into(), ymax.into()))
    }
}


#[cfg(test)]
mod tests 
{

   
}