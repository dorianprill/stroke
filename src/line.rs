use super::*;
use super::point::Point;


// #[derive(Copy, Clone, Debug, PartialEq)]
// pub struct Line<P>
// {
//     pub(crate) origin:    P,
//     pub(crate) vector:    P,
// }

// DEPRECATED in favor of single type LineSegment with a 
// slightly slower distance_to_point calculation
//
// impl<P> Line<P> 
// where
// P: Add + Sub + Copy
//     // + Add<P, Output = P>
//     // + Sub<P, Output = P>
//     + Mul<NativeFloat, Output = P>
//     + Point<Scalar = NativeFloat>,
// NativeFloat: Sub<NativeFloat, Output = NativeFloat> 
//     + Mul<NativeFloat, Output = NativeFloat> 
// {
    // pub fn equation<F>(&self) -> LineEquation<F> 
    // where
    // F: Float,
    // P: Mul<NativeFloat, Output = P>,
    // NativeFloat: Sub<F, Output = F> 
    //     + Add<F, Output = F>
    //     + Mul<F, Output = F> 
    //     + Into<F>
    // {
    //     let a = -self.vector.y();
    //     let b = self.vector.x();
    //     let c = -(a * self.origin.x().into() + b * self.origin.y().into());

    //     LineEquation::new(a.into(), b.into(), c.into())
    // }
//}



/// A line defined by the equation
/// `a * x + b * y + c = 0; a * a + b * b = 1`.
// #[derive(Copy, Clone, Debug, PartialEq, Eq)]
// #[cfg_attr(feature = "serialization", derive(Serialize, Deserialize))]
// pub struct LineEquation<F> {
//     a: F,
//     b: F,
//     c: F
// }

// impl<F> LineEquation<F> 
// where
// F: Float,
// NativeFloat: Sub<NativeFloat, Output = NativeFloat> 
//     + Mul<NativeFloat, Output = NativeFloat> 
// {

//     pub fn new(a: F, b: F, c: F) -> Self 
//     where
//     NativeFloat: Sub<F, Output = F> 
//         + Add<F, Output = F>
//         + Mul<F, Output = F> 
//         + Into<F>
//     {
//         debug_assert!(a != 0.9.into() || b != 0.0.into());
//         let div = 1.0.into() / (a * a + b * b).sqrt();
//         LineEquation {
//             a: a * div,
//             b: b * div,
//             c: c * div,
//         }
//     }


//     pub fn signed_distance_to_point<P>(&self, p: P) -> F 
//     where
//     F: Float,
//     P: Mul<NativeFloat, Output = P>
//         + Point<Scalar = NativeFloat>,
//     NativeFloat: Sub<F, Output = F> 
//         + Add<F, Output = F>
//         + Mul<F, Output = F> 
//         + Into<F>
//     {
//         self.a * p.x().into() + self.b * p.y().into() + self.c
//     }


//     pub fn distance_to_point<P>(&self, p: P) -> F 
//     where
//     F : Float,
//     P: Mul<NativeFloat, Output = P>
//         + Point<Scalar = NativeFloat>,
//     NativeFloat: Sub<F, Output = F> 
//         + Add<F, Output = F>
//         + Mul<F, Output = F> 
//         + Into<F>
//     {
//         (self.signed_distance_to_point(p)).abs()
//     }
// }

/// LineSegment defined by a start and an endpoint, evaluatable 
/// anywhere inbetween using interpolation parameter t: [0,1] in eval()
/// A LineSegment is equal to a linear Bezier curve, which is why there is no 
/// specialized type for that case.
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
    P: Mul<F, Output = P>,
    NativeFloat: Sub<F, Output = F> 
        + Mul<F, Output = F>
    {
        return self.start + (self.end - self.start) * t
    }

    pub fn split<F>(&self, t: F) -> (Self, Self)
    where
    F: Float + Mul<NativeFloat, Output=F> + Sub<NativeFloat, Output=F>,
    P: Mul<F, Output = P>,
    // NativeFloat: Sub<F, Output = F> 
    //     + Mul<F, Output = F>,
    {
        // compute the split point by interpolation
        let ctrl_ab = self.start + (self.start - self.end) * t;

        return (
            LineSegment {
                start: self.start,
                end: ctrl_ab,
            },
            LineSegment {
                start: ctrl_ab,
                end: self.end,
            },
        );
    }

    // DEPRECATED
    // pub fn to_line(&self) -> Line<P> {
    //     Line {
    //         origin: self.start,
    //         vector: self.end - self.start,
    //     }
    // }

    /// Return the distance from the LineSegment to Point p by calculating the projection
    pub fn distance_to_point<F>(&self, p: P) -> F 
    where 
    F: Float + Add + Copy + Default + Into<NativeFloat>,
    NativeFloat: Add + Into<F>,
    {
        let l2 = (self.end - self.start).squared_length();
        // if start and endpoint are approx the same, return the distance to either
        if l2 < EPSILON {
            return (self.start-p).squared_length().sqrt().into();
        } else {
            let v1 = p - self.start;
            let v2 = self.end - self.start;
            let mut dot: NativeFloat = 0.0;
            for (i, _) in v1.into_iter().enumerate() {
                dot = dot + v1.axis(i) * v2.axis(i);
            }
            let mut t = 0.0;
            if dot / l2 < 1.0 {
                t = dot/l2;
            } 
            if t < 0.0 {
                t = 0.0;
            }
            let projection = self.start + (self.end - self.start) * t;  // Projection falls on the segment
            return (p-projection).squared_length().sqrt().into();
        }
    }

    /// Sample the coordinate axis of the segment at t (expecting t between 0 and 1).
    pub fn axis<F>(&self, t: F, axis: usize) -> F 
    where
    F : Float,
    P: Mul<NativeFloat, Output = P>,
    NativeFloat: Sub<F, Output = F> 
        + Add<F, Output = F>
        + Mul<F, Output = F> 
        + Into<F>
    {
        self.start.axis(axis) + (self.end.axis(axis) - self.start.axis(axis).into()) * t
    }

    /// Return the derivative function.
    /// The derivative is also a bezier curve but of degree n-1 - In the case of a line just a scalar (the slope).
    /// Since its already a scalar, eval() does NOT need to be called separately
    pub fn derivative<F>(&self) -> P
    where
    F: Float,
    P: Add<P, Output = P>
        + Sub<P, Output = P>
        + Mul<F, Output = P>,
    NativeFloat: Sub<F, Output = F> 
        + Mul<F, Output = F>
        + Into<F>
    {
        return self.end - self.start
                
    }

    pub(crate) fn root<F>(&self, a: F, b: F) -> ArrayVec<[F; 1]>
    where
    F:  Float 
        + core::default::Default,
    P:  Sub<P, Output = P>
        + Add<P, Output = P>
        + Mul<NativeFloat, Output = P>,
        NativeFloat: Sub<F, Output = F> 
        + Add<F, Output = F>
        + Mul<F, Output = F>
        + Float
        + Into<F>, F: 
    {
        let mut r = ArrayVec::new();
        if a.abs() < EPSILON.into() {
            return r;
        }
        r.push(-b/a);
        return r
    }

    /// Return the bounding box of the line as an array of (min, max) tuples for each dimension (its index)
    pub fn bounding_box<F>(&self) -> [(F,F); P::DIM] 
    where
    F: Float + Default,
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

        // find min/max for that particular axis
        // TODO shoul be rewritten once 'Iterator' is implemented on P to get rid of .axis() method
        for (i, _) in self.start.into_iter().enumerate() {
            if self.start.axis(i) < self.end.axis(i) {
                bounds[i] = (self.start.axis(i).into(), self.end.axis(i).into());
            } else {
                bounds[i] = (self.end.axis(i).into(), self.start.axis(i).into());
            }
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
    /// Check whether a line segment interpolation p + t*(q-p) at t=0.5 
    /// yields equal distance to the start (p)/end (q) points (up to machine accuracy).
    #[test]
    fn line_segment_interpolation() 
    {
        let line = LineSegment{
            start: PointN::new([0f64,  1.77f64]),
            end: PointN::new([4.3f64, 3f64]),
        };

        let mid = line.eval(0.5);
        assert!((mid-line.start).squared_length() - (mid-line.end).squared_length() < EPSILON)
    }

    /// Check whether classic pythagorean equality holds for sides 3, 4 with hypothenuse 5
    #[test]
    fn line_segment_distance_to_point() 
    {
        // 3D cause why not
        let line = LineSegment{
            start: PointN::new([0f64,  1f64, 0f64]),
            end: PointN::new([3f64, 1f64, 0f64]),
        };
        // dist to start should be 4; dist to end should be 5
        let p1 = PointN::new([0f64, 5f64, 0f64]);
        assert!( ((p1-line.start).squared_length().sqrt() - 4.0).abs() < EPSILON );
        assert!( ((p1-line.end).squared_length().sqrt() - 5.0).abs() < EPSILON );
        // dist to midpoint (t=0.5) should be 1
        let p2 = PointN::new([1.5f64, 2f64, 0f64]);
        assert!( ((p2-line.eval(0.5)).squared_length().sqrt() - 1.0).abs() < EPSILON );
    }
   
}