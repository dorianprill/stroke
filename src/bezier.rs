use super::*;
use super::point::Point;

/// General implementation of a Bezier curve of arbitrary degree.
/// The curve is solely defined by an array of 'control_points'. The degree is defined as degree = control_points.len() - 1.
/// Points on the curve can be evaluated with an interpolation parameter 't' in interval [0,1] using the eval() and eval_casteljau() methods.
/// Generic parameters:
/// P: Generic points 'P' as defined by there Point trait
/// const generic parameters:
/// N: Number of control points
#[derive(Clone, Copy)]
pub struct Bezier<P, const N: usize> 
where 
P: Point + Copy,
{
    /// Control points which define the curve and hence its degree
    control_points: [P; N],
}

impl<P, const N: usize> Bezier<P, {N}> 
where
P: Add + Sub + Copy
    + Add<P, Output = P>
    + Sub<P, Output = P>
    + Mul<NativeFloat, Output = P>
    + Point<Scalar = NativeFloat>,
{
    /// Create a new Bezier curve that interpolates the `control_points`. The degree is defined as degree = control_points.len() - 1.
    /// Desired curve must have a valid number of control points and knots in relation to its degree or the constructor will return None. 
    /// A B-Spline curve requires at least one more control point than the degree (`control_points.len() >
    /// degree`) and the number of knots should be equal to `control_points.len() + degree + 1`.
    pub fn new(control_points: [P; N]) -> Bezier<P, {N}> {
        Bezier{
            control_points
        }
    }


    /// Evaluate a point on the curve at point 't' which should be in the interval [0,1]
    /// This is implemented using De Casteljau's algorithm (over a temporary array with const generic sizing)
    pub fn eval<F>(&self, t: F) -> P 
    where
    F: Float,
    P: Add<P, Output = P>
        + Sub<P, Output = P>
        + Mul<F, Output = P>,
    NativeFloat: Sub<F, Output = F> 
        + Mul<F, Output = F>
        + Into<F>
    {
        // start with a copy of the original control points array and succesively use it for evaluation
        let mut p: [P; N] = self.control_points;
        // loop up to degree = control_points.len() -1
        for i in 1..=p.len() {
            for j in 0..p.len() - i {
                p[j] = p[j] * (1.0 - t) + p[j+1] * t;
            }
        }
        p[0]
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
        // start with a copy of the original control points for now
        // TODO how to initialize const generic array without using unsafe?
        let mut left: [P; N] = self.control_points.clone();
        let mut right: [P; N] = self.control_points.clone();
        // these points get overriden each iteration; we save the intermediate results to 'left' and 'right'
        let mut casteljau_points: [P; N] = self.control_points.clone();

        for i in 1..=casteljau_points.len() 
        {
            // save start point of level
            left[i-1] = casteljau_points[0];
            // save end point of level
            right[right.len()-i] = casteljau_points[right.len()-i];
            // calculate next level of points (one less point each level until we reach one point, the one at t)
            for j in 0..casteljau_points.len() - i {
                casteljau_points[j] = casteljau_points[j] * (1.0-t) + casteljau_points[j+1] * t; 
            }
        }
        return ( Bezier{ control_points: left }, Bezier{ control_points: right })
    }

    // TODO expressions in const generics are not yet supported
    // pub fn derivative<F>(&self) -> Bezier<P, {N-1}>
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
    //     let new_points: [P; {N-1}]; 
    //     return Bezier::new(new_points)
    // }

}

#[cfg(test)]
mod tests 
{
    use super::*;
    use super::point2::Point2;
    //use crate::num_traits::{Pow};
    #[test]
    fn eval_endpoints() {
        let points = [
                Point2::new(0f64,  1.77f64),
                Point2::new(1.1f64, -1f64),
                Point2::new(4.3f64,3f64),
                Point2::new(3.2f64, -4f64),
                Point2::new(7.3f64, 2.7f64),
                Point2::new(8.9f64, 1.7f64)];
        // try to initialize an object
        let curve: Bezier<Point2<f64>, 6> = Bezier::new(points);

        // let nsteps: usize = 100;                                
        // for t in 0 ..= nsteps {
        //     let t = (t as f64) * 1f64/(nsteps as f64);
        //     curve.eval(t);
        // }
        // check if start/end points match
        let max_err = 1e-14;
        let mut err: Point2<f64>;
        let start = curve.eval(0.0);
        err = start - points[0];
        assert!( (err.x.abs() < max_err) && (err.y.abs() < max_err) );
        let end = curve.eval(1.0);
        err = end - points[points.len() - 1 ];
        assert!( (err.x.abs() < max_err) && (err.y.abs() < max_err) );
    }

    #[test]
    fn split_equivalence() {
        // chose some arbitrary control points and construct a cubic bezier
        let bezier = Bezier{control_points: 
            [Point2{x:0f64,  y:1.77f64},
            Point2{x:2.9f64, y:0f64},
            Point2{x:4.3f64, y:3f64},
            Point2{x:3.2f64, y:-4f64}]
        };
        // split it at an arbitrary point
        let at = 0.5;
        let (left, right) = bezier.split(at);
        // compare left and right subcurves with parent curve
        // this is tricky as we have to map t->t/2 (for left) which will 
        // inevitably contain rounding errors from floating point ops.
        // instead, take the difference of the two points which must not exceed the absolute error
        // TODO update test to use norm() instead, once implemented for Point
        let max_err = 1e-14;
        let nsteps: usize =  1000;                                      
        for t in 0..=nsteps {
        let t = t as f64 * 1f64/(nsteps as f64);
        // dbg!(t);
        // dbg!(bezier.eval(t/2.0));
        // dbg!(left.eval(t));
        // dbg!(bezier.eval((t*0.5)+0.5));
        // dbg!(right.eval(t));
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
}