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
        let mut p: [P; {N}] = self.control_points;
        // loop up to degree = control_points.len() -1
        for i in 1..p.len() {
            for j in 0..p.len() - i {
                p[j] = p[j] * (1.0 - t) + p[j+1] * t;
            }
        }
        p[0]
    }
}

#[cfg(test)]
mod tests 
{
    use super::*;
    use super::point2::Point2;
    //use crate::num_traits::{Pow};
    #[test]
    fn construct_and_eval() {
        let points = [
                Point2::new(0f64,  1.77f64),
                Point2::new(1.1f64, -1f64),
                Point2::new(4.3f64,3f64),
                Point2::new(3.2f64, -4f64),
                Point2::new(7.3f64, 2.7f64),
                Point2::new(8.9f64, 1.7f64)];
        // try to initialize an object
        let curve: Bezier<Point2<f64>, 6> = Bezier::new(points);

        let nsteps: usize = 100;                                
        for t in 0 ..= nsteps {
            let t = (t as f64) * 1f64/(nsteps as f64);
            curve.eval(t);
        }
    }
}