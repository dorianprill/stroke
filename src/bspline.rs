use core::slice::*;

use super::point::Point;
use super::*;

/// General Implementation of a BSpline with choosable degree, control points and knots.
/// Generic parameters:
/// P: const generic points array 'P' as defined by the Point trait
/// F: Any float value used for the knots and interpolation (usually the same as the internal generic parameter within P<F>).
/// const generic parameters:
/// C: Number of control points
/// K: Number of Knots
/// D: Degree of the piecewise function used for interpolation degree = order - 1
/// While C, K, O relate to each other in the following manner
///     K = C + O where O = D + 1
#[derive(Clone)]
pub struct BSpline<P, const K: usize, const C: usize, const D: usize>
where
    P: Point,
{
    /// Knot vector
    knots: [P::Scalar; K],
    /// Control points
    control_points: [P; C],
}

impl<P, const K: usize, const C: usize, const D: usize> BSpline<P, { K }, { C }, { D }>
where
    P: Point,
{
    /// Create a new B-spline curve that interpolates
    /// the `control_points` using a piecewise polynomial of `degree` within intervals specified by the `knots`.
    /// The knots _must_ be sorted in non-decreasing order, the constructor enforces this which may yield undesired results.
    /// The degree is defined as `curve_order - 1`.
    /// Desired curve must have a valid number of control points and knots in relation to its degree or the constructor will return None.
    /// A B-Spline curve requires at least one more control point than the degree (`control_points.len() >
    /// degree`) and the number of knots should be equal to `control_points.len() + degree + 1`.
    #[allow(clippy::if_same_then_else)] // allow until a proper Error type is defined for bsplines
    pub fn new(
        knots: [P::Scalar; K],
        control_points: [P; C],
        //degree: usize,
    ) -> Option<BSpline<P, { K }, { C }, { D }>> {
        if control_points.len() <= {D} {
            //panic!("Too few control points for curve");
            None
        } else if knots.len() != control_points.len() + {D} + 1 {
            // panic!(format!("Invalid number of knots, got {}, expected {}", knots.len(),
            //     control_points.len() + degree + 1));
            None
        } else {
            // TODO force sorting of the knots required for binary search (knot span) -> mutable reference required
            // FIX maybe dont sort and just use linear search for knot span, as knot vectors wont be large anyway
            //.sort_by(|a, b| a.partial_cmp(b).unwrap());
            Some(BSpline {
                control_points,
                knots,
            })
        }
    }

    /// Compute a point on the curve at `t` using iterative de boor algorithm. 
    /// The parameter **must** be in the inclusive range of values returned 
    /// by `knot_domain`. If `t` is out of bounds this function will assert
    /// on debug builds and on release builds you'll likely get an out of bounds crash.
    pub fn eval(&self, t: P::Scalar) -> P 
    where [(); D+1]: Sized {
        debug_assert!(t >= self.knot_domain().0 && t <= self.knot_domain().1);
        // Find the knot span that contains t i.e. the first index with a knot value greater than the t we're searching for.
        // We need to find the start of the knot span t is in, such that: knots[span] <= t < knots[span + 1]
        // Note: A custom function is used to exploit binary search (knots are sorted)
        let span = match self.upper_bounds(&self.knots[..], t) {
            Some(x) if x == 0 => {D}, // degree
            Some(x) if x >= self.knots.len() - {D} - 1 => {
                self.knots.len() - {D} - 1
            }
            Some(x) => x,
            None => self.knots.len() - {D} - 1,
        };
        self.de_boor_iterative(t, span)
    }

    /// Returns an iterator over the control points.
    pub fn control_points(&self) -> Iter<'_, P> {
        self.control_points.iter()
    }

    /// Returns an iterator over the knots.
    pub fn knots(&self) -> Iter<'_, P::Scalar> {
        self.knots.iter()
    }

    /// Get the min and max knot domain values for finding the `t` range to compute
    /// the curve over. The curve is only defined over the inclusive range `[min, max]`,
    /// passing a `t` value outside of this range will result in an assert on debug builds
    /// and likely a crash on release builds.
    pub fn knot_domain(&self) -> (P::Scalar, P::Scalar) {
        (
            self.knots[D],
            self.knots[self.knots.len() - 1 - D],
        )
    }

    /// Iteratively compute de Boor's B-spline algorithm, this computes the recursive
    /// de Boor algorithm tree from the bottom up. At each level we use the results
    /// from the previous one to compute this level and store the results in the
    /// array indices we no longer need to compute the current level (the left one
    /// used computing node j).
    fn de_boor_iterative(&self, t: P::Scalar, start_knot: usize) -> P 
    where [(); D+1]: Sized 
    {
        let mut tmp: [P; D+1] = [P::default(); D+1];
        for j in 0..=D {
            let p = j + start_knot - D - 1;
            tmp[p] = self.control_points[p];
        }
        for lvl in 0..D {
            let k = lvl + 1;
            for j in 0..D - lvl {
                let i = j + k + start_knot - D;
                let alpha =
                    (t - self.knots[i - 1]) / (self.knots[i + D - k] - self.knots[i - 1]);
                debug_assert!(!alpha.is_nan());
                tmp[j] = tmp[j] * (-alpha + 1.0) + tmp[j + 1] * alpha;
            }
        }
        tmp[0]
    }

    /// Return the index of the first element greater than the value passed.
    /// Becaus the knot vector is sorted, this function uses binary search.
    /// If no element greater than the value passed is found, the function returns None.
    fn upper_bounds(&self, data: &[P::Scalar], value: P::Scalar) -> Option<usize> {
        let mut first = 0usize;
        let mut step;
        let mut count = data.len() as isize;
        while count > 0 {
            step = count / 2;
            let it = first + step as usize;
            if !value.lt(&data[it]) {
                first = it + 1;
                count -= step + 1;
            } else {
                count = step;
            }
        }
        // If we didn't find an element greater than value
        if first == data.len() {
            None
        } else {
            Some(first)
        }
    }

    /// Approximates the arc length of the curve by flattening it with straight line segments.
    /// This approximation is unfeasable if desired accuracy is greater than ~2 decimal places
    pub fn arclen(&self, nsteps: usize) -> P::Scalar 
    where [(); D+1]: Sized 
    {
        let stepsize = P::Scalar::from(1.0 / (nsteps as NativeFloat));
        let mut arclen: P::Scalar = 0.0.into();
        // evaluate the curve, t needs to be inside the knot domain!
        // we need to map [0...1] to kmin..kmax
        let (kmin, kmax) = self.knot_domain();
        let nsteps: usize = 100;
        for t in 0..=nsteps {
            let t =
                kmin + (P::Scalar::from(t as NativeFloat) / kmax) * 1.0 / (nsteps as NativeFloat);
            //dbg!(kmin, kmax, t);
            let p1 = self.eval(t);
            let p2 = self.eval(t + stepsize);
            arclen = (p1 - p2).squared_length().sqrt();
        }
        arclen
    }


    /// Returns the derivative curve of self which has N-1 control points.
    /// The derivative of an nth degree B-Spline curve is an (n-1)th degree (d) B-Spline curve,
    /// with the same knot vector, and new control points Q0...Qn-1 derived from the
    /// original control points Pi as: 
    ///                 d 
    /// Qi =    ----------------- (P(i+1)-P(i))
    ///         k[i+d+1] - k[i+1].
    /// with degree = curve_order - 1
    /// TODO test & verify function!
    pub fn derivative(&self) -> BSpline<P, K, {C-1}, {D-1}> {
        let mut new_points: [P; C-1] = [P::default(); C-1];
        for (i, _) in self.control_points.iter().enumerate() {
            new_points[i] =
                (self.control_points[i + 1] - self.control_points[i]) * ((D) as NativeFloat / (self.knots[i+D+1] - self.knots[i+1]).into());
            if i == self.control_points.len() - 2 {
                break;
            }
        }
        BSpline{knots: self.knots, control_points: new_points }
    }

}

#[cfg(test)]
mod tests {
    //use std;
    use super::PointN;
    use super::*;
    //use crate::num_traits::{Pow};
    #[test]
    fn construct_and_eval() {
        // degree 3, 4 control points => 4+3+1=8 knots
        let points = [
            PointN::new([0f64, 1.77f64]),
            PointN::new([1.1f64, -1f64]),
            PointN::new([4.3f64, 3f64]),
            PointN::new([3.2f64, -4f64]),
        ];
        let knots: [f64; 8] = [0., 0., 0., 1., 2., 3., 3., 3.];
        // try to make a b-spline with the given parameters
        let b: Option<BSpline<PointN<f64, 2>, 8, 4, 4>> = BSpline::new(knots, points);
        let curve = match b {
            None => return,
            Some(b) => b,
        };
        // evaluate the curve, t needs to be inside the knot domain!
        // we need to map [0...1] to kmin..kmax
        let (kmin, kmax) = curve.knot_domain();
        let nsteps: usize = 100;
        for t in 0..=nsteps {
            let t = kmin + (t as f64 / kmax) * 1f64 / (nsteps as f64);
            //dbg!(kmin, kmax, t);
            curve.eval(t);
        }
        // for now only check if has positive arclen
        assert!(curve.arclen(100) > 0.);
    }
}
