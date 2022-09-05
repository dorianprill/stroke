use core::slice::*;

use crate::bezier::Scalar;

use super::*;
use nalgebra::{ComplexField, SVector};

pub struct If<const B: bool>;
pub trait True {}
impl True for If<true> {}

/// General Implementation of a BSpline with choosable degree, control points and knots.
/// Generic parameters:
/// SVector<T, DIM>: const generic points array 'SVector<T, DIM>' as defined by the Point trait
/// F: Any float value used for the knots and interpolation (usually the same as the internal generic parameter within SVector<T, DIM><F>).
/// const generic parameters:D
/// C: Number of control points
/// K: Number of Knots
/// D: Degree of the piecewise function used for interpolation degree = order - 1
/// While C, K, O relate to each other in the following manner
///     K = C + O where O = D + 1
/// K = C + D + 1
/// D = K - C - 1
/// C = K - D - 1
#[derive(Clone)]
pub struct BSpline<T: Scalar, const DIM: usize, const K: usize, const C: usize>
where
    [(); K - C]: Sized,
{
    /// Knot vector
    knots: [T; K],

    /// Control points
    control_points: [SVector<T, DIM>; C],
}

impl<T: Scalar, const DIM: usize, const K: usize, const C: usize> BSpline<T, DIM, { K }, { C }>
where
    [(); K - C]: Sized,
    // where
    //     If<{ K > D - 1 }>: True,
    //     [(); K - D - 1]: Sized,
    //     [(); D +1]: Sized,
{
    const fn D() -> usize {
        K - C - 1
    }

    /// Create a new B-spline curve that interpolates
    /// the `control_points` using a piecewise polynomial of `degree` within intervals specified by the `knots`.
    /// The knots _must_ be sorted in non-decreasing order, the constructor enforces this which may yield undesired results.
    /// The degree is defined as `curve_order - 1`.
    /// Desired curve must have a valid number of control points and knots in relation to its degree or the constructor will return None.
    /// A B-Spline curve requires at least one more control point than the degree (`control_points.len() >
    /// degree`) and the number of knots should be equal to `control_points.len() + degree + 1`.
    #[allow(clippy::if_same_then_else)] // allow until a proper Error type is defined for bsplines
    pub fn new(
        knots: [T; K],
        control_points: [SVector<T, DIM>; C],
        //degree: usize,
    ) -> Option<BSpline<T, DIM, K, C>> {
        if control_points.len() <= Self::D() {
            //panic!("Too few control points for curve");
            None
        } else if knots.len() != control_points.len() + Self::D() + 1 {
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
    pub fn eval(&self, t: T) -> SVector<T, DIM> {
        debug_assert!(t >= self.knot_domain().0 && t <= self.knot_domain().1);

        let mut span = self.knots.partition_point(|&knot| knot <= t);
        if span > 0 { 
            span -= 1;
        }

        // Find the knot span that contains t i.e. the first index with a knot value greater than the t we're searching for.
        // We need to find the start of the knot span t is in, such that: knots[span] <= t < knots[span + 1]
        // Note: A custom function is used to exploit binary search (knots are sorted)
        // let span = match self.upper_bounds(&self.knots[..], t) {
        //     Some(x) if x == 0 => Self::D(), // degree
        //     Some(x) if x >= self.knots.len() - Self::D() - 1 => self.knots.len() - Self::D() - 1,
        //     Some(x) => x,
        //     None => self.knots.len() - Self::D() - 1,
        // };
        self.de_boor_iterative(t, span)
    }

    /// Returns an iterator over the control points.
    pub fn control_points(&self) -> Iter<'_, SVector<T, DIM>> {
        self.control_points.iter()
    }

    /// Returns an iterator over the knots.
    pub fn knots(&self) -> Iter<'_, T> {
        self.knots.iter()
    }

    /// Get the min and max knot domain values for finding the `t` range to compute
    /// the curve over. The curve is only defined over the inclusive range `[min, max]`,
    /// passing a `t` value outside of this range will result in an assert on debug builds
    /// and likely a crash on release builds.
    pub fn knot_domain(&self) -> (T, T) {
        (
            self.knots[Self::D()],
            self.knots[self.knots.len() - 1 - Self::D()],
        )
    }

    /// Iteratively compute de Boor's B-spline algorithm, this computes the recursive
    /// de Boor algorithm tree from the bottom up. At each level we use the results
    /// from the previous one to compute this level and store the results in the
    /// array indices we no longer need to compute the current level (the left one
    /// used computing node j).
    fn de_boor_iterative(&self, t: T, start_knot: usize) -> SVector<T, DIM> {
        let mut tmp = [[T::zero(); DIM].into(); K - C];
        #[cfg(test)]
        {
            dbg!(t, start_knot, tmp.len(), Self::D(), K - C, self.control_points.len());
        }

        for j in 0..Self::D() {
            let i = j + start_knot - Self::D();
            #[cfg(test)]
            {
                dbg!(i);
            }
            let cp = self.control_points[i];
            tmp[j] = cp;
        }

        for r in 1..=Self::D() {
            for j in (r..=Self::D()).rev() {
                let alpha = (t - self.knots[j + start_knot - Self::D()])
                    / (self.knots[j + 1 + start_knot - r] - self.knots[j + start_knot - Self::D()]);

                debug_assert!(!alpha.is_nan());
                tmp[j] = (tmp[j - 1].scale(T::from(1.0).unwrap() - alpha)) + (tmp[j].scale(alpha));
            }
        }
        #[cfg(test)]
        {
            dbg!(tmp);
        }
        tmp[Self::D()]
    }

    /// Return the index of the first element greater than the value passed.
    /// Becaus the knot vector is sorted, this function uses binary search.
    /// If no element greater than the value passed is found, the function returns None.
    fn upper_bounds(&self, data: &[T], value: T) -> Option<usize> {
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
    pub fn arclen(&self, nsteps: usize) -> T {
        let stepsize = T::from(1.0).unwrap() / T::from(nsteps).unwrap();
        let mut arclen = T::zero();
        // evaluate the curve, t needs to be inside the knot domain!
        // we need to map [0...1] to kmin..kmax
        let (kmin, kmax) = self.knot_domain();
        let nsteps: usize = 100;
        for t in 0..=nsteps {
            let t = kmin
                + (T::from(t as NativeFloat).unwrap() / kmax) * T::from(1.0).unwrap()
                    / T::from(nsteps).unwrap();
            //dbg!(kmin, kmax, t);
            let p1 = self.eval(t);
            let p2 = self.eval(t + stepsize);
            arclen = ComplexField::sqrt((p1 - p2).magnitude_squared());
        }
        arclen
    }

    /// Returns the derivative curve of self which has N-1 control points.
    /// The derivative of an nth degree B-Spline curve is an (n-1)th degree (d) B-Spline curve,
    /// with the same knot vector, and new control points Q0...Qn-1 derived from the
    /// original control points Pi as:
    ///                 d
    /// Qi =    ----------------- (SVector<T, DIM>(i+1)-SVector<T, DIM>(i))
    ///         k[i+d+1] - k[i+1].
    /// with degree = curve_order - 1
    /// TODO test & verify function!
    pub fn derivative(&self) -> BSpline<T, DIM, K, { C - 1 }>
    where
        [(); K - (C - 1)]: Sized,
    {
        let mut new_points: [SVector<T, DIM>; C - 1] = [[T::zero(); DIM].into(); C - 1];
        for (i, _) in self.control_points.iter().enumerate() {
            new_points[i] = (self.control_points[i + 1] - self.control_points[i])
                * (T::from(Self::D()).unwrap()
                    / (self.knots[i + Self::D() + 1] - self.knots[i + 1]).into());
            if i == self.control_points.len() - 2 {
                break;
            }
        }
        BSpline {
            knots: self.knots,
            control_points: new_points,
        }
    }
}
/*
impl<T: Scalar, const DIM: usize, const K: usize, const D: usize> BSpline<T, DIM, { K }, { D }>
where
    If<{ K > D - 1 }>: True,
    [(); K - D - 1]: Sized,
    [(); { K - D - 2 }]: Sized,
    If<{ K > { D - 1 } - 1 }>: True,
{

}
*/
#[cfg(test)]
mod tests {
    use std::dbg;

    use nalgebra::Vector2;

    //use std;
    use super::*;
    //use crate::num_traits::{Pow};

    #[test]
    fn construct_and_eval() {
        // degree 3, 4 control points => 4+3+1=8 knots
        let points = [
            [0f64, 1.77f64].into(),
            [1.1f64, -1f64].into(),
            [4.3f64, 3f64].into(),
            [3.2f64, -4f64].into(),
        ];
        let knots: [f64; 8] = [0., 0., 0., 1., 2., 3., 3., 3.];
        // try to make a b-spline with the given parameters
        let b = BSpline::new(knots, points);
        let curve = match b {
            None => return,
            Some(b) => b,
        };

        // evaluate the curve, t needs to be inside the knot domain!
        // we need to map [0...1] to kmin..kmax
        let (kmin, kmax) = curve.knot_domain();
        dbg!((kmin, kmax));
        dbg!(curve.eval(kmin));
        dbg!(curve.eval(kmax));
        dbg!(curve.eval(1.5));
        assert!((curve.eval(kmin) - Vector2::new(1.358333, 0.35916666)).magnitude_squared() < EPSILON);
        assert!((curve.eval(1.5) - Vector2::new(2.63125, 0.8678125)).magnitude_squared() < EPSILON);
        assert!((curve.eval(kmax) - Vector2::new(3.49166666666, 0.583333333)).magnitude_squared() < EPSILON);

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
