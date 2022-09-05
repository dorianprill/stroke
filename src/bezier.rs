use core::fmt::Debug;
use core::iter::IntoIterator;
use core::slice;
use nalgebra::{
    ClosedMul, ComplexField,
    RealField, SVector,
};
use num_traits::{NumCast, Zero};

//use crate::roots::RootFindingError;

use super::*;
use crate::spline::Spline;

//pub trait Scalar : nalgebra::Scalar + Clone + PartialOrd + SimdValue<Element = Self> + SimdPartialOrd + ClosedAdd + ClosedDiv + ClosedMul + ClosedSub + Clone + Field {}

pub trait Scalar: RealField + Copy + NumCast + Default + ClosedMul + Debug + Float {}

impl<T: RealField + Copy + NumCast + Default + ClosedMul + Debug + Float> Scalar for T {}

/// General implementation of a Bezier curve of arbitrary degree (= number of control points - 1).
/// The curve is solely defined by an array of 'control_points'. The degree is defined as degree = control_points.len() - 1.
/// Points on the curve can be evaluated with an interpolation parameter 't' in interval [0,1] using the eval() and eval_casteljau() methods.
/// Generic parameters:
/// T: Generic points 'T' as defined by there Point trait
/// const generic parameters:
/// N: Number of control points
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct Bezier<T: Scalar, const DIM: usize, const N: usize> {
    /// Control points which define the curve and hence its degree
    pub control_points: [SVector<T, DIM>; N],
}

impl<T: Scalar, const DIM: usize, const N: usize> Spline<T, DIM> for Bezier<T, DIM, N> {
    fn eval(&self, t: T) -> SVector<T, DIM> {
        self.eval(t)
    }
}

impl<T: Scalar, const DIM: usize, const N: usize> IntoIterator for Bezier<T, DIM, N> {
    type Item = SVector<T, DIM>;
    type IntoIter = core::array::IntoIter<Self::Item, N>;

    fn into_iter(self) -> Self::IntoIter {
        IntoIterator::into_iter(self.control_points)
    }
}
impl<'a, T: Scalar, const DIM: usize, const N: usize> IntoIterator for &'a mut Bezier<T, DIM, N> {
    type Item = &'a mut SVector<T, DIM>;
    type IntoIter = slice::IterMut<'a, SVector<T, DIM>>;

    fn into_iter(self) -> slice::IterMut<'a, SVector<T, DIM>> {
        self.control_points.iter_mut()
    }
}

impl<T: Scalar, const DIM: usize, const N: usize> Bezier<T, DIM, N> {
    /// Create a new Bezier curve that interpolates the `control_points`. The degree is defined as degree = control_points.len() - 1.
    /// Desired curve must have a valid number of control points and knots in relation to its degree or the constructor will return None.
    /// A B-Spline curve requires at least one more control point than the degree (`control_points.len() >
    /// degree`) and the number of knots should be equal to `control_points.len() + degree + 1`.
    pub fn new(control_points: [SVector<T, DIM>; N]) -> Self {
        Bezier { control_points }
    }

    pub fn control_points(&self) -> [SVector<T, DIM>; N] {
        self.control_points.clone()
    }

    /// Evaluate a point on the curve at point 't' which should be in the interval [0,1]
    /// This is implemented using De Casteljau's algorithm (over a temporary array with const generic sizing)
    pub fn eval(&self, t: T) -> SVector<T, DIM> {
        //let t = t.into();
        // start with a copy of the original control points array and succesively use it for evaluation
        let mut p: [SVector<T, DIM>; N] = self.control_points;
        p[0].add_scalar(t);
        // loop up to degree = control_points.len() -1
        for i in 1..=p.len() {
            for j in 0..p.len() - i {
                p[j] = p[j] * (-t + T::from(1.0).unwrap()) + p[j + 1] * t;
            }
        }
        p[0]
    }

    pub fn split(&self, t: T) -> (Self, Self) {
        // start with a copy of the original control points for now
        // TODO how to initialize const generic array without using unsafe?
        let mut left: [SVector<T, DIM>; N] = self.control_points;
        let mut right: [SVector<T, DIM>; N] = self.control_points;
        // these points get overriden each iteration; we save the intermediate results to 'left' and 'right'
        let mut casteljau_points: [SVector<T, DIM>; N] = self.control_points;

        for i in 1..=casteljau_points.len() {
            // save start point of level
            left[i - 1] = casteljau_points[0];
            // save end point of level
            right[right.len() - i] = casteljau_points[right.len() - i];
            // calculate next level of points (one less point each level until we reach one point, the one at t)
            for j in 0..casteljau_points.len() - i {
                casteljau_points[j] = casteljau_points[j] * (-t + T::from(1.0).unwrap())
                    + casteljau_points[j + 1] * t;
            }
        }

        (
            Bezier {
                control_points: left,
            },
            Bezier {
                control_points: right,
            },
        )
    }

    /// Returns the derivative curve of self which has N-1 control points.
    /// The derivative of an nth degree Bézier curve is an (n-1)th degree Bézier curve,
    /// with one fewer term, and new weights w0...wn-1 derived from the
    /// original weights as n(wi+1 - wi). So for a 3rd degree curve, with four weights,
    /// the derivative has three new weights:
    ///     w0 = 3(w1-w0), w'1 = 3(w2-w1) and w'2 = 3(w3-w2).
    pub fn derivative(&self) -> Bezier<T, DIM, { N - 1 }> {
        let _x = self.control_points[0];
        let mut new_points = [[T::zero(); DIM].into(); N - 1];
        for (i, _) in self.control_points.iter().enumerate() {
            new_points[i] =
                (self.control_points[i + 1] - self.control_points[i]) * (T::from(N - 1).unwrap());
            if i == self.control_points.len() - 2 {
                break;
            }
        }
        Bezier::new(new_points)
    }

    // /// Returns the real roots of the Bezier curve along one of its coordinate
    // /// axes (i.e. the control points' axes) or a specific RootFindingError.
    // /// There are the same number of roots as the degree of the curve nroots = degree = N_points-1
    // fn real_roots(&self,
    //     axis: usize,
    //     eps: Option<T>,
    //     max_iter: Option<usize>
    // ) -> Result<ArrayVec<[T; N-1]>, RootFindingError>
    // {
    //     todo!();
    //     // Compute the axis-wise polynomial coefficients e.g. quadratic has N coefs a,b,c in at^2 + bt + c
    //     // to do this generically, we need to find the coefs of the bezier of degree n by binomial expansion
    //     // B_n(t) = sum_1_to_n ( binom(n,i) * s^(n-i) * t^i * p[i])
    //     let mut res: ArrayVec<[T; N-1]> = ArrayVec::new();
    //     let mut npascal:    [T; N] = [ T::from(0.0); N];
    //     let poly_coefs:     [T; N] = [ T::from(0.0); N];

    //     // 1. calculate the n-th row of pascals triangle on a zero-based index (all values for i in the binom(n,i) part)
    //     //    1      N = 0 (wouldn't compile due to index out of bounds)
    //     //  1  (1)   N = 1 (last 1 is always omitted)
    //     // 1  2  (1) N = 2
    //     npascal[0] = 1.0.into();
    //     for i in 1usize..N {
    //         npascal[i] = npascal[i-1] * (N - i + 1) as NativeFloat / i as NativeFloat;
    //     }
    //     // 2. calculate the coefficients to binom(n,i) and p[i] (the s^(n-i) part)

    //     // 3. find candidate points for roots of that curve (compare zero crossings)

    //     // 4. search for roots using newton-raphson algo
    //     let eps = eps.unwrap_or(1e-3.into());
    //     let max_iter = max_iter.unwrap_or(128);

    //     let mut x = T::from(0.0);

    //     let mut iter = 0;
    //     loop {
    //         let f = f(x);
    //         let d = d(x);
    //         if f < eps {
    //             return Ok(x);
    //         }
    //         // if derivative is 0
    //         if d < EPSILON {
    //             // either try to choose a better starting point
    //             if iter == 0 {
    //                 x = x + 1.0;
    //                 iter = iter + 1;
    //                 continue;
    //             // or fail
    //             } else {
    //                 return Err(RootFindingError::ZeroDerivative);
    //             }
    //         }

    //         let x1 = x - f / d;
    //         if (x - x1).abs() < eps {
    //             return Ok(x1);
    //         }

    //         x = x1;
    //         iter = iter + 1;

    //         if iter == max_iter {
    //             return Err(RootFindingError::FailedToConverge);
    //         }
    //     }
    //     res
    // }

    /// Approximates the arc length of the curve by flattening it with straight line segments.
    /// This works quite well, at ~32 segments it should already provide an error in the decimal places
    /// The accuracy gain falls off with more steps so this approximation is unfeasable if desired accuracy is greater than 1-2 decimal places
    pub fn arclen(&self, nsteps: usize) -> T {
        let stepsize = T::from(1.0).unwrap() / (T::from(nsteps).unwrap());
        let mut arclen: T = T::zero();
        for t in 1..nsteps {
            let t = T::from(t).unwrap() / (T::from(nsteps).unwrap());
            let p1 = self.eval(t);
            let p2 = self.eval(t + stepsize);

            arclen = arclen + ComplexField::sqrt((p1 - p2).magnitude_squared());
        }
        return arclen;
    }

    pub fn start(&self) -> SVector<T, DIM> {
        self.control_points[0]
    }

    pub fn end(&self) -> SVector<T, DIM> {
        self.control_points[self.control_points.len() - 1]
    }

    pub fn baseline(&self) -> Bezier<T, DIM, 2> {
        Bezier::new([self.start(), self.end()])
    }

    pub fn is_linear(&self, tolerance: T) -> bool {
        if ({ DIM } == 2) {
            true
        } else {
            let line = self.baseline();
            self.control_points[1..self.control_points.len() - 1]
                .into_iter()
                .any(|point| line.distance_to_point(point) > tolerance)
        }
    }

    /// Return the bounding box of the curve as an array of (min, max) tuples for each dimension (its index)
    pub fn bounding_box(&self) -> (SVector<T, DIM>, SVector<T, DIM>) {
        ([T::zero(); DIM].into(), [T::zero(); DIM].into())
    }
}

impl<T: Scalar, const DIM: usize> Bezier<T, DIM, 2> {
    pub fn distance_to_point(&self, _point: &SVector<T, DIM>) -> T {
        T::zero()
    }
}

#[cfg(test)]
mod tests {
    use nalgebra::Vector2;

    
    
    use super::*;

    //use crate::num_traits::{Pow};
    #[test]
    fn eval_endpoints() {
        let points = [
            Vector2::new(0f64, 1.77f64),
            [1.1f64, -1f64].into(),
            [4.3f64, 3f64].into(),
            [3.2f64, -4f64].into(),
            [7.3f64, 2.7f64].into(),
            [8.9f64, 1.7f64].into(),
        ];

        let curve = Bezier::new(points);
        let _b = curve;

        // check if start/end points match
        let start = curve.eval(0.0);
        let err_start = start - points[0];
        assert!(err_start.magnitude_squared() < EPSILON);

        let end = curve.eval(1.0);
        let err_end = end - points[points.len() - 1];
        assert!(err_end.magnitude_squared() < EPSILON);
    }

    #[test]
    fn split_equivalence() {
        // chose some arbitrary control points and construct a cubic bezier
        let bezier = Bezier {
            control_points: [
                [0f64, 1.77f64].into(),
                [2.9f64, 0f64].into(),
                [4.3f64, 3f64].into(),
                [3.2f64, -4f64].into(),
            ],
        };
        // split it at an arbitrary point
        let at = 0.5;
        let (left, right) = bezier.split(at);
        // compare left and right subcurves with parent curve
        // take the difference of the two points which must not exceed the absolute error
        let nsteps: usize = 1000;
        for t in 0..=nsteps {
            let t = t as f64 * 1f64 / (nsteps as f64);
            // check the left part of the split curve
            let mut err = bezier.eval(t / 2.0) - left.eval(t);
            assert!(err.magnitude_squared() < EPSILON);
            // check the right part of the split curve
            err = bezier.eval((t * 0.5) + 0.5) - right.eval(t);
            assert!(err.magnitude_squared() < EPSILON);
        }
    }
}
