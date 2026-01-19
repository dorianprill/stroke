use core::iter::IntoIterator;
use core::slice;

//use crate::roots::RootFindingError;

use super::*;
use crate::find_root::FindRoot;
use crate::point::Point;
use crate::roots::RootFindingError;
use crate::spline::Spline;

const MAX_ROOT_DEPTH: usize = 32;
const ROOT_TOLERANCE: NativeFloat = 1e-6;

/// General implementation of a Bezier curve of arbitrary degree (= number of control points - 1).
///
/// The curve is solely defined by an array of 'control_points'. The degree is defined as degree = control_points.len() - 1.
/// Points on the curve can be evaluated with an interpolation parameter 't' in interval [0,1] using the eval() and eval_casteljau() methods.
/// Generic parameters:
/// P: Generic points 'P' as defined by there Point trait
/// const generic parameters:
/// N: Number of control points
#[derive(Clone, Copy)]
pub struct Bezier<P, const N: usize>
where
    P: Point,
{
    /// Control points which define the curve and hence its degree
    control_points: [P; N],
}

impl<P, const N: usize> Spline<P> for Bezier<P, { N }>
where
    P: Point,
{
    fn eval(&self, t: P::Scalar) -> P {
        self.eval(t)
    }
}

impl<P: Point, const N: usize> IntoIterator for Bezier<P, { N }> {
    type Item = P;
    type IntoIter = core::array::IntoIter<Self::Item, N>;

    fn into_iter(self) -> Self::IntoIter {
        IntoIterator::into_iter(self.control_points)
    }
}

impl<'a, P: Point, const N: usize> IntoIterator for &'a mut Bezier<P, { N }> {
    type Item = &'a mut P;
    type IntoIter = slice::IterMut<'a, P>;

    fn into_iter(self) -> slice::IterMut<'a, P> {
        self.control_points.iter_mut()
    }
}

impl<P, const N: usize> Bezier<P, { N }>
where
    P: Point,
{
    /// Create a new Bezier curve that interpolates the `control_points`. The degree is defined as degree = control_points.len() - 1.
    /// Desired curve must have a valid number of control points and knots in relation to its degree or the constructor will return None.
    /// A B-Spline curve requires at least one more control point than the degree (`control_points.len() >
    /// degree`) and the number of knots should be equal to `control_points.len() + degree + 1`.
    pub fn new(control_points: [P; N]) -> Bezier<P, { N }> {
        Bezier { control_points }
    }

    pub fn control_points(&self) -> [P; N] {
        self.control_points
    }

    /// Evaluate a point on the curve at point 't' which should be in the interval [0,1]
    /// This is implemented using De Casteljau's algorithm (over a temporary array with const generic sizing)
    pub fn eval(&self, t: P::Scalar) -> P {
        //let t = t.into();
        // start with a copy of the original control points array and succesively use it for evaluation
        let mut p: [P; N] = self.control_points;
        // loop up to degree = control_points.len() -1
        for i in 1..=p.len() {
            for j in 0..p.len() - i {
                p[j] = p[j] * (-t + 1.0) + p[j + 1] * t;
            }
        }
        p[0]
    }

    /// Calculates the minimum distance between given 'point' and the curve.
    /// Uses two passes with the same amount of steps in t:
    /// 1. coarse search over the whole curve
    /// 2. fine search around the minimum yielded by the coarse search
    pub fn distance_to_point(&self, point: P) -> P::Scalar {
        let nsteps: usize = 64;
        let mut tmin: P::Scalar = 0.5.into();
        let mut dmin: P::Scalar = (point - self.control_points[0]).squared_length();
        // 1. coarse pass
        for i in 0..nsteps {
            // calculate next step value
            let t: P::Scalar =
                (i as NativeFloat * 1.0 as NativeFloat / (nsteps as NativeFloat)).into();
            // calculate distance to candidate
            let candidate = self.eval(t);
            if (candidate - point).squared_length() < dmin {
                tmin = t;
                dmin = (candidate - point).squared_length();
            }
        }
        // 2. fine pass
        for i in 0..nsteps {
            // calculate next step value ( a 64th of a 64th from first step)
            let t: P::Scalar =
                (i as NativeFloat * 1.0 as NativeFloat / ((nsteps * nsteps) as NativeFloat)).into();
            // calculate distance to candidate centered around tmin from before
            let candidate: P = self.eval(tmin + t - t * (nsteps as NativeFloat / 2.0));
            if (candidate - point).squared_length() < dmin {
                tmin = t;
                dmin = (candidate - point).squared_length();
            }
        }
        dmin.sqrt()
    }

    pub fn split(&self, t: P::Scalar) -> (Self, Self) {
        // start with a copy of the original control points for now
        // TODO how to initialize const generic array without using unsafe?
        let mut left: [P; N] = self.control_points;
        let mut right: [P; N] = self.control_points;
        // these points get overriden each iteration; we save the intermediate results to 'left' and 'right'
        let mut casteljau_points: [P; N] = self.control_points;

        for i in 1..=casteljau_points.len() {
            // save start point of level
            left[i - 1] = casteljau_points[0];
            // save end point of level
            right[right.len() - i] = casteljau_points[right.len() - i];
            // calculate next level of points (one less point each level until we reach one point, the one at t)
            for j in 0..casteljau_points.len() - i {
                casteljau_points[j] =
                    casteljau_points[j] * (-t + 1.0) + casteljau_points[j + 1] * t;
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
    pub fn derivative(&self) -> Bezier<P, { N - 1 }> {
        let mut new_points: [P; N - 1] = [P::default(); N - 1];
        for (i, _) in self.control_points.iter().enumerate() {
            new_points[i] =
                (self.control_points[i + 1] - self.control_points[i]) * ((N - 1) as NativeFloat);
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
    //     eps: Option<P::Scalar>,
    //     max_iter: Option<usize>
    // ) -> Result<ArrayVec<[P::Scalar; N-1]>, RootFindingError>
    // {
    //     todo!();
    //     // Compute the axis-wise polynomial coefficients e.g. quadratic has N coefs a,b,c in at^2 + bt + c
    //     // to do this generically, we need to find the coefs of the bezier of degree n by binomial expansion
    //     // B_n(t) = sum_1_to_n ( binom(n,i) * s^(n-i) * t^i * p[i])
    //     let mut res: ArrayVec<[P::Scalar; N-1]> = ArrayVec::new();
    //     let mut npascal:    [P::Scalar; N] = [ P::Scalar::from(0.0); N];
    //     let poly_coefs:     [P::Scalar; N] = [ P::Scalar::from(0.0); N];

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

    //     let mut x = P::Scalar::from(0.0);

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
    pub fn arclen(&self, nsteps: usize) -> P::Scalar {
        let stepsize = P::Scalar::from(1.0 / (nsteps as NativeFloat));
        let mut arclen: P::Scalar = 0.0.into();
        for t in 1..nsteps {
            let t = P::Scalar::from(t as NativeFloat * 1.0 / (nsteps as NativeFloat));
            let p1 = self.eval(t);
            let p2 = self.eval(t + stepsize);

            arclen = arclen + (p1 - p2).squared_length().sqrt();
        }
        arclen
    }

    /// Find parameter values where the derivative crosses zero for a given axis.
    pub fn derivative_roots(&self, axis: usize) -> ArrayVec<[P::Scalar; N - 1]>
    where
        [(); N - 1]: Sized,
        [P::Scalar; N - 1]: tinyvec::Array<Item = P::Scalar>,
    {
        self.derivative_roots_with_tolerance(axis, P::Scalar::from(ROOT_TOLERANCE))
    }

    /// Find parameter values where the derivative crosses zero for a given axis using a tolerance.
    pub fn derivative_roots_with_tolerance(
        &self,
        axis: usize,
        tolerance: P::Scalar,
    ) -> ArrayVec<[P::Scalar; N - 1]>
    where
        [(); N - 1]: Sized,
        [P::Scalar; N - 1]: tinyvec::Array<Item = P::Scalar>,
    {
        if N < 2 {
            return ArrayVec::new();
        }
        let control = self.derivative_axis_control_points(axis);
        let mut roots = ArrayVec::<[P::Scalar; N - 1]>::new();
        let t0 = P::Scalar::from(0.0 as NativeFloat);
        let t1 = P::Scalar::from(1.0 as NativeFloat);
        Self::find_roots_1d(control, t0, t1, 0, tolerance, &mut roots);

        roots.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap());
        let mut unique = ArrayVec::<[P::Scalar; N - 1]>::new();
        for &root in roots.iter() {
            if unique
                .last()
                .map(|v| (root - *v).abs() > tolerance)
                .unwrap_or(true)
            {
                unique.push(root);
            }
        }
        unique
    }

    /// Return the bounding box of the curve as an array of (min, max) tuples for each dimension.
    pub fn bounding_box(&self) -> [(P::Scalar, P::Scalar); P::DIM]
    where
        [(); N - 1]: Sized,
        [P::Scalar; N - 1]: tinyvec::Array<Item = P::Scalar>,
    {
        let tolerance = P::Scalar::from(ROOT_TOLERANCE);
        let mut bounds = [(P::Scalar::default(), P::Scalar::default()); P::DIM];
        let zero = P::Scalar::from(0.0 as NativeFloat);
        let one = P::Scalar::from(1.0 as NativeFloat);

        for dim in 0..P::DIM {
            let mut min = self.control_points[0].axis(dim);
            let mut max = min;
            let end = self.control_points[N - 1].axis(dim);
            if end < min {
                min = end;
            }
            if end > max {
                max = end;
            }

            let roots = self.derivative_roots_with_tolerance(dim, tolerance);
            for &t in roots.iter() {
                if t > zero && t < one {
                    let value = self.eval(t).axis(dim);
                    if value < min {
                        min = value;
                    }
                    if value > max {
                        max = value;
                    }
                }
            }
            bounds[dim] = (min, max);
        }

        bounds
    }

    /// Find a root for a particular axis using Newton-Raphson on the scalar axis function.
    pub fn root_newton_axis(
        &self,
        value: P::Scalar,
        axis: usize,
        start: P::Scalar,
        eps: Option<P::Scalar>,
        max_iter: Option<usize>,
    ) -> Result<P::Scalar, RootFindingError>
    where
        [(); N - 1]: Sized,
    {
        FindRoot::root_newton_axis(self, value, axis, start, eps, max_iter)
    }

    fn derivative_axis_control_points(&self, axis: usize) -> [P::Scalar; N - 1]
    where
        [(); N - 1]: Sized,
    {
        let scale = (N - 1) as NativeFloat;
        core::array::from_fn(|i| {
            (self.control_points[i + 1].axis(axis) - self.control_points[i].axis(axis)) * scale
        })
    }

    fn find_roots_1d<const M: usize>(
        control: [P::Scalar; M],
        t0: P::Scalar,
        t1: P::Scalar,
        depth: usize,
        tolerance: P::Scalar,
        results: &mut ArrayVec<[P::Scalar; M]>,
    ) where
        [P::Scalar; M]: tinyvec::Array<Item = P::Scalar>,
    {
        let mut min = control[0];
        let mut max = control[0];
        for value in &control[1..] {
            if *value < min {
                min = *value;
            }
            if *value > max {
                max = *value;
            }
        }

        if min > P::Scalar::from(0.0 as NativeFloat)
            || max < P::Scalar::from(0.0 as NativeFloat)
        {
            return;
        }

        if depth >= MAX_ROOT_DEPTH || (max - min).abs() <= tolerance {
            let denom = control[M - 1] - control[0];
            let root = if denom.abs() > tolerance {
                t0 + (P::Scalar::from(0.0 as NativeFloat) - control[0]) * (t1 - t0) / denom
            } else {
                (t0 + t1) * P::Scalar::from(0.5 as NativeFloat)
            };
            if results.len() < results.capacity() {
                results.push(root);
            }
            return;
        }

        let (left, right) = Self::subdivide_1d(control);
        let mid = (t0 + t1) * P::Scalar::from(0.5 as NativeFloat);
        Self::find_roots_1d(left, t0, mid, depth + 1, tolerance, results);
        Self::find_roots_1d(right, mid, t1, depth + 1, tolerance, results);
    }

    fn subdivide_1d<const M: usize>(
        control: [P::Scalar; M],
    ) -> ([P::Scalar; M], [P::Scalar; M]) {
        let half = P::Scalar::from(0.5 as NativeFloat);
        let mut left = [P::Scalar::default(); M];
        let mut right = [P::Scalar::default(); M];
        let mut temp = control;

        for i in 0..M {
            left[i] = temp[0];
            right[M - 1 - i] = temp[M - 1 - i];
            for j in 0..(M - 1 - i) {
                temp[j] = (temp[j] + temp[j + 1]) * half;
            }
        }

        (left, right)
    }
}

impl<P, const N: usize> FindRoot<P> for Bezier<P, { N }>
where
    P: Point,
    [(); N - 1]: Sized,
{
    fn parameter_domain(&self) -> (P::Scalar, P::Scalar) {
        (
            P::Scalar::from(0.0 as NativeFloat),
            P::Scalar::from(1.0 as NativeFloat),
        )
    }

    fn axis_value(&self, t: P::Scalar, axis: usize) -> Result<P::Scalar, RootFindingError> {
        Ok(self.eval(t).axis(axis))
    }

    fn axis_derivative(&self, t: P::Scalar, axis: usize) -> Result<P::Scalar, RootFindingError> {
        Ok(self.derivative().eval(t).axis(axis))
    }
}

#[cfg(test)]
mod tests {
    use super::CubicBezier;
    use super::PointN;
    use super::QuadraticBezier;
    use super::*;

    //use crate::num_traits::{Pow};
    #[test]
    fn eval_endpoints() {
        let points = [
            PointN::new([0f64, 1.77f64]),
            PointN::new([1.1f64, -1f64]),
            PointN::new([4.3f64, 3f64]),
            PointN::new([3.2f64, -4f64]),
            PointN::new([7.3f64, 2.7f64]),
            PointN::new([8.9f64, 1.7f64]),
        ];

        let curve: Bezier<PointN<f64, 2>, 6> = Bezier::new(points);

        // check if start/end points match
        let start = curve.eval(0.0);
        let err_start = start - points[0];
        assert!(err_start.squared_length() < EPSILON);

        let end = curve.eval(1.0);
        let err_end = end - points[points.len() - 1];
        assert!(err_end.squared_length() < EPSILON);
    }

    #[test]
    fn distance_to_point() {
        // chose some arbitrary control points and construct a cubic bezier
        let bezier = Bezier {
            control_points: [
                PointN::new([0f64, 1.77f64]),
                PointN::new([2.9f64, 0f64]),
                PointN::new([4.3f64, 3f64]),
                PointN::new([3.2f64, -4f64]),
            ],
        };
        assert!(
            bezier.distance_to_point(PointN::new([-5.1, -5.6]))
                > bezier.distance_to_point(PointN::new([5.1, 5.6]))
        );
    }

    #[test]
    fn split_equivalence() {
        // chose some arbitrary control points and construct a cubic bezier
        let bezier = Bezier {
            control_points: [
                PointN::new([0f64, 1.77f64]),
                PointN::new([2.9f64, 0f64]),
                PointN::new([4.3f64, 3f64]),
                PointN::new([3.2f64, -4f64]),
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
            assert!(err.squared_length() < EPSILON);
            // check the right part of the split curve
            err = bezier.eval((t * 0.5) + 0.5) - right.eval(t);
            assert!(err.squared_length() < EPSILON);
        }
    }

    #[test]
    /// Check whether the generic implementation is
    /// equivalent to the specialized cubic implementation
    fn equivalence_cubic_specialization() {
        let cubic_bezier = CubicBezier::new(
            PointN::new([0f64, 1.77f64]),
            PointN::new([1.1f64, -1f64]),
            PointN::new([4.3f64, 3f64]),
            PointN::new([3.2f64, -4f64]),
        );

        let generic_bezier = Bezier {
            control_points: [
                PointN::new([0f64, 1.77f64]),
                PointN::new([1.1f64, -1f64]),
                PointN::new([4.3f64, 3f64]),
                PointN::new([3.2f64, -4f64]),
            ],
        };

        let nsteps: usize = 1000;
        for t in 0..=nsteps {
            let t = t as f64 * 1f64 / (nsteps as f64);
            let err = cubic_bezier.eval(t) - generic_bezier.eval(t);
            assert!(err.squared_length() < EPSILON);
        }
    }

    #[test]
    /// Check whether the generic implementation is
    /// equivalent to the specialized quadratic implementation
    fn equivalence_quadratic_specialization() {
        let quadratic_bezier = QuadraticBezier::new(
            PointN::new([0f64, 1.77f64]),
            PointN::new([1.1f64, -1f64]),
            PointN::new([3.2f64, -4f64]),
        );

        let generic_bezier = Bezier {
            control_points: [
                PointN::new([0f64, 1.77f64]),
                PointN::new([1.1f64, -1f64]),
                PointN::new([3.2f64, -4f64]),
            ],
        };

        let nsteps: usize = 1000;
        for t in 0..=nsteps {
            let t = t as f64 * 1f64 / (nsteps as f64);
            let err = quadratic_bezier.eval(t) - generic_bezier.eval(t);
            assert!(err.squared_length() < EPSILON);
        }
    }

    #[test]
    fn derivative_root_quadratic_1d() {
        let bezier = Bezier {
            control_points: [
                PointN::new([0f64]),
                PointN::new([1f64]),
                PointN::new([0f64]),
            ],
        };

        let roots = bezier.derivative_roots(0);
        assert_eq!(roots.len(), 1);
        assert!((roots[0] - 0.5).abs() < 1e-6);
    }

    #[test]
    fn derivative_roots_cubic_two_roots() {
        let bezier = Bezier {
            control_points: [
                PointN::new([0f64]),
                PointN::new([4f64]),
                PointN::new([-4f64]),
                PointN::new([0f64]),
            ],
        };

        let roots = bezier.derivative_roots(0);
        assert_eq!(roots.len(), 2);

        let expected = 0.5;
        let offset = (3.0f64).sqrt() / 6.0;
        let r0 = expected - offset;
        let r1 = expected + offset;
        assert!((roots[0] - r0).abs() < 1e-4);
        assert!((roots[1] - r1).abs() < 1e-4);
    }

    #[test]
    fn derivative_roots_no_roots() {
        let bezier = Bezier {
            control_points: [
                PointN::new([0f64]),
                PointN::new([1f64]),
                PointN::new([2f64]),
            ],
        };

        let roots = bezier.derivative_roots(0);
        assert!(roots.is_empty());
    }

    #[test]
    fn bounding_box_cubic_1d_extrema() {
        let bezier = Bezier {
            control_points: [
                PointN::new([0f64]),
                PointN::new([4f64]),
                PointN::new([-4f64]),
                PointN::new([0f64]),
            ],
        };

        let offset = (3.0f64).sqrt() / 6.0;
        let t0 = 0.5 - offset;
        let t1 = 0.5 + offset;
        let v0 = bezier.eval(t0).axis(0);
        let v1 = bezier.eval(t1).axis(0);
        let expected_min = v0.min(v1);
        let expected_max = v0.max(v1);

        let bounds = bezier.bounding_box();
        assert!((bounds[0].0 - expected_min).abs() < 1e-4);
        assert!((bounds[0].1 - expected_max).abs() < 1e-4);
    }

    #[test]
    fn root_newton_axis_linear() {
        let bezier = Bezier {
            control_points: [PointN::new([0f64]), PointN::new([2f64])],
        };

        let root = bezier
            .root_newton_axis(1.0, 0, 0.25, None, Some(64))
            .unwrap();
        assert!((root - 0.5).abs() < 1e-6);
    }

    #[test]
    fn root_newton_axis_zero_derivative() {
        let bezier = Bezier {
            control_points: [PointN::new([1f64]), PointN::new([1f64])],
        };

        let result = bezier.root_newton_axis(0.0, 0, 0.5, None, Some(8));
        assert!(matches!(result, Err(RootFindingError::ZeroDerivative)));
    }

    #[test]
    fn bounding_box_matches_cubic_specialization() {
        let generic_bezier = Bezier {
            control_points: [
                PointN::new([0f64, 1.77f64]),
                PointN::new([1.1f64, -1f64]),
                PointN::new([4.3f64, 3f64]),
                PointN::new([3.2f64, -4f64]),
            ],
        };
        let cubic_bezier = CubicBezier::new(
            PointN::new([0f64, 1.77f64]),
            PointN::new([1.1f64, -1f64]),
            PointN::new([4.3f64, 3f64]),
            PointN::new([3.2f64, -4f64]),
        );

        let generic_bounds = generic_bezier.bounding_box();
        let cubic_bounds = cubic_bezier.bounding_box();
        let tol = 1e-5;
        for i in 0..2 {
            assert!((generic_bounds[i].0 - cubic_bounds[i].0).abs() < tol);
            assert!((generic_bounds[i].1 - cubic_bounds[i].1).abs() < tol);
        }
    }
}
