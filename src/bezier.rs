//! Const-generic Bezier curves of arbitrary degree.

use core::iter::IntoIterator;
use core::slice;

use num_traits::{Float, NumCast};
use super::*;
use crate::find_root::FindRoot;
use crate::point::{Point, PointIndex, PointNorm};
use crate::roots::RootFindingError;

const MAX_ROOT_DEPTH: usize = 32;
const ROOT_TOLERANCE: f64 = 1e-6;
const DEFAULT_DISTANCE_STEPS: usize = 64;
const DEFAULT_LENGTH_STEPS: usize = 64;
const MAX_LENGTH_ITERS: usize = 32;
const LOCAL_SEARCH_ITERS: usize = 16;

/// General implementation of a Bezier curve of arbitrary degree (`N - 1`).
///
/// The curve is defined by its control points. Points on the curve can be
/// evaluated with an interpolation parameter `t` in `[0, 1]` using `eval()`.
/// Methods that need component access or norms add `PointIndex`/`PointNorm`
/// bounds as required.
///
/// # Parameters
/// - `P`: point type implementing [`Point`](crate::point::Point).
/// - `N`: number of control points.
#[derive(Clone, Copy)]
pub struct Bezier<P, const N: usize>
where
    P: Point,
{
    /// Control points which define the curve and hence its degree
    control_points: [P; N],
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
    /// Create a new Bezier curve from the provided control points.
    /// The degree is `N - 1`.
    pub fn new(control_points: [P; N]) -> Bezier<P, { N }> {
        Bezier { control_points }
    }

    /// Return the control points array.
    pub fn control_points(&self) -> [P; N] {
        self.control_points
    }

    /// Return the start point of the curve.
    pub fn start(&self) -> P {
        let zero = <P::Scalar as NumCast>::from(0.0).unwrap();
        self.eval(zero)
    }

    /// Return the end point of the curve.
    pub fn end(&self) -> P {
        let one = <P::Scalar as NumCast>::from(1.0).unwrap();
        self.eval(one)
    }

    /// Return a curve with reversed direction.
    pub fn reverse(&self) -> Self {
        let control_points = core::array::from_fn(|i| self.control_points[N - 1 - i]);
        Bezier { control_points }
    }

    /// Evaluate a point on the curve at point 't' which should be in the interval [0,1] (unchecked!)
    /// This is implemented using De Casteljau's algorithm (over a temporary array with const generic sizing)
    pub fn eval(&self, t: P::Scalar) -> P {
        //let t = t.into();
        let one = <P::Scalar as NumCast>::from(1.0).unwrap();
        // start with a copy of the original control points array and succesively use it for evaluation
        let mut p: [P; N] = self.control_points;
        // loop up to degree = control_points.len() -1
        for i in 1..p.len() {
            for j in 0..p.len() - i {
                p[j] = p[j] * (one - t) + p[j + 1] * t;
            }
        }
        p[0]
    }

    /// Approximate the minimum distance between given `point` and the curve.
    /// Uses a coarse sampling pass over the full domain and a local search around
    /// the best sample. `nsteps` is the number of coarse samples.
    pub fn distance_to_point_approx(&self, point: P, nsteps: usize) -> P::Scalar
    where
        P: PointNorm,
    {
        let nsteps = nsteps.max(1);
        let nsteps_scalar = <P::Scalar as NumCast>::from(nsteps as f64).unwrap();
        let zero = <P::Scalar as NumCast>::from(0.0).unwrap();
        let half = <P::Scalar as NumCast>::from(0.5).unwrap();
        let three = <P::Scalar as NumCast>::from(3.0).unwrap();

        let mut best_i = 0usize;
        let mut best_d = (self.eval(zero) - point).squared_norm();
        for i in 1..=nsteps {
            let t = <P::Scalar as NumCast>::from(i as f64).unwrap() / nsteps_scalar;
            let candidate = self.eval(t);
            let d = (candidate - point).squared_norm();
            if d < best_d {
                best_d = d;
                best_i = i;
            }
        }

        let left_i = if best_i == 0 { 0 } else { best_i - 1 };
        let right_i = if best_i == nsteps { nsteps } else { best_i + 1 };
        let mut left = <P::Scalar as NumCast>::from(left_i as f64).unwrap() / nsteps_scalar;
        let mut right = <P::Scalar as NumCast>::from(right_i as f64).unwrap() / nsteps_scalar;

        for _ in 0..LOCAL_SEARCH_ITERS {
            let third = (right - left) / three;
            let t1 = left + third;
            let t2 = right - third;
            let d1 = (self.eval(t1) - point).squared_norm();
            let d2 = (self.eval(t2) - point).squared_norm();
            if d1 < d2 {
                right = t2;
            } else {
                left = t1;
            }
        }

        let t = (left + right) * half;
        (self.eval(t) - point).squared_norm().sqrt()
    }

    /// Approximate the minimum distance between given `point` and the curve using
    /// a default sampling resolution.
    pub fn distance_to_point(&self, point: P) -> P::Scalar
    where
        P: PointNorm,
    {
        self.distance_to_point_approx(point, DEFAULT_DISTANCE_STEPS)
    }

    /// Approximate parameter `t` at arc length `s`.
    pub fn t_at_length_approx(&self, s: P::Scalar, nsteps: usize) -> P::Scalar
    where
        P: PointNorm,
    {
        let zero = <P::Scalar as NumCast>::from(0.0).unwrap();
        let one = <P::Scalar as NumCast>::from(1.0).unwrap();
        let half = <P::Scalar as NumCast>::from(0.5).unwrap();

        let total = self.arclen(nsteps);
        if total <= P::Scalar::epsilon() {
            return zero;
        }
        let target = s.clamp(zero, total);

        let mut lo = zero;
        let mut hi = one;
        for _ in 0..MAX_LENGTH_ITERS {
            let mid = (lo + hi) * half;
            let len = self.arclen_partial(mid, nsteps);
            if len < target {
                lo = mid;
            } else {
                hi = mid;
            }
        }

        (lo + hi) * half
    }

    /// Approximate parameter `t` at arc length `s` using a default resolution.
    pub fn t_at_length(&self, s: P::Scalar) -> P::Scalar
    where
        P: PointNorm,
    {
        self.t_at_length_approx(s, DEFAULT_LENGTH_STEPS)
    }

    /// Evaluate the point at arc length `s`.
    pub fn point_at_length_approx(&self, s: P::Scalar, nsteps: usize) -> P
    where
        P: PointNorm,
    {
        self.eval(self.t_at_length_approx(s, nsteps))
    }

    /// Evaluate the point at arc length `s` using a default resolution.
    pub fn point_at_length(&self, s: P::Scalar) -> P
    where
        P: PointNorm,
    {
        self.point_at_length_approx(s, DEFAULT_LENGTH_STEPS)
    }

    /// Split the curve at `t` into two sub-curves.
    pub fn split(&self, t: P::Scalar) -> (Self, Self) {
        let one = <P::Scalar as NumCast>::from(1.0).unwrap();
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
                    casteljau_points[j] * (one - t) + casteljau_points[j + 1] * t;
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
    /// The derivative of an nth degree Bezier curve is an (n-1)th degree Bezier curve,
    /// with one fewer term, and new weights w0...wn-1 derived from the
    /// original weights as n(wi+1 - wi). So for a 3rd degree curve, with four weights,
    /// the derivative has three new weights:
    ///     w0 = 3(w1-w0), w'1 = 3(w2-w1) and w'2 = 3(w3-w2).
    pub fn derivative(&self) -> Bezier<P, { N - 1 }> {
        let scale = <P::Scalar as NumCast>::from((N - 1) as f64).unwrap();
        let new_points: [P; N - 1] = core::array::from_fn(|i| {
            (self.control_points[i + 1] - self.control_points[i]) * scale
        });
        Bezier::new(new_points)
    }

    /// Return the unit tangent direction at `t`.
    pub fn tangent(&self, t: P::Scalar) -> P
    where
        P: PointNorm,
        [(); N - 1]: Sized,
    {
        let one = <P::Scalar as NumCast>::from(1.0).unwrap();
        let dir = self.derivative().eval(t);
        let len = dir.squared_norm().sqrt();
        if len <= P::Scalar::epsilon() {
            dir
        } else {
            dir * (one / len)
        }
    }

    /// Return the curvature magnitude at `t`.
    pub fn curvature(&self, t: P::Scalar) -> P::Scalar
    where
        P: PointNorm + PointDot,
        [(); N - 1]: Sized,
        [(); N - 2]: Sized,
        [(); (N - 1) - 1]: Sized,
    {
        let zero = <P::Scalar as NumCast>::from(0.0).unwrap();
        let first = self.derivative();
        let v = first.eval(t);
        let a = if N < 3 {
            v - v
        } else {
            first.derivative().eval(t)
        };
        let v2 = v.squared_norm();
        if v2 <= P::Scalar::epsilon() {
            return zero;
        }
        let a2 = a.squared_norm();
        let dot = PointDot::dot(&v, &a);
        let mut num = v2 * a2 - dot * dot;
        if num < zero {
            num = zero;
        }
        let denom = v2 * v2.sqrt();
        if denom <= P::Scalar::epsilon() {
            zero
        } else {
            num.sqrt() / denom
        }
    }

    /// Return the principal normal direction at `t`.
    ///
    /// Returns `None` if the velocity is zero or curvature is undefined.
    pub fn normal(&self, t: P::Scalar) -> Option<P>
    where
        P: PointNorm + PointDot,
        [(); N - 1]: Sized,
        [(); N - 2]: Sized,
        [(); (N - 1) - 1]: Sized,
    {
        let first = self.derivative();
        let v = first.eval(t);
        let v2 = v.squared_norm();
        if v2 <= P::Scalar::epsilon() {
            return None;
        }
        let a = if N < 3 {
            v - v
        } else {
            first.derivative().eval(t)
        };
        let dot = PointDot::dot(&v, &a);
        let a_perp = a - v * (dot / v2);
        let n2 = a_perp.squared_norm();
        if n2 <= P::Scalar::epsilon() {
            return None;
        }
        let one = <P::Scalar as NumCast>::from(1.0).unwrap();
        Some(a_perp * (one / n2.sqrt()))
    }

    /// Approximates the arc length of the curve by flattening it with straight line segments.
    /// This works quite well, at ~32 segments it should already provide an error in the decimal places
    /// The accuracy gain falls off with more steps so this approximation is unfeasable if desired accuracy is greater than 1-2 decimal places
    pub fn arclen(&self, nsteps: usize) -> P::Scalar
    where
        P: PointNorm,
    {
        let nsteps = nsteps.max(1);
        let nsteps_scalar = <P::Scalar as NumCast>::from(nsteps as f64).unwrap();
        let one = <P::Scalar as NumCast>::from(1.0).unwrap();
        let stepsize = one / nsteps_scalar;
        let mut arclen: P::Scalar = <P::Scalar as NumCast>::from(0.0).unwrap();
        for i in 0..nsteps {
            let t0 = <P::Scalar as NumCast>::from(i as f64).unwrap() / nsteps_scalar;
            let t1 = if i + 1 == nsteps {
                one
            } else {
                t0 + stepsize
            };
            let p1 = self.eval(t0);
            let p2 = self.eval(t1);

            arclen = arclen + (p1 - p2).squared_norm().sqrt();
        }
        arclen
    }

    fn arclen_partial(&self, t: P::Scalar, nsteps: usize) -> P::Scalar
    where
        P: PointNorm,
    {
        let zero = <P::Scalar as NumCast>::from(0.0).unwrap();
        let one = <P::Scalar as NumCast>::from(1.0).unwrap();
        let t = t.clamp(zero, one);
        if t <= zero {
            return zero;
        }

        let nsteps = nsteps.max(1);
        let nsteps_scalar = <P::Scalar as NumCast>::from(nsteps as f64).unwrap();
        let dt = t / nsteps_scalar;
        let mut arclen = zero;
        let mut prev = self.eval(zero);

        for i in 1..=nsteps {
            let i_scalar = <P::Scalar as NumCast>::from(i as f64).unwrap();
            let ti = if i == nsteps { t } else { dt * i_scalar };
            let p = self.eval(ti);
            arclen = arclen + (p - prev).squared_norm().sqrt();
            prev = p;
        }

        arclen
    }

    /// Find parameter values where the derivative crosses zero for a given axis.
    pub fn derivative_roots(&self, axis: usize) -> ArrayVec<[P::Scalar; N - 1]>
    where
        P: PointIndex,
        [(); N - 1]: Sized,
        [P::Scalar; N - 1]: tinyvec::Array<Item = P::Scalar>,
    {
        self.derivative_roots_with_tolerance(axis, <P::Scalar as NumCast>::from(ROOT_TOLERANCE).unwrap())
    }

    /// Find parameter values where the derivative crosses zero for a given axis using a tolerance.
    pub fn derivative_roots_with_tolerance(
        &self,
        axis: usize,
        tolerance: P::Scalar,
    ) -> ArrayVec<[P::Scalar; N - 1]>
    where
        P: PointIndex,
        [(); N - 1]: Sized,
        [P::Scalar; N - 1]: tinyvec::Array<Item = P::Scalar>,
    {
        if N < 2 {
            return ArrayVec::new();
        }
        let control = self.derivative_axis_control_points(axis);
        let mut roots = ArrayVec::<[P::Scalar; N - 1]>::new();
        let t0 = <P::Scalar as NumCast>::from(0.0).unwrap();
        let t1 = <P::Scalar as NumCast>::from(1.0).unwrap();
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
        P: PointIndex,
        [(); N - 1]: Sized,
        [P::Scalar; N - 1]: tinyvec::Array<Item = P::Scalar>,
    {
        let tolerance = <P::Scalar as NumCast>::from(ROOT_TOLERANCE).unwrap();
        let mut bounds = [(<P::Scalar as NumCast>::from(0.0).unwrap(), <P::Scalar as NumCast>::from(0.0).unwrap()); P::DIM];
        let zero = <P::Scalar as NumCast>::from(0.0).unwrap();
        let one = <P::Scalar as NumCast>::from(1.0).unwrap();

        for dim in 0..P::DIM {
            let mut min = self.control_points[0][dim];
            let mut max = min;
            let end = self.control_points[N - 1][dim];
            if end < min {
                min = end;
            }
            if end > max {
                max = end;
            }

            let roots = self.derivative_roots_with_tolerance(dim, tolerance);
            for &t in roots.iter() {
                if t > zero && t < one {
                    let value = self.eval(t)[dim];
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
        P: PointIndex,
        [(); N - 1]: Sized,
    {
        FindRoot::root_newton_axis(self, value, axis, start, eps, max_iter)
    }

    fn derivative_axis_control_points(&self, axis: usize) -> [P::Scalar; N - 1]
    where
        P: PointIndex,
        [(); N - 1]: Sized,
    {
        let scale = <P::Scalar as NumCast>::from((N - 1) as f64).unwrap();
        core::array::from_fn(|i| {
            (self.control_points[i + 1][axis] - self.control_points[i][axis]) * scale
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

        if min > <P::Scalar as NumCast>::from(0.0).unwrap()
            || max < <P::Scalar as NumCast>::from(0.0).unwrap()
        {
            return;
        }

        if depth >= MAX_ROOT_DEPTH || (max - min).abs() <= tolerance {
            let denom = control[M - 1] - control[0];
            let root = if denom.abs() > tolerance {
                t0 + (<P::Scalar as NumCast>::from(0.0).unwrap() - control[0]) * (t1 - t0) / denom
            } else {
                (t0 + t1) * <P::Scalar as NumCast>::from(0.5).unwrap()
            };
            if results.len() < results.capacity() {
                results.push(root);
            }
            return;
        }

        let (left, right) = Self::subdivide_1d(control);
        let mid = (t0 + t1) * <P::Scalar as NumCast>::from(0.5).unwrap();
        Self::find_roots_1d(left, t0, mid, depth + 1, tolerance, results);
        Self::find_roots_1d(right, mid, t1, depth + 1, tolerance, results);
    }

    fn subdivide_1d<const M: usize>(
        control: [P::Scalar; M],
    ) -> ([P::Scalar; M], [P::Scalar; M]) {
        let half = <P::Scalar as NumCast>::from(0.5).unwrap();
        let mut left = core::array::from_fn(|_| <P::Scalar as NumCast>::from(0.0).unwrap());
        let mut right = core::array::from_fn(|_| <P::Scalar as NumCast>::from(0.0).unwrap());
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
    P: PointIndex,
    [(); N - 1]: Sized,
{
    fn parameter_domain(&self) -> (P::Scalar, P::Scalar) {
        (
            <P::Scalar as NumCast>::from(0.0).unwrap(),
            <P::Scalar as NumCast>::from(1.0).unwrap(),
        )
    }

    fn axis_value(&self, t: P::Scalar, axis: usize) -> Result<P::Scalar, RootFindingError> {
        Ok(self.eval(t)[axis])
    }

    fn axis_derivative(&self, t: P::Scalar, axis: usize) -> Result<P::Scalar, RootFindingError> {
        Ok(self.derivative().eval(t)[axis])
    }
}

#[cfg(test)]
mod tests {
    use super::CubicBezier;
    use super::QuadraticBezier;
    use super::*;
    use crate::{PointN, EPSILON};

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
        assert!(err_start.squared_norm() < EPSILON);

        let end = curve.eval(1.0);
        let err_end = end - points[points.len() - 1];
        assert!(err_end.squared_norm() < EPSILON);
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
    fn arclen_line_approx() {
        let bezier: Bezier<PointN<f64, 2>, 2> = Bezier::new([
            PointN::new([0.0, 0.0]),
            PointN::new([10.0, 0.0]),
        ]);
        let expected = (bezier.eval(1.0) - bezier.eval(0.0)).squared_norm().sqrt();
        let length = bezier.arclen(32);
        let tolerance = expected * 0.05;
        assert!((length - expected).abs() <= tolerance);
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
            assert!(err.squared_norm() < EPSILON);
            // check the right part of the split curve
            err = bezier.eval((t * 0.5) + 0.5) - right.eval(t);
            assert!(err.squared_norm() < EPSILON);
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
            assert!(err.squared_norm() < EPSILON);
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
            assert!(err.squared_norm() < EPSILON);
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
        let v0 = bezier.eval(t0)[0];
        let v1 = bezier.eval(t1)[0];
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

    #[test]
    fn bezier_api_parity() {
        let bezier = Bezier::<PointN<f64, 2>, 2>::new([
            PointN::new([0f64, 0f64]),
            PointN::new([2f64, 0f64]),
        ]);

        assert_eq!(bezier.start(), bezier.eval(0.0));
        assert_eq!(bezier.end(), bezier.eval(1.0));

        let reversed = bezier.reverse();
        let p0 = bezier.eval(0.25);
        let p1 = reversed.eval(0.75);
        assert!((p0 - p1).squared_norm() < EPSILON);

        let tangent = bezier.tangent(0.3);
        assert!((tangent[0] - 1.0).abs() < EPSILON);
        assert!(tangent[1].abs() < EPSILON);

        let curvature = bezier.curvature(0.3);
        assert!(curvature.abs() < EPSILON);
        assert!(bezier.normal(0.3).is_none());

        let t = bezier.t_at_length(1.0);
        assert!((t - 0.5).abs() < 1e-3);

        let p = bezier.point_at_length(1.0);
        assert!((p - PointN::new([1.0, 0.0])).squared_norm() < EPSILON);
    }
}
