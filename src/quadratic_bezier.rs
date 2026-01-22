//! Quadratic Bezier curve specialization.

use super::{ArrayVec, LineSegment, Point, PointDot, PointIndex, PointNorm};
use num_traits::{Float, NumCast};

const DEFAULT_DISTANCE_STEPS: usize = 64;
const DEFAULT_LENGTH_STEPS: usize = 64;
const MAX_LENGTH_ITERS: usize = 32;
const LOCAL_SEARCH_ITERS: usize = 16;

/// Quadratic Bezier curve defined by start, control, and end points.
///
/// Methods that need component access or norms add `PointIndex`/`PointNorm`
/// bounds as required.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct QuadraticBezier<P: Point> {
    pub(crate) start: P,
    pub(crate) ctrl: P,
    pub(crate) end: P,
}

impl<P: Point> QuadraticBezier<P>
where
    P: Point,
{
    /// Creates a new instance of QuadraticBezier from the given control points
    pub fn new(start: P, ctrl: P, end: P) -> Self {
        QuadraticBezier { start, ctrl, end }
    }

    /// Evaluates the quadratic bezier curve at 't' using direct evaluation, which may not be numerically stable
    pub fn eval(&self, t: P::Scalar) -> P {
        let t2 = t * t;
        let one = <P::Scalar as NumCast>::from(1.0).unwrap();
        let one_t = one - t;
        let one_t2 = one_t * one_t;

        self.start * one_t2
            + self.ctrl * <P::Scalar as NumCast>::from(2.0).unwrap() * one_t * t
            + self.end * t2
    }

    /// Evaluates the quadratic bezier curve at t using the numerically stable De Casteljau algorithm
    pub fn eval_casteljau(&self, t: P::Scalar) -> P {
        // unrolled de casteljau algorithm
        // _1ab is the first iteration from first (a) to second (b) control point and so on
        let ctrl_1ab = self.start + (self.ctrl - self.start) * t;
        let ctrl_1bc = self.ctrl + (self.end - self.ctrl) * t;
        // second iteration, return final point on the curve ctrl_2ab
        ctrl_1ab + (ctrl_1bc - ctrl_1ab) * t
    }

    /// Return the control points array.
    pub fn control_points(&self) -> [P; 3] {
        [self.start, self.ctrl, self.end]
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
        QuadraticBezier {
            start: self.end,
            ctrl: self.ctrl,
            end: self.start,
        }
    }

    /// Split the curve at `t` into two sub-curves.
    pub fn split(&self, t: P::Scalar) -> (Self, Self) {
        // unrolled de casteljau algorithm
        // _1ab is the first iteration from first (a) to second (b) control point and so on
        let ctrl_1ab = self.start + (self.ctrl - self.start) * t;
        let ctrl_1bc = self.ctrl + (self.end - self.ctrl) * t;
        // second iteration
        let ctrl_2ab = ctrl_1ab + (ctrl_1bc - ctrl_1ab) * t;

        (
            QuadraticBezier {
                start: self.start,
                ctrl: ctrl_1ab,
                end: ctrl_2ab,
            },
            QuadraticBezier {
                start: ctrl_2ab,
                ctrl: ctrl_1bc,
                end: self.end,
            },
        )
    }

    /// Return the derivative curve.
    /// The derivative is also a bezier curve but of degree n-1.
    /// In the case of a quadratic derivative it is just a line segment
    /// which also implementes eval(), as it is just a linear bezier curve.
    pub fn derivative(&self) -> LineSegment<P> {
        let two = <P::Scalar as NumCast>::from(2.0).unwrap();
        LineSegment {
            start: (self.ctrl - self.start) * two,
            end: (self.end - self.ctrl) * two,
        }
    }

    /// Return the unit tangent direction at `t`.
    pub fn tangent(&self, t: P::Scalar) -> P
    where
        P: PointNorm,
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
    {
        let zero = <P::Scalar as NumCast>::from(0.0).unwrap();
        let v = self.derivative().eval(t);
        let a = self.derivative().derivative();
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
    {
        let v = self.derivative().eval(t);
        let v2 = v.squared_norm();
        if v2 <= P::Scalar::epsilon() {
            return None;
        }
        let a = self.derivative().derivative();
        let dot = PointDot::dot(&v, &a);
        let a_perp = a - v * (dot / v2);
        let n2 = a_perp.squared_norm();
        if n2 <= P::Scalar::epsilon() {
            return None;
        }
        let one = <P::Scalar as NumCast>::from(1.0).unwrap();
        Some(a_perp * (one / n2.sqrt()))
    }

    /// Direct Derivative - Sample the axis coordinate at 'axis' of the curve's derivative at t
    /// without creating a new curve. This is a convenience function for `.derivative().eval(t)[n]`.
    /// Parameters:
    ///   t: the sampling parameter on the curve interval [0..1]
    ///   axis: the index of the coordinate axis [0..N]
    /// Returns:
    ///   Scalar value of the points own type type F  
    /// May be deprecated in the future.  
    /// This function can cause out of bounds panic when axis is larger than dimension of P
    pub fn dd(&self, t: P::Scalar, axis: usize) -> P::Scalar
    where
        P: PointIndex,
    {
        let two = <P::Scalar as NumCast>::from(2.0).unwrap();
        let four = <P::Scalar as NumCast>::from(4.0).unwrap();
        let c0 = t * two - two;
        let c1 = two - four * t;
        let c2 = two * t;

        self.start[axis] * c0 + self.ctrl[axis] * c1 + self.end[axis] * c2
    }

    /// Approximates the arc length of the curve by flattening it with straight line segments.
    /// This works quite well, at ~32 segments it should already provide an error < 0.5
    /// Remember arclen also works by linear approximation, not the integral, so we have to accept error!
    /// This approximation is unfeasable if desired accuracy is greater than 2 decimal places
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
            let t1 = if i + 1 == nsteps { one } else { t0 + stepsize };
            let p1 = self.eval_casteljau(t0);
            let p2 = self.eval_casteljau(t1);

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
        let mut prev = self.eval_casteljau(zero);

        for i in 1..=nsteps {
            let i_scalar = <P::Scalar as NumCast>::from(i as f64).unwrap();
            let ti = if i == nsteps { t } else { dt * i_scalar };
            let p = self.eval_casteljau(ti);
            arclen = arclen + (p - prev).squared_norm().sqrt();
            prev = p;
        }

        arclen
    }

    /// Solve for the roots of the polynomial at^2 + bt + c
    /// Returns an ArrayVec of roots in the order
    /// needs to be called for x and y components separately
    pub(crate) fn real_roots(
        &self,
        a: P::Scalar,
        b: P::Scalar,
        c: P::Scalar,
    ) -> ArrayVec<[P::Scalar; 2]>
    where
        [P::Scalar; 2]: tinyvec::Array<Item = P::Scalar>,
    {
        let mut result = ArrayVec::new();
        let eps = P::Scalar::epsilon();
        let zero = <P::Scalar as NumCast>::from(0.0).unwrap();
        let two = <P::Scalar as NumCast>::from(2.0).unwrap();
        let four = <P::Scalar as NumCast>::from(4.0).unwrap();

        // check if can be handled below quadratic order
        if a.abs() < eps {
            if b.abs() < eps {
                // no solutions
                return result;
            }
            // is linear equation
            result.push(-c / b);
            return result;
        }
        // is quadratic equation
        let delta = b * b - a * c * four;
        if delta > zero {
            let sqrt_delta = delta.sqrt();
            result.push((-b - sqrt_delta) / (a * two));
            result.push((-b + sqrt_delta) / (a * two));
        } else if delta.abs() < eps {
            result.push(-b / (a * two));
        }
        result
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
        self.eval_casteljau(self.t_at_length_approx(s, nsteps))
    }

    /// Evaluate the point at arc length `s` using a default resolution.
    pub fn point_at_length(&self, s: P::Scalar) -> P
    where
        P: PointNorm,
    {
        self.point_at_length_approx(s, DEFAULT_LENGTH_STEPS)
    }

    /// Returns the line segment formed by the curve's start and endpoint
    pub fn baseline(&self) -> LineSegment<P> {
        LineSegment {
            start: self.start,
            end: self.end,
        }
    }

    /// Checks if, given some tolerance, the curve can be considered equal to a line segment
    pub fn is_linear(&self, tolerance: P::Scalar) -> bool
    where
        P: PointNorm,
    {
        // if start and end are (nearly) the same
        // TODO using squared norm vs machine epsilon OK?
        if (self.start - self.end).squared_norm() < P::Scalar::epsilon() {
            return false;
        }
        // else check if ctrl points lie on baseline i.e. all points are colinear
        self.are_points_colinear(tolerance)
    }

    /// Determines if, given some tolerance, all of the control points are colinear
    /// This private function is wrapped publically by is_linear()
    fn are_points_colinear(&self, tolerance: P::Scalar) -> bool
    where
        P: PointNorm,
    {
        let line = self.baseline();
        line.distance_to_point(self.ctrl) <= tolerance
    }

    /// Determines if, given some tolerance, the control points of the curve can be considered equal.
    /// If true, the curve is just a singular point
    pub fn is_a_point(&self, tolerance: P::Scalar) -> bool
    where
        P: PointNorm,
    {
        let tolerance_squared = tolerance * tolerance;
        // Use <= so that tolerance can be zero.
        (self.start - self.end).squared_norm() <= tolerance_squared
            && (self.start - self.ctrl).squared_norm() <= tolerance_squared
    }

    /// Solves the quadratic bezier function given a particular coordinate axis value
    /// by solving the roots for the axis functions
    /// Parameters:
    /// value: the coordinate value on the particular axis
    /// axis: the index of the axis
    /// Returns those roots of the function that are in the interval [0.0, 1.0].
    #[allow(dead_code)]
    fn solve_t_for_axis(&self, value: P::Scalar, axis: usize) -> ArrayVec<[P::Scalar; 3]>
    where
        P: PointIndex + PointNorm,
        [P::Scalar; 2]: tinyvec::Array<Item = P::Scalar>,
        [P::Scalar; 3]: tinyvec::Array<Item = P::Scalar>,
    {
        let mut result = ArrayVec::new();
        if self.is_a_point(P::Scalar::epsilon())
            || (self.are_points_colinear(<P::Scalar as NumCast>::from(0.0).unwrap())
                && (self.start - self.end).squared_norm() < P::Scalar::epsilon())
        {
            return result;
        }
        // these are just the x or y components of the points
        let a = self.start[axis]
            + self.ctrl[axis] * <P::Scalar as NumCast>::from(-2.0).unwrap()
            + self.end[axis];
        let b = self.start[axis] * <P::Scalar as NumCast>::from(-2.0).unwrap()
            + self.ctrl[axis] * <P::Scalar as NumCast>::from(2.0).unwrap();
        let c = self.start[axis] - value;

        let roots = self.real_roots(a, b, c);
        let zero = <P::Scalar as NumCast>::from(0.0).unwrap();
        let one = <P::Scalar as NumCast>::from(1.0).unwrap();
        for &root in roots.iter() {
            if root > zero && root < one {
                result.push(root);
            }
        }

        result
    }

    /// Return the bounding box of the curve as an array of (min, max) tuples for each dimension (its index)
    pub fn bounding_box(&self) -> [(P::Scalar, P::Scalar); P::DIM]
    where
        P: PointIndex,
        [P::Scalar; 1]: tinyvec::Array<Item = P::Scalar>,
        [P::Scalar; 3]: tinyvec::Array<Item = P::Scalar>,
    {
        let zero = <P::Scalar as NumCast>::from(0.0).unwrap();
        let one = <P::Scalar as NumCast>::from(1.0).unwrap();
        let mut bounds = [(zero, zero); P::DIM];
        let derivative = self.derivative();
        // calculate coefficients for the derivative as a function of t: at + b
        // po: [1, -1]
        // p1: [0,  1]
        //      b   a
        let a = derivative.start * <P::Scalar as NumCast>::from(-1.0).unwrap() + derivative.end;
        let b = derivative.start;

        for dim in 0..P::DIM {
            // calculate roots for t over x axis and plug them into the bezier function
            //  to get x,y values (make vec 2 bigger for t=0,t=1 values)
            let mut extrema: ArrayVec<[P::Scalar; 3]> = ArrayVec::new();
            extrema.extend(derivative.root(a[dim], b[dim]).iter().copied());
            // only retain roots for which t is in [0..1]
            extrema.retain(|root| -> bool { *root > zero && *root < one });
            // evaluates roots in original function
            for t in extrema.iter_mut() {
                *t = self.eval_casteljau(*t)[dim];
            }
            // add y-values for start and end point as candidates
            extrema.push(self.start[dim]);
            extrema.push(self.end[dim]);
            // sort to get min and max values for bounding box
            extrema.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap());

            // determine xmin, xmax, ymin, ymax, from the set {B(xroots), B(yroots), B(0), B(1)}
            // (Intermediate control points can't form a boundary)
            // unwrap() is ok as it always at least contains the endpoints
            bounds[dim] = (extrema[0], *extrema.last().unwrap());
        }
        bounds
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{EPSILON, PointN};

    #[test]
    fn eval_endpoints() {
        let bezier = QuadraticBezier::new(
            PointN::new([0f64, 0f64]),
            PointN::new([1f64, 1f64]),
            PointN::new([2f64, 0f64]),
        );

        assert_eq!(bezier.eval(0.0), PointN::new([0.0, 0.0]));
        assert_eq!(bezier.eval(1.0), PointN::new([2.0, 0.0]));
    }

    #[test]
    fn arclen_line_approx() {
        let bezier = QuadraticBezier::new(
            PointN::new([0f64, 0f64]),
            PointN::new([1f64, 0f64]),
            PointN::new([2f64, 0f64]),
        );
        let length = bezier.arclen(32);
        assert!((length - 2.0).abs() < 1e-3);
    }

    #[test]
    fn eval_equivalence() {
        // all eval methods should be approximately equivalent for well defined test cases
        // and not equivalent where numerical stability becomes an issue for normal eval
        let bezier = QuadraticBezier::new(
            PointN::new([0f64, 1.77f64]),
            PointN::new([4.3f64, 3f64]),
            PointN::new([3.2f64, -4f64]),
        );

        let nsteps: usize = 1000;
        for t in 0..=nsteps {
            let t = t as f64 * 1f64 / (nsteps as f64);
            let p1 = bezier.eval(t);
            let p2 = bezier.eval_casteljau(t);
            let err = p2 - p1;
            assert!(err.squared_norm() < EPSILON);
        }
    }

    #[test]
    fn split_equivalence() {
        // chose some arbitrary control points and construct a cubic bezier
        let bezier = QuadraticBezier {
            start: PointN::new([0f64, 1.77f64]),
            ctrl: PointN::new([4.3f64, 3f64]),
            end: PointN::new([3.2f64, -4f64]),
        };
        // split it at an arbitrary point
        let at = 0.5;
        let (left, right) = bezier.split(at);
        // compare left and right subcurves with parent curve
        // take the difference of the two points which must not exceed the absolute error
        let nsteps: usize = 1000;
        for t in 0..=nsteps {
            let t = t as f64 * 1f64 / (nsteps as f64);
            // left
            let mut err = bezier.eval(t / 2.0) - left.eval(t);
            assert!(err.squared_norm() < EPSILON);
            // right
            err = bezier.eval((t * 0.5) + 0.5) - right.eval(t);
            assert!(err.squared_norm() < EPSILON);
        }
    }

    #[test]
    fn bounding_box_contains() {
        // check if bounding box for a curve contains all points (with some approximation error)
        let bezier = QuadraticBezier {
            start: PointN::new([0f64, 1.77f64]),
            ctrl: PointN::new([4.3f64, -3f64]),
            end: PointN::new([3.2f64, 4f64]),
        };

        let bounds = bezier.bounding_box();

        let max_err = 1e-2;

        let nsteps: usize = 100;
        for t in 0..=nsteps {
            let t = t as f64 * 1f64 / (nsteps as f64);
            let p = bezier.eval_casteljau(t);

            for (idx, axis) in p.into_iter().enumerate() {
                assert!((axis >= (bounds[idx].0 - max_err)) && (axis <= (bounds[idx].1 + max_err)))
            }
        }
    }

    #[test]
    fn distance_to_point() {
        // degree 3, 4 control points => 4+3+1=8 knots
        let curve = QuadraticBezier {
            start: PointN::new([0f64, 1.77f64]),
            ctrl: PointN::new([4.3f64, 3f64]),
            end: PointN::new([3.2f64, -4f64]),
        };
        assert!(
            curve.distance_to_point(PointN::new([-5.1, -5.6]))
                > curve.distance_to_point(PointN::new([5.1, 5.6]))
        );
    }

    #[test]
    fn quadratic_api_parity() {
        let bezier = QuadraticBezier::new(
            PointN::new([0f64, 0f64]),
            PointN::new([1f64, 0f64]),
            PointN::new([2f64, 0f64]),
        );

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
        assert!((t - 0.5).abs() < 1e-6);

        let p = bezier.point_at_length(1.0);
        assert!((p - PointN::new([1.0, 0.0])).squared_norm() < EPSILON);
    }

    #[test]
    fn quadratic_curvature_nonzero() {
        let bezier = QuadraticBezier::new(
            PointN::new([0f64, 0f64]),
            PointN::new([1f64, 1f64]),
            PointN::new([2f64, 0f64]),
        );

        let curvature = bezier.curvature(0.5);
        assert!(curvature > 0.0);

        let normal = bezier.normal(0.5).unwrap();
        assert!((normal.squared_norm() - 1.0).abs() < 1e-6);
    }
}
