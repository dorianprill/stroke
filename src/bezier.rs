//! Const-generic Bezier curves of arbitrary degree.

use core::iter::IntoIterator;
use core::slice;

use super::*;
use crate::arc_length::integrate;
use crate::error::{CurveConstructionError, validate_control_points};
use crate::find_root::FindRoot;
use crate::point::{Point, PointIndex, PointNorm};
use crate::roots::RootFindingError;
use num_traits::{Float, NumCast};

const MAX_ROOT_DEPTH: usize = 64;
const DEFAULT_DISTANCE_STEPS: usize = 64;
const LOCAL_SEARCH_ITERS: usize = 16;

/// General implementation of a Bezier curve of arbitrary degree (`N - 1`).
///
/// The curve is defined by its control points. Points on the curve can be
/// evaluated with an interpolation parameter `t` in `[0, 1]` using `eval()`.
/// Methods that need component access or norms add `PointIndex`/`PointNorm`
/// bounds as required.
///
/// # Parameters
/// - `P`: point type implementing [`Point`].
/// - `N`: number of control points.
#[derive(Clone, Copy)]
pub struct Bezier<P, const N: usize>
where
    P: Point,
{
    /// Control points which define the curve and hence its degree
    pub(crate) control_points: [P; N],
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
    ///
    /// # Errors
    ///
    /// Returns [`CurveConstructionError::NoControlPoints`] when `N` is zero, or
    /// [`CurveConstructionError::NonFiniteControlPoint`] with the index of the
    /// first control point containing a non-finite coordinate.
    pub fn new(control_points: [P; N]) -> Result<Bezier<P, { N }>, CurveConstructionError> {
        if control_points.is_empty() {
            return Err(CurveConstructionError::NoControlPoints);
        }
        validate_control_points(&control_points)?;
        Ok(Bezier { control_points })
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

    /// Evaluate a point on the curve at point 't' which should be in the interval \[0, 1\] (unchecked!)
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

    /// Return parameter `t` at arc length `s`, clamped to the curve's endpoints.
    /// `accuracy` is an absolute length tolerance, subject to floating-point roundoff.
    pub fn t_at_length(
        &self,
        s: P::Scalar,
        accuracy: P::Scalar,
    ) -> Result<P::Scalar, ArcLengthError>
    where
        P: PointNorm,
    {
        let total = self.arclen(accuracy)?;
        crate::arc_length::inverse(s, accuracy, total, |a, b, tolerance| {
            self.arclen_between(a, b, tolerance)
        })
    }

    /// Evaluate the point at arc length `s` using adaptive quadrature.
    pub fn point_at_length(&self, s: P::Scalar, accuracy: P::Scalar) -> Result<P, ArcLengthError>
    where
        P: PointNorm,
    {
        Ok(self.eval(self.t_at_length(s, accuracy)?))
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
                casteljau_points[j] = casteljau_points[j] * (one - t) + casteljau_points[j + 1] * t;
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
        let new_points: [P; N - 1] =
            core::array::from_fn(|i| (self.control_points[i + 1] - self.control_points[i]) * scale);
        Bezier {
            control_points: new_points,
        }
    }

    /// Return the unit tangent direction at `t`.
    pub fn tangent(&self, t: P::Scalar) -> P
    where
        P: PointNorm,
        [(); N - 1]: Sized,
    {
        let one = <P::Scalar as NumCast>::from(1.0).unwrap();
        let zero = <P::Scalar as NumCast>::from(0.0).unwrap();
        let dir = self.velocity(t);
        let len = dir.norm();
        if len == zero { dir } else { dir * (one / len) }
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
            // N < 3 has no second derivative; use a typed zero vector.
            v * zero
        } else {
            first.derivative().eval(t)
        };
        let v2 = v.squared_norm();
        if v2 == zero {
            return zero;
        }
        let a2 = a.squared_norm();
        let dot = PointDot::dot(&v, &a);
        let mut num = v2 * a2 - dot * dot;
        if num < zero {
            num = zero;
        }
        let denom = v2 * v2.sqrt();
        if denom == zero {
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
        let zero = <P::Scalar as NumCast>::from(0.0).unwrap();
        let first = self.derivative();
        let v = first.eval(t);
        let v2 = v.squared_norm();
        if v2 == zero {
            return None;
        }
        let a = if N < 3 {
            // N < 3 has no second derivative; use a typed zero vector.
            let zero = <P::Scalar as NumCast>::from(0.0).unwrap();
            v * zero
        } else {
            first.derivative().eval(t)
        };
        let dot = PointDot::dot(&v, &a);
        let a_perp = a - v * (dot / v2);
        let n2 = a_perp.squared_norm();
        if n2 == zero {
            return None;
        }
        let one = <P::Scalar as NumCast>::from(1.0).unwrap();
        Some(a_perp * (one / n2.sqrt()))
    }

    /// Approximate the arc length by flattening the curve into `segments` lines.
    pub fn arclen_polyline(&self, segments: usize) -> P::Scalar
    where
        P: PointNorm,
    {
        let segments = segments.max(1);
        let count = <P::Scalar as NumCast>::from(segments).unwrap();
        let zero = <P::Scalar as NumCast>::from(0.0).unwrap();
        let one = <P::Scalar as NumCast>::from(1.0).unwrap();
        let mut length = zero;
        let mut previous = self.control_points[0];
        for i in 1..=segments {
            let t = if i == segments {
                one
            } else {
                <P::Scalar as NumCast>::from(i).unwrap() / count
            };
            let point = self.eval(t);
            length = length + (point - previous).norm();
            previous = point;
        }
        length
    }

    /// Return the arc length using adaptive Gauss-Legendre quadrature.
    pub fn arclen(&self, accuracy: P::Scalar) -> Result<P::Scalar, ArcLengthError>
    where
        P: PointNorm,
    {
        self.arclen_to(<P::Scalar as NumCast>::from(1.0).unwrap(), accuracy)
    }

    fn arclen_to(&self, t: P::Scalar, accuracy: P::Scalar) -> Result<P::Scalar, ArcLengthError>
    where
        P: PointNorm,
    {
        self.arclen_between(<P::Scalar as NumCast>::from(0.0).unwrap(), t, accuracy)
    }

    /// Integrate arc length over `start..end`, where `0 <= start <= end <= 1`.
    /// The positive finite `accuracy` is an absolute length tolerance.
    pub fn arclen_between(
        &self,
        start: P::Scalar,
        end: P::Scalar,
        accuracy: P::Scalar,
    ) -> Result<P::Scalar, ArcLengthError>
    where
        P: PointNorm,
    {
        crate::arc_length::validate_accuracy(accuracy)?;
        let zero = <P::Scalar as NumCast>::from(0.0).unwrap();
        let one = <P::Scalar as NumCast>::from(1.0).unwrap();
        if !start.is_finite() || !end.is_finite() || start < zero || end > one || start > end {
            return Err(ArcLengthError::InvalidParameter);
        }
        if N < 2 {
            return Ok(zero);
        }
        if N == 2 {
            return LineSegment {
                start: self.control_points[0],
                end: self.control_points[1],
            }
            .arclen_between(start, end, accuracy);
        }
        // Every zero-velocity point is also a root of any nonzero velocity
        // component. Split at one such component's roots to expose cusps.
        let mut selected = [zero; N];
        let mut largest = zero;
        for axis in 0..P::DIM {
            let control = self.elevated_velocity_axis(axis);
            if control.iter().any(|value| !value.is_finite()) {
                return Err(ArcLengthError::NonFiniteComputation);
            }
            let scale = control
                .iter()
                .fold(zero, |scale, value| scale.max(value.abs()));
            if scale > largest {
                largest = scale;
                selected = control;
            }
        }
        if largest == zero {
            return Ok(zero);
        }
        // Share the accuracy budget across all roots plus the final interval.
        let pieces = <P::Scalar as NumCast>::from(N + 1).unwrap();
        let tolerance =
            (accuracy / pieces).max(P::Scalar::min_positive_value() * P::Scalar::epsilon());
        let root_tolerance = P::Scalar::epsilon() * <P::Scalar as NumCast>::from(8_u8).unwrap();
        let mut previous = start;
        let mut sum = zero;
        let mut result = Ok(());
        let add_interval = |a, b| {
            integrate(a, b, tolerance, |parameter| {
                let speed = self.velocity(parameter).norm();
                if speed.is_finite() {
                    Ok(speed)
                } else {
                    Err(ArcLengthError::NonFiniteComputation)
                }
            })
        };
        crate::numeric::bernstein_roots(selected, root_tolerance, |root| {
            if root > previous && root < end && result.is_ok() {
                match add_interval(previous, root) {
                    Ok(length) => sum = sum + length,
                    Err(error) => result = Err(error),
                }
                previous = root;
            }
        });
        result?;
        sum = sum + add_interval(previous, end)?;
        if sum.is_finite() {
            Ok(sum)
        } else {
            Err(ArcLengthError::NonFiniteComputation)
        }
    }

    // Degree-elevate the unscaled derivative to N coefficients. Keeping its
    // size N avoids leaking N-1 const-expression bounds into the length API.
    fn elevated_velocity_axis(&self, axis: usize) -> [P::Scalar; N]
    where
        P: PointIndex,
    {
        let zero = <P::Scalar as NumCast>::from(0.0).unwrap();
        let one = <P::Scalar as NumCast>::from(1.0).unwrap();
        let degree = <P::Scalar as NumCast>::from(N - 1).unwrap();
        core::array::from_fn(|i| {
            let weight = <P::Scalar as NumCast>::from(i).unwrap() / degree;
            let left = if i > 0 {
                self.control_points[i][axis] - self.control_points[i - 1][axis]
            } else {
                zero
            };
            let right = if i + 1 < N {
                self.control_points[i + 1][axis] - self.control_points[i][axis]
            } else {
                zero
            };
            left * weight + right * (one - weight)
        })
    }

    fn velocity(&self, t: P::Scalar) -> P {
        let zero = <P::Scalar as NumCast>::from(0.0).unwrap();
        if N < 2 {
            return self.control_points[0] * zero;
        }
        let degree = <P::Scalar as NumCast>::from(N - 1).unwrap();
        let one = <P::Scalar as NumCast>::from(1.0).unwrap();
        let mut points = self.control_points;
        for (i, point) in points.iter_mut().enumerate().take(N - 1) {
            *point = (self.control_points[i + 1] - self.control_points[i]) * degree;
        }
        for level in 1..N - 1 {
            for i in 0..N - 1 - level {
                points[i] = points[i] * (one - t) + points[i + 1] * t;
            }
        }
        points[0]
    }

    /// Find isolated derivative zeros in `(0, 1)` for a coordinate axis.
    /// Tangential roots are included; an identically constant axis returns no roots.
    pub fn derivative_roots(&self, axis: usize) -> ArrayVec<[P::Scalar; N - 1]>
    where
        P: PointIndex,
        [(); N - 1]: Sized,
        [P::Scalar; N - 1]: tinyvec::Array<Item = P::Scalar>,
    {
        self.derivative_roots_with_tolerance(
            axis,
            P::Scalar::epsilon() * <P::Scalar as NumCast>::from(64_u8).unwrap(),
        )
    }

    /// Find isolated derivative zeros using a tolerance in parameter space.
    /// The absolute tolerance is floored at eight scalar epsilons.
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
        let tolerance = tolerance
            .abs()
            .max(P::Scalar::epsilon() * <P::Scalar as NumCast>::from(8_u8).unwrap());
        let control = self.derivative_axis_control_points(axis);
        let mut roots = ArrayVec::<[P::Scalar; N - 1]>::new();
        let t0 = <P::Scalar as NumCast>::from(0.0).unwrap();
        let t1 = <P::Scalar as NumCast>::from(1.0).unwrap();
        crate::numeric::bernstein_roots(control, tolerance, |root| {
            if root > t0 && root < t1 {
                Self::push_root(root, tolerance, &mut roots);
            }
        });

        roots
    }

    /// Return a tight axis-aligned box of the curve on `[0, 1]`.
    ///
    /// Bernstein convex-hull subdivision encloses every coordinate extremum,
    /// including repeated roots, to floating-point resolution. This does not
    /// depend on sampled sign changes or an absolute derivative-value threshold.
    pub fn bounding_box(&self) -> BoundingBox<P::Scalar, { P::DIM }>
    where
        P: PointIndex,
    {
        let axes = core::array::from_fn(|axis| {
            let original: [P::Scalar; N] = core::array::from_fn(|i| self.control_points[i][axis]);
            let control_min = original.iter().copied().fold(original[0], Float::min);
            let control_max = original.iter().copied().fold(original[0], Float::max);
            // Subdivision in a local coordinate system makes the refinement
            // criterion translation independent. Restore the origin only once,
            // when the final extrema are rounded to the output coordinates.
            let mut origin = original[0];
            let mut control = original.map(|value| value - origin);
            if control.iter().any(|value| !value.is_finite()) {
                // Opposite finite extremes can have an unrepresentable range.
                // Their absolute coordinates remain safe for convex subdivision.
                origin = <P::Scalar as NumCast>::from(0.0).unwrap();
                control = original;
            }
            let start = control[0];
            let end = control[N - 1];
            let mut min = start.min(end);
            let mut max = start.max(end);
            let lo = control.iter().copied().fold(control[0], Float::min);
            let hi = control.iter().copied().fold(control[0], Float::max);
            let half = <P::Scalar as NumCast>::from(0.5).unwrap();
            let tolerance = ((hi * half - lo * half)
                * P::Scalar::epsilon()
                * <P::Scalar as NumCast>::from(32_u8).unwrap())
            .max(lo.abs().max(hi.abs()) * P::Scalar::epsilon());
            Self::bound_axis(control, tolerance, 0, &mut min, &mut max);
            (
                (min + origin)
                    .min(original[0])
                    .min(original[N - 1])
                    .max(control_min),
                (max + origin)
                    .max(original[0])
                    .max(original[N - 1])
                    .min(control_max),
            )
        });
        BoundingBox::from_axes_unchecked(axes)
    }

    fn bound_axis(
        control: [P::Scalar; N],
        tolerance: P::Scalar,
        depth: usize,
        min: &mut P::Scalar,
        max: &mut P::Scalar,
    ) {
        let lo = control.iter().copied().fold(control[0], Float::min);
        let hi = control.iter().copied().fold(control[0], Float::max);
        if lo >= *min && hi <= *max {
            return;
        }
        if depth >= MAX_ROOT_DEPTH || (*min - lo <= tolerance && hi - *max <= tolerance) {
            *min = min.min(lo);
            *max = max.max(hi);
            return;
        }
        if control.windows(2).all(|pair| pair[0] <= pair[1])
            || control.windows(2).all(|pair| pair[0] >= pair[1])
        {
            *min = min.min(control[0].min(control[N - 1]));
            *max = max.max(control[0].max(control[N - 1]));
            return;
        }
        let (left, right) = Self::subdivide_1d(control);
        Self::bound_axis(left, tolerance, depth + 1, min, max);
        Self::bound_axis(right, tolerance, depth + 1, min, max);
    }

    /// Return a conservative box around all control points.
    pub fn control_box(&self) -> BoundingBox<P::Scalar, { P::DIM }>
    where
        P: PointIndex,
    {
        BoundingBox::from_points(&self.control_points)
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

    fn push_root<const M: usize>(
        root: P::Scalar,
        tolerance: P::Scalar,
        results: &mut ArrayVec<[P::Scalar; M]>,
    ) where
        [P::Scalar; M]: tinyvec::Array<Item = P::Scalar>,
    {
        if results
            .iter()
            .any(|previous| (root - *previous).abs() <= tolerance)
        {
            return;
        }
        if results.len() < results.capacity() {
            let position = results
                .iter()
                .position(|previous| root < *previous)
                .unwrap_or(results.len());
            results.insert(position, root);
        }
    }

    fn subdivide_1d<const M: usize>(control: [P::Scalar; M]) -> ([P::Scalar; M], [P::Scalar; M]) {
        let half = <P::Scalar as NumCast>::from(0.5).unwrap();
        let mut left = core::array::from_fn(|_| <P::Scalar as NumCast>::from(0.0).unwrap());
        let mut right = core::array::from_fn(|_| <P::Scalar as NumCast>::from(0.0).unwrap());
        let mut temp = control;

        for i in 0..M {
            left[i] = temp[0];
            right[M - 1 - i] = temp[M - 1 - i];
            for j in 0..(M - 1 - i) {
                temp[j] = temp[j] * half + temp[j + 1] * half;
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
    use crate::CurveConstructionError;
    use crate::PointN;
    use crate::test_support::{assert_approx_equal, assert_point_approx_equal};

    const ABSOLUTE_TOLERANCE: f64 = 1e-12;
    const RELATIVE_TOLERANCE: f64 = 1e-12;

    #[test]
    fn rejects_invalid_control_points() {
        let empty = Bezier::<PointN<f64, 1>, 0>::new([]);
        assert!(matches!(
            empty,
            Err(CurveConstructionError::NoControlPoints)
        ));

        let non_finite = Bezier::new([
            PointN::new([0.0_f64]),
            PointN::new([f64::NAN]),
            PointN::new([f64::INFINITY]),
        ]);
        assert!(matches!(
            non_finite,
            Err(CurveConstructionError::NonFiniteControlPoint { index: 1 })
        ));
    }

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

        let curve: Bezier<PointN<f64, 2>, 6> = Bezier::new(points).unwrap();

        // check if start/end points match
        let start = curve.eval(0.0);
        assert_point_approx_equal(start, points[0], ABSOLUTE_TOLERANCE, RELATIVE_TOLERANCE);

        let end = curve.eval(1.0);
        assert_point_approx_equal(
            end,
            points[points.len() - 1],
            ABSOLUTE_TOLERANCE,
            RELATIVE_TOLERANCE,
        );
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
        let bezier: Bezier<PointN<f64, 2>, 2> =
            Bezier::new([PointN::new([0.0, 0.0]), PointN::new([10.0, 0.0])]).unwrap();
        let expected = (bezier.eval(1.0) - bezier.eval(0.0)).squared_norm().sqrt();
        let length = bezier.arclen_polyline(32);
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
            assert_point_approx_equal(
                bezier.eval(t / 2.0),
                left.eval(t),
                ABSOLUTE_TOLERANCE,
                RELATIVE_TOLERANCE,
            );
            // check the right part of the split curve
            assert_point_approx_equal(
                bezier.eval((t * 0.5) + 0.5),
                right.eval(t),
                ABSOLUTE_TOLERANCE,
                RELATIVE_TOLERANCE,
            );
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
        )
        .unwrap();

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
            assert_point_approx_equal(
                cubic_bezier.eval(t),
                generic_bezier.eval(t),
                ABSOLUTE_TOLERANCE,
                RELATIVE_TOLERANCE,
            );
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
        )
        .unwrap();

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
            assert_point_approx_equal(
                quadratic_bezier.eval(t),
                generic_bezier.eval(t),
                ABSOLUTE_TOLERANCE,
                RELATIVE_TOLERANCE,
            );
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
        )
        .unwrap();

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
        ])
        .unwrap();

        assert_eq!(bezier.start(), bezier.eval(0.0));
        assert_eq!(bezier.end(), bezier.eval(1.0));

        let reversed = bezier.reverse();
        let p0 = bezier.eval(0.25);
        let p1 = reversed.eval(0.75);
        assert_point_approx_equal(p0, p1, ABSOLUTE_TOLERANCE, RELATIVE_TOLERANCE);

        let tangent = bezier.tangent(0.3);
        assert_approx_equal(tangent[0], 1.0, ABSOLUTE_TOLERANCE, RELATIVE_TOLERANCE);
        assert_approx_equal(tangent[1], 0.0, ABSOLUTE_TOLERANCE, RELATIVE_TOLERANCE);

        let curvature = bezier.curvature(0.3);
        assert_approx_equal(curvature, 0.0, ABSOLUTE_TOLERANCE, RELATIVE_TOLERANCE);
        assert!(bezier.normal(0.3).is_none());

        let t = bezier.t_at_length(1.0, 1e-12).unwrap();
        assert!((t - 0.5).abs() < 1e-3);

        let p = bezier.point_at_length(1.0, 1e-12).unwrap();
        assert_point_approx_equal(p, PointN::new([1.0, 0.0]), 1e-9, 1e-9);
    }

    #[test]
    fn small_scale_bezier_remains_non_degenerate() {
        let scale = 1e-20_f64;
        let bezier =
            Bezier::<PointN<f64, 2>, 2>::new([PointN::new([0.0, 0.0]), PointN::new([scale, 0.0])])
                .unwrap();

        assert_point_approx_equal(
            bezier.tangent(0.5),
            PointN::new([1.0, 0.0]),
            ABSOLUTE_TOLERANCE,
            RELATIVE_TOLERANCE,
        );
        assert_approx_equal(
            bezier.t_at_length(scale * 0.5, 1e-12).unwrap(),
            0.5,
            1e-9,
            1e-9,
        );
    }
}
