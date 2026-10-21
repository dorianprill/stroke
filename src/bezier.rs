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
#[derive(Clone, Copy, Debug, PartialEq)]
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

impl<'a, P: Point, const N: usize> IntoIterator for &'a Bezier<P, { N }> {
    type Item = &'a P;
    type IntoIter = slice::Iter<'a, P>;

    fn into_iter(self) -> Self::IntoIter {
        self.control_points.iter()
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

    /// Borrow the control points in curve order.
    pub fn control_points(&self) -> slice::Iter<'_, P> {
        self.control_points.iter()
    }

    /// Consume the curve, returning its control-point array.
    pub fn into_control_points(self) -> [P; N] {
        self.control_points
    }

    /// Replace all control points atomically, leaving the curve unchanged on error.
    /// Mutable control-point iteration is deliberately unavailable:
    /// ```compile_fail
    /// use stroke::{Bezier, ArrayPoint};
    /// let mut curve = Bezier::new([ArrayPoint::new([0.0_f64]), ArrayPoint::new([1.0])]).unwrap();
    /// for point in &mut curve { *point = ArrayPoint::new([f64::NAN]); }
    /// ```
    pub fn set_control_points(&mut self, points: [P; N]) -> Result<(), CurveConstructionError> {
        *self = Self::new(points)?;
        Ok(())
    }

    /// Map control points and validate the result without modifying this curve.
    /// This maps the entire curve geometrically only when `map` is affine.
    pub fn try_map_control_points(
        &self,
        map: impl FnMut(P) -> P,
    ) -> Result<Self, CurveConstructionError> {
        Self::new(self.control_points.map(map))
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

    /// Return the `N` Bernstein basis weights in control-point order.
    ///
    /// For degree `n = N - 1`, weight `i` is `binomial(n, i) * t^i * (1-t)^(n-i)`.
    /// The recurrence uses no factorials and allocates no heap storage.
    /// Finite parameters outside `[0, 1]` extrapolate the polynomial; weights
    /// there need not be positive. In `[0, 1]` they are nonnegative and sum to
    /// one, up to roundoff. Endpoint weights are exact.
    ///
    /// Underflow may round small weights to zero. Extrapolation can magnify
    /// rounding errors, especially at high degrees. To evaluate a curve point,
    /// prefer [`Self::eval`]; summing weighted controls can overflow or lose
    /// precision even when de Casteljau evaluation remains usable.
    ///
    /// # Errors
    ///
    /// Return [`crate::BezierBasisError::NonFiniteParameter`] for NaN or infinity,
    /// or [`crate::BezierBasisError::NonFiniteComputation`] if a weight cannot
    /// be represented as a finite scalar.
    ///
    /// ```
    /// use stroke::{ArrayPoint, Bezier};
    /// let curve = Bezier::new([0.0_f64, 2.0, 4.0].map(|x| ArrayPoint::from([x])))?;
    /// assert_eq!(curve.basis_functions(0.5)?, [0.25, 0.5, 0.25]);
    /// assert_eq!(curve.basis_functions(2.0)?, [1.0, -4.0, 4.0]);
    /// # Ok::<(), Box<dyn core::error::Error>>(())
    /// ```
    pub fn basis_functions(&self, t: P::Scalar) -> Result<[P::Scalar; N], crate::BezierBasisError> {
        Ok(crate::bernstein::basis::<_, N, 0>(t)?[0])
    }

    /// Return Bernstein basis values and parameter derivatives through order `R`.
    ///
    /// Entry `[r][i]` is the `r`th derivative of weight `i` with respect to `t`.
    /// Row zero matches [`Self::basis_functions`]. Orders above the degree are
    /// exactly zero, including all positive orders for a constant curve.
    /// Parameter derivatives use the same `[0, 1]` units as [`Self::derivative_at`].
    ///
    /// Finite extrapolation, roundoff and errors follow [`Self::basis_functions`].
    /// Any unrepresentable requested derivative returns
    /// [`crate::BezierBasisError::NonFiniteComputation`]. Workspace uses
    /// `O(N * (R + 1))` stack storage, with no heap allocation.
    ///
    /// ```
    /// #![allow(incomplete_features)]
    /// #![feature(generic_const_exprs)]
    /// use stroke::{ArrayPoint, Bezier};
    /// let curve = Bezier::new([0.0_f64, 2.0, 4.0].map(|x| ArrayPoint::from([x])))?;
    /// let basis = curve.basis_functions_with_derivatives::<3>(0.5)?;
    /// assert_eq!(basis[1], [-1.0, 0.0, 1.0]);
    /// assert_eq!(basis[2], [2.0, -4.0, 2.0]);
    /// assert_eq!(basis[3], [0.0; 3]);
    /// # Ok::<(), Box<dyn core::error::Error>>(())
    /// ```
    pub fn basis_functions_with_derivatives<const R: usize>(
        &self,
        t: P::Scalar,
    ) -> Result<[[P::Scalar; N]; R + 1], crate::BezierBasisError>
    where
        [(); R + 1]: Sized,
    {
        crate::bernstein::basis::<_, N, R>(t)
    }

    /// Evaluate a point on the curve at point 't' which should be in the interval \[0, 1\] (unchecked!)
    /// This is implemented using De Casteljau's algorithm (over a temporary array with const generic sizing)
    pub fn eval(&self, t: P::Scalar) -> P {
        // start with a copy of the original control points array and succesively use it for evaluation
        let mut p: [P; N] = self.control_points;
        // loop up to degree = control_points.len() -1
        for i in 1..p.len() {
            for j in 0..p.len() - i {
                p[j] = crate::numeric::lerp(p[j], p[j + 1], t);
            }
        }
        p[0]
    }

    /// Approximate the minimum distance between given `point` and the curve.
    /// Uses a coarse sampling pass over the full domain and a local search around
    /// the best sample. `nsteps` counts coarse subdivisions (at least one).
    /// Every sampled/refined candidate is retained if better; the search has no
    /// global-minimum or accuracy guarantee. Non-finite query points return NaN;
    /// distances outside the scalar range are infinite.
    pub fn distance_to_point_approx(&self, point: P, nsteps: usize) -> P::Scalar
    where
        P: PointNorm,
    {
        crate::distance::approximate(point, nsteps, |t| {
            Ok::<_, core::convert::Infallible>(self.eval(t))
        })
        .unwrap_or_else(|never| match never {})
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
                    crate::numeric::lerp(casteljau_points[j], casteljau_points[j + 1], t);
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

    /// Extract `start..end`, reparameterized to `[0, 1]`.
    /// Requires finite `0 <= start <= end <= 1`; equal endpoints yield a constant curve.
    pub fn subcurve(
        &self,
        range: core::ops::Range<P::Scalar>,
    ) -> Result<Self, crate::CurveParameterError> {
        let zero = <P::Scalar as NumCast>::from(0.0).unwrap();
        let one = <P::Scalar as NumCast>::from(1.0).unwrap();
        let (start, end) = (range.start, range.end);
        if !start.is_finite() || !end.is_finite() || start < zero || end > one || start > end {
            return Err(crate::CurveParameterError::InvalidRange);
        }
        let result = if start == zero && end == one {
            *self
        } else if start == end {
            Self {
                control_points: [self.eval(start); N],
            }
        } else if start == zero {
            self.split(end).0
        } else {
            self.split(start).1.split((end - start) / (one - start)).0
        };
        if result.control_points.iter().all(Point::is_finite) {
            Ok(result)
        } else {
            Err(crate::CurveParameterError::NonFiniteComputation)
        }
    }

    /// Split at a distance from the start, clamping finite distances to the endpoints.
    pub fn split_at_length(
        &self,
        length: P::Scalar,
        accuracy: P::Scalar,
    ) -> Result<(Self, Self), ArcLengthError>
    where
        P: PointNorm,
    {
        Ok(self.split(self.t_at_length(length, accuracy)?))
    }

    /// Evaluate the derivative with respect to `t` without constructing another curve.
    /// Constant curves return zero; `t` is expected to be in `[0, 1]`.
    pub fn derivative_at(&self, t: P::Scalar) -> P {
        self.velocity(t)
    }

    /// Return the derivative magnitude with respect to `t` in `[0, 1]`.
    pub fn speed(&self, t: P::Scalar) -> P::Scalar
    where
        P: PointNorm,
    {
        self.derivative_at(t).norm()
    }

    /// Returns the derivative curve of self which has N-1 control points.
    /// The derivative of an nth degree Bezier curve is an (n-1)th degree Bezier curve,
    /// with one fewer term, and new weights w0...wn-1 derived from the
    /// original weights as n(wi+1 - wi). So for a 3rd degree curve, with four weights,
    /// the derivative has three new weights:
    ///     w0 = 3(w1-w0), w'1 = 3(w2-w1) and w'2 = 3(w3-w2).
    /// Requires at least two control points. For constant curves use
    /// [`Self::derivative_at`], which returns zero.
    /// ```compile_fail
    /// #![feature(generic_const_exprs)]
    /// #![allow(incomplete_features)]
    /// use stroke::{ArrayPoint, Bezier};
    /// let constant = Bezier::new([ArrayPoint::from([1.0_f64])]).unwrap();
    /// let _ = constant.derivative();
    /// ```
    pub fn derivative(&self) -> Bezier<P, { N - 1 }>
    where
        [(); N - 2]: Sized,
    {
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
    {
        crate::differential::Differential::bezier(self.control_points, t, false).tangent()
    }

    /// Return the curvature magnitude at `t`.
    pub fn curvature(&self, t: P::Scalar) -> P::Scalar
    where
        P: PointNorm,
    {
        crate::differential::Differential::bezier(self.control_points, t, true).curvature()
    }

    /// Return the principal normal direction at `t`.
    ///
    /// Returns `None` if the velocity is zero or curvature is undefined.
    pub fn normal(&self, t: P::Scalar) -> Option<P>
    where
        P: PointNorm,
    {
        crate::differential::Differential::bezier(self.control_points, t, true).normal()
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
        crate::differential::Differential::bezier(self.control_points, t, false).derivative()
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
    /// `eps` is a positive finite absolute coordinate residual tolerance (default
    /// `1e-6`). This local search can fail even when a root exists.
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
{
    fn parameter_domain(&self) -> (P::Scalar, P::Scalar) {
        (
            <P::Scalar as NumCast>::from(0.0).unwrap(),
            <P::Scalar as NumCast>::from(1.0).unwrap(),
        )
    }

    fn axis_value(&self, t: P::Scalar, axis: usize) -> Result<P::Scalar, RootFindingError> {
        if axis >= P::DIM {
            return Err(RootFindingError::InvalidAxis);
        }
        if !t.is_finite()
            || t < <P::Scalar as NumCast>::from(0.0).unwrap()
            || t > <P::Scalar as NumCast>::from(1.0).unwrap()
        {
            return Err(RootFindingError::InvalidParameter);
        }
        Ok(self.eval(t)[axis])
    }

    fn axis_derivative(&self, t: P::Scalar, axis: usize) -> Result<P::Scalar, RootFindingError> {
        if axis >= P::DIM {
            return Err(RootFindingError::InvalidAxis);
        }
        if !t.is_finite()
            || t < <P::Scalar as NumCast>::from(0.0).unwrap()
            || t > <P::Scalar as NumCast>::from(1.0).unwrap()
        {
            return Err(RootFindingError::InvalidParameter);
        }
        Ok(self.derivative_at(t)[axis])
    }
}

#[cfg(test)]
mod tests {
    use super::CubicBezier;
    use super::QuadraticBezier;
    use super::*;
    use crate::ArrayPoint;
    use crate::CurveConstructionError;
    use crate::test_support::{assert_approx_equal, assert_point_approx_equal};

    const ABSOLUTE_TOLERANCE: f64 = 1e-12;
    const RELATIVE_TOLERANCE: f64 = 1e-12;

    #[test]
    fn rejects_invalid_control_points() {
        let empty = Bezier::<ArrayPoint<f64, 1>, 0>::new([]);
        assert!(matches!(
            empty,
            Err(CurveConstructionError::NoControlPoints)
        ));

        let non_finite = Bezier::new([
            ArrayPoint::new([0.0_f64]),
            ArrayPoint::new([f64::NAN]),
            ArrayPoint::new([f64::INFINITY]),
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
            ArrayPoint::new([0f64, 1.77f64]),
            ArrayPoint::new([1.1f64, -1f64]),
            ArrayPoint::new([4.3f64, 3f64]),
            ArrayPoint::new([3.2f64, -4f64]),
            ArrayPoint::new([7.3f64, 2.7f64]),
            ArrayPoint::new([8.9f64, 1.7f64]),
        ];

        let curve: Bezier<ArrayPoint<f64, 2>, 6> = Bezier::new(points).unwrap();

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
                ArrayPoint::new([0f64, 1.77f64]),
                ArrayPoint::new([2.9f64, 0f64]),
                ArrayPoint::new([4.3f64, 3f64]),
                ArrayPoint::new([3.2f64, -4f64]),
            ],
        };
        assert!(
            bezier.distance_to_point(ArrayPoint::new([-5.1, -5.6]))
                > bezier.distance_to_point(ArrayPoint::new([5.1, 5.6]))
        );
    }

    #[test]
    fn arclen_line_approx() {
        let bezier: Bezier<ArrayPoint<f64, 2>, 2> =
            Bezier::new([ArrayPoint::new([0.0, 0.0]), ArrayPoint::new([10.0, 0.0])]).unwrap();
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
                ArrayPoint::new([0f64, 1.77f64]),
                ArrayPoint::new([2.9f64, 0f64]),
                ArrayPoint::new([4.3f64, 3f64]),
                ArrayPoint::new([3.2f64, -4f64]),
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
            ArrayPoint::new([0f64, 1.77f64]),
            ArrayPoint::new([1.1f64, -1f64]),
            ArrayPoint::new([4.3f64, 3f64]),
            ArrayPoint::new([3.2f64, -4f64]),
        )
        .unwrap();

        let generic_bezier = Bezier {
            control_points: [
                ArrayPoint::new([0f64, 1.77f64]),
                ArrayPoint::new([1.1f64, -1f64]),
                ArrayPoint::new([4.3f64, 3f64]),
                ArrayPoint::new([3.2f64, -4f64]),
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
            ArrayPoint::new([0f64, 1.77f64]),
            ArrayPoint::new([1.1f64, -1f64]),
            ArrayPoint::new([3.2f64, -4f64]),
        )
        .unwrap();

        let generic_bezier = Bezier {
            control_points: [
                ArrayPoint::new([0f64, 1.77f64]),
                ArrayPoint::new([1.1f64, -1f64]),
                ArrayPoint::new([3.2f64, -4f64]),
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
                ArrayPoint::new([0f64]),
                ArrayPoint::new([1f64]),
                ArrayPoint::new([0f64]),
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
                ArrayPoint::new([0f64]),
                ArrayPoint::new([4f64]),
                ArrayPoint::new([-4f64]),
                ArrayPoint::new([0f64]),
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
                ArrayPoint::new([0f64]),
                ArrayPoint::new([1f64]),
                ArrayPoint::new([2f64]),
            ],
        };

        let roots = bezier.derivative_roots(0);
        assert!(roots.is_empty());
    }

    #[test]
    fn bounding_box_cubic_1d_extrema() {
        let bezier = Bezier {
            control_points: [
                ArrayPoint::new([0f64]),
                ArrayPoint::new([4f64]),
                ArrayPoint::new([-4f64]),
                ArrayPoint::new([0f64]),
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
            control_points: [ArrayPoint::new([0f64]), ArrayPoint::new([2f64])],
        };

        let root = bezier
            .root_newton_axis(1.0, 0, 0.25, None, Some(64))
            .unwrap();
        assert!((root - 0.5).abs() < 1e-6);
    }

    #[test]
    fn root_newton_axis_zero_derivative() {
        let bezier = Bezier {
            control_points: [ArrayPoint::new([1f64]), ArrayPoint::new([1f64])],
        };

        let result = bezier.root_newton_axis(0.0, 0, 0.5, None, Some(8));
        assert!(matches!(result, Err(RootFindingError::ZeroDerivative)));
    }

    #[test]
    fn bounding_box_matches_cubic_specialization() {
        let generic_bezier = Bezier {
            control_points: [
                ArrayPoint::new([0f64, 1.77f64]),
                ArrayPoint::new([1.1f64, -1f64]),
                ArrayPoint::new([4.3f64, 3f64]),
                ArrayPoint::new([3.2f64, -4f64]),
            ],
        };
        let cubic_bezier = CubicBezier::new(
            ArrayPoint::new([0f64, 1.77f64]),
            ArrayPoint::new([1.1f64, -1f64]),
            ArrayPoint::new([4.3f64, 3f64]),
            ArrayPoint::new([3.2f64, -4f64]),
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
        let bezier = Bezier::<ArrayPoint<f64, 2>, 2>::new([
            ArrayPoint::new([0f64, 0f64]),
            ArrayPoint::new([2f64, 0f64]),
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
        assert_point_approx_equal(p, ArrayPoint::new([1.0, 0.0]), 1e-9, 1e-9);
    }

    #[test]
    fn small_scale_bezier_remains_non_degenerate() {
        let scale = 1e-20_f64;
        let bezier = Bezier::<ArrayPoint<f64, 2>, 2>::new([
            ArrayPoint::new([0.0, 0.0]),
            ArrayPoint::new([scale, 0.0]),
        ])
        .unwrap();

        assert_point_approx_equal(
            bezier.tangent(0.5),
            ArrayPoint::new([1.0, 0.0]),
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
