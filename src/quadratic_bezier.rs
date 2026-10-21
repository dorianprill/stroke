//! Quadratic Bezier curve specialization.

use super::{ArcLengthError, ArrayVec, BoundingBox, LineSegment, Point, PointIndex, PointNorm};
use crate::error::{CurveConstructionError, validate_control_points};
use crate::numeric::is_negligible;
use num_traits::{Float, NumCast};

const DEFAULT_DISTANCE_STEPS: usize = 64;

/// Quadratic Bezier curve defined by start, control, and end points.
///
/// Bounding boxes require `PointIndex`. Distance and differential geometry
/// methods require `PointNorm`.
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
    /// Create a quadratic Bezier curve from its three control points.
    ///
    /// # Errors
    ///
    /// Returns [`CurveConstructionError::NonFiniteControlPoint`] with the index
    /// of the first control point containing a non-finite coordinate.
    pub fn new(start: P, ctrl: P, end: P) -> Result<Self, CurveConstructionError> {
        validate_control_points(&[start, ctrl, end])?;
        Ok(QuadraticBezier { start, ctrl, end })
    }

    /// Evaluate at `t` using endpoint-preserving de Casteljau interpolation.
    pub fn eval(&self, t: P::Scalar) -> P {
        crate::Bezier::from(*self).eval(t)
    }

    /// Evaluate at `t` using de Casteljau interpolation, as in [`Self::eval`].
    pub fn eval_casteljau(&self, t: P::Scalar) -> P {
        self.eval(t)
    }

    /// Return the 3 Bernstein basis weights in control-point order.
    ///
    /// Finite parameters outside `[0, 1]` extrapolate the polynomial.
    /// See [`Bezier::basis_functions`](crate::Bezier::basis_functions) for
    /// numerical limits and [`BezierBasisError`](crate::BezierBasisError).
    pub fn basis_functions(&self, t: P::Scalar) -> Result<[P::Scalar; 3], crate::BezierBasisError> {
        Ok(crate::bernstein::basis::<_, 3, 0>(t)?[0])
    }

    /// Return basis values and parameter derivatives through order `R`.
    ///
    /// Entry `[r][i]` corresponds to derivative order `r` and control `i`.
    /// Orders above 2 are zero. Finite parameters may extrapolate.
    /// See [`Bezier::basis_functions_with_derivatives`](crate::Bezier::basis_functions_with_derivatives)
    /// for numerical limits, errors and workspace requirements.
    pub fn basis_functions_with_derivatives<const R: usize>(
        &self,
        t: P::Scalar,
    ) -> Result<[[P::Scalar; 3]; R + 1], crate::BezierBasisError>
    where
        [(); R + 1]: Sized,
    {
        crate::bernstein::basis::<_, 3, R>(t)
    }

    /// Borrow the control points in curve order.
    pub fn control_points(&self) -> core::array::IntoIter<&P, 3> {
        [&self.start, &self.ctrl, &self.end].into_iter()
    }

    /// Return the control-point array, consuming the curve.
    pub fn into_control_points(self) -> [P; 3] {
        [self.start, self.ctrl, self.end]
    }

    /// Replace all control points. Leave the curve unchanged if validation fails.
    pub fn set_control_points(&mut self, points: [P; 3]) -> Result<(), CurveConstructionError> {
        *self = Self::new(points[0], points[1], points[2])?;
        Ok(())
    }

    /// Return a curve with mapped control points. Reject non-finite results.
    /// Affine maps, such as translations, also map every point of the curve.
    pub fn try_map_control_points(
        &self,
        map: impl FnMut(P) -> P,
    ) -> Result<Self, CurveConstructionError> {
        let points = (*self).into_control_points().map(map);
        Self::new(points[0], points[1], points[2])
    }

    /// Evaluate the derivative with respect to `t` in `[0, 1]`.
    pub fn derivative_at(&self, t: P::Scalar) -> P {
        crate::Bezier::from(*self).derivative_at(t)
    }

    /// Return the speed, the magnitude of the derivative with respect to `t`.
    pub fn speed(&self, t: P::Scalar) -> P::Scalar
    where
        P: PointNorm,
    {
        self.derivative_at(t).norm()
    }

    /// Extract `start..end` and give the result a new parameter range `[0, 1]`.
    /// Requires finite `0 <= start <= end <= 1`. Equal endpoints give a constant curve.
    pub fn subcurve(
        &self,
        range: core::ops::Range<P::Scalar>,
    ) -> Result<Self, crate::CurveParameterError> {
        crate::Bezier::from(*self).subcurve(range).map(Self::from)
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
        let (left, right) = crate::Bezier::from(*self).split_at_length(length, accuracy)?;
        Ok((Self::from(left), Self::from(right)))
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
        let (left, right) = crate::Bezier::from(*self).split(t);
        (Self::from(left), Self::from(right))
    }

    /// Return the derivative as a line segment.
    /// Evaluate it with `derivative().eval(t)` or use [`Self::derivative_at`].
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
        crate::Bezier::from(*self).tangent(t)
    }

    /// Return the curvature magnitude at `t`.
    pub fn curvature(&self, t: P::Scalar) -> P::Scalar
    where
        P: PointNorm,
    {
        crate::Bezier::from(*self).curvature(t)
    }

    /// Return the principal normal direction at `t`.
    ///
    /// Returns `None` if the velocity is zero or curvature is undefined.
    pub fn normal(&self, t: P::Scalar) -> Option<P>
    where
        P: PointNorm,
    {
        crate::Bezier::from(*self).normal(t)
    }

    /// Approximate arc length with `segments` straight lines at equal parameter steps.
    pub fn arclen_polyline(&self, segments: usize) -> P::Scalar
    where
        P: PointNorm,
    {
        crate::Bezier {
            control_points: [self.start, self.ctrl, self.end],
        }
        .arclen_polyline(segments)
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
    /// `accuracy` must be finite and positive, in the same units as the coordinates.
    pub fn arclen_between(
        &self,
        start: P::Scalar,
        end: P::Scalar,
        accuracy: P::Scalar,
    ) -> Result<P::Scalar, ArcLengthError>
    where
        P: PointNorm,
    {
        crate::Bezier {
            control_points: [self.start, self.ctrl, self.end],
        }
        .arclen_between(start, end, accuracy)
    }

    /// Return the sorted real roots of `a*t^2 + b*t + c`.
    /// Call once for each coordinate polynomial.
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
        let zero = <P::Scalar as NumCast>::from(0.0).unwrap();
        let two = <P::Scalar as NumCast>::from(2.0).unwrap();
        let four = <P::Scalar as NumCast>::from(4.0).unwrap();
        let coefficient_scale = a.abs().max(b.abs()).max(c.abs());
        if coefficient_scale == zero || !coefficient_scale.is_finite() {
            return result;
        }
        let a = a / coefficient_scale;
        let b = b / coefficient_scale;
        let c = c / coefficient_scale;
        let coefficient_scale = <P::Scalar as NumCast>::from(1.0).unwrap();

        // Treat negligible leading coefficients as lower-degree polynomials.
        if is_negligible(a, coefficient_scale) {
            if is_negligible(b, coefficient_scale) {
                // The polynomial is constant and has no isolated roots.
                return result;
            }
            result.push(-c / b);
            return result;
        }
        let b_squared = b * b;
        let four_ac = a * c * four;
        let delta = b_squared - four_ac;
        if is_negligible(delta, b_squared.abs() + four_ac.abs()) {
            result.push(-b / (a * two));
        } else if delta > zero {
            let sqrt_delta = delta.sqrt();
            // The q form avoids cancellation in one of the two roots.
            let q = -(b + b.signum() * sqrt_delta) / two;
            result.push(q / a);
            result.push(c / q);
        }
        result.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap());
        result
    }

    /// Estimate the minimum distance to `point`.
    /// Sample at `nsteps` equal parameter intervals, then refine near the closest
    /// sample. Use at least one interval and retain the best distance found.
    /// This local search can miss a closer point and has no error bound.
    /// Non-finite query points return NaN. Distances outside the scalar range are infinite.
    pub fn distance_to_point_approx(&self, point: P, nsteps: usize) -> P::Scalar
    where
        P: PointNorm,
    {
        crate::Bezier::from(*self).distance_to_point_approx(point, nsteps)
    }

    /// Estimate the minimum distance using the default number of parameter intervals.
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
        Ok(self.eval_casteljau(self.t_at_length(s, accuracy)?))
    }

    /// Return the line joining the curve's endpoints.
    pub fn baseline(&self) -> LineSegment<P> {
        LineSegment {
            start: self.start,
            end: self.end,
        }
    }

    /// Return whether the controls lie within `abs(tolerance)` of the endpoint segment.
    /// The endpoints must be farther apart than the tolerance.
    /// Non-finite tolerances return `false`. Negative tolerances use their absolute value.
    pub fn is_linear(&self, tolerance: P::Scalar) -> bool
    where
        P: PointNorm,
    {
        let tolerance = tolerance.abs();
        tolerance.is_finite()
            && crate::distance::between(self.start, self.end) > tolerance
            && self.are_points_colinear(tolerance)
    }

    /// Check the control point's distance to the endpoint segment.
    fn are_points_colinear(&self, tolerance: P::Scalar) -> bool
    where
        P: PointNorm,
    {
        let line = self.baseline();
        line.distance_to_point(self.ctrl) <= tolerance
    }

    /// Return whether every control point is within `abs(tolerance)` of the start.
    /// Zero requires equal coordinates. Non-finite tolerances return `false`.
    pub fn is_a_point(&self, tolerance: P::Scalar) -> bool
    where
        P: PointNorm,
    {
        crate::distance::is_a_point(&[self.start, self.ctrl, self.end], tolerance)
    }

    /// Return parameters in `[0, 1]` where coordinate `axis` equals `value`.
    /// Solve the coordinate polynomial analytically.
    #[allow(dead_code)]
    fn solve_t_for_axis(&self, value: P::Scalar, axis: usize) -> ArrayVec<[P::Scalar; 3]>
    where
        P: PointIndex + PointNorm,
        [P::Scalar; 2]: tinyvec::Array<Item = P::Scalar>,
        [P::Scalar; 3]: tinyvec::Array<Item = P::Scalar>,
    {
        let mut result = ArrayVec::new();
        let zero = <P::Scalar as NumCast>::from(0.0).unwrap();
        if self.is_a_point(zero)
            || (self.are_points_colinear(zero) && (self.start - self.end).squared_norm() == zero)
        {
            return result;
        }
        let a = self.start[axis]
            + self.ctrl[axis] * <P::Scalar as NumCast>::from(-2.0).unwrap()
            + self.end[axis];
        let b = self.start[axis] * <P::Scalar as NumCast>::from(-2.0).unwrap()
            + self.ctrl[axis] * <P::Scalar as NumCast>::from(2.0).unwrap();
        let c = self.start[axis] - value;

        let roots = self.real_roots(a, b, c);
        let one = <P::Scalar as NumCast>::from(1.0).unwrap();
        for &root in roots.iter() {
            if root > zero && root < one {
                result.push(root);
            }
        }

        result
    }

    /// Return the tight axis-aligned box using analytic coordinate extrema.
    pub fn bounding_box(&self) -> BoundingBox<P::Scalar, { P::DIM }>
    where
        P: PointIndex,
        [P::Scalar; 1]: tinyvec::Array<Item = P::Scalar>,
        [P::Scalar; 3]: tinyvec::Array<Item = P::Scalar>,
    {
        let zero = <P::Scalar as NumCast>::from(0.0).unwrap();
        let one = <P::Scalar as NumCast>::from(1.0).unwrap();
        let mut bounds = [(zero, zero); P::DIM];
        let derivative = self.derivative();
        // Form the linear derivative a*t + b for each coordinate.
        let a = derivative.start * <P::Scalar as NumCast>::from(-1.0).unwrap() + derivative.end;
        let b = derivative.start;
        if !a.is_finite() || !b.is_finite() {
            return crate::Bezier {
                control_points: [self.start, self.ctrl, self.end],
            }
            .bounding_box();
        }

        for dim in 0..P::DIM {
            // Coordinate extrema occur at derivative roots or curve endpoints.
            let mut extrema: ArrayVec<[P::Scalar; 3]> = ArrayVec::new();
            extrema.extend(derivative.root(a[dim], b[dim]).iter().copied());
            extrema.retain(|root| -> bool { *root > zero && *root < one });
            for t in extrema.iter_mut() {
                *t = self.eval_casteljau(*t)[dim];
            }
            // Include both endpoints as bounds candidates.
            extrema.push(self.start[dim]);
            extrema.push(self.end[dim]);
            extrema.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap());

            // The candidate set contains both endpoints, so it is nonempty.
            bounds[dim] = (extrema[0], *extrema.last().unwrap());
        }
        BoundingBox::from_axes_unchecked(bounds)
    }

    /// Fit an oriented box around the curve on `[0, 1]`, testing up to 256 orientations.
    ///
    /// Delegates to the generic Bézier search without degree-specific optimizations.
    /// The heuristic search guarantees neither a local nor a global minimum.
    /// See [`Bezier::bounding_box_oriented_approx`](crate::Bezier::bounding_box_oriented_approx)
    /// for numerical limits and errors.
    pub fn bounding_box_oriented_approx(
        &self,
    ) -> Result<crate::OrientedBoundingBox<P::Scalar, { P::DIM }>, crate::OrientedBoundingBoxError>
    where
        P: PointIndex,
    {
        crate::Bezier::from(*self).bounding_box_oriented_approx()
    }

    /// Fit an oriented box with at most `max_evaluations` full-curve trials.
    ///
    /// Delegates to the generic Bézier search. See
    /// [`Bezier::bounding_box_oriented_approx_with_budget`](crate::Bezier::bounding_box_oriented_approx_with_budget)
    /// for work limits, stack usage and errors. The budget does not bound volume error.
    pub fn bounding_box_oriented_approx_with_budget(
        &self,
        max_evaluations: usize,
    ) -> Result<crate::OrientedBoundingBox<P::Scalar, { P::DIM }>, crate::OrientedBoundingBoxError>
    where
        P: PointIndex,
    {
        crate::Bezier::from(*self).bounding_box_oriented_approx_with_budget(max_evaluations)
    }

    /// Return a conservative box around the control points.
    pub fn control_box(&self) -> BoundingBox<P::Scalar, { P::DIM }>
    where
        P: PointIndex,
    {
        BoundingBox::from_points(&[self.start, self.ctrl, self.end])
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ArrayPoint;
    use crate::CurveConstructionError;
    use crate::test_support::{assert_approx_equal, assert_point_approx_equal};

    const ABSOLUTE_TOLERANCE: f64 = 1e-12;
    const RELATIVE_TOLERANCE: f64 = 1e-12;

    #[test]
    fn rejects_non_finite_control_points() {
        let result = QuadraticBezier::new(
            ArrayPoint::new([0.0_f64]),
            ArrayPoint::new([f64::NAN]),
            ArrayPoint::new([f64::INFINITY]),
        );

        assert_eq!(
            result,
            Err(CurveConstructionError::NonFiniteControlPoint { index: 1 })
        );
    }

    #[test]
    fn eval_endpoints() {
        let bezier = QuadraticBezier::new(
            ArrayPoint::new([0f64, 0f64]),
            ArrayPoint::new([1f64, 1f64]),
            ArrayPoint::new([2f64, 0f64]),
        )
        .unwrap();

        assert_eq!(bezier.eval(0.0), ArrayPoint::new([0.0, 0.0]));
        assert_eq!(bezier.eval(1.0), ArrayPoint::new([2.0, 0.0]));
    }

    #[test]
    fn arclen_line_approx() {
        let bezier = QuadraticBezier::new(
            ArrayPoint::new([0f64, 0f64]),
            ArrayPoint::new([1f64, 0f64]),
            ArrayPoint::new([2f64, 0f64]),
        )
        .unwrap();
        let length = bezier.arclen_polyline(32);
        assert!((length - 2.0).abs() < 1e-3);
    }

    #[test]
    fn eval_equivalence() {
        // Both evaluation entry points use de Casteljau interpolation.
        let bezier = QuadraticBezier::new(
            ArrayPoint::new([0f64, 1.77f64]),
            ArrayPoint::new([4.3f64, 3f64]),
            ArrayPoint::new([3.2f64, -4f64]),
        )
        .unwrap();

        let nsteps: usize = 1000;
        for t in 0..=nsteps {
            let t = t as f64 * 1f64 / (nsteps as f64);
            let p1 = bezier.eval(t);
            let p2 = bezier.eval_casteljau(t);
            assert_point_approx_equal(p2, p1, ABSOLUTE_TOLERANCE, RELATIVE_TOLERANCE);
        }
    }

    #[test]
    fn eval_equivalence_f32() {
        let bezier = QuadraticBezier::new(
            ArrayPoint::new([0.0_f32, 1.77]),
            ArrayPoint::new([4.3, 3.0]),
            ArrayPoint::new([3.2, -4.0]),
        )
        .unwrap();

        let nsteps: usize = 1000;
        for step in 0..=nsteps {
            let t = step as f32 / nsteps as f32;
            assert_point_approx_equal(bezier.eval(t), bezier.eval_casteljau(t), 1e-5, 1e-5);
        }
    }

    #[test]
    fn split_equivalence() {
        let bezier = QuadraticBezier {
            start: ArrayPoint::new([0f64, 1.77f64]),
            ctrl: ArrayPoint::new([4.3f64, 3f64]),
            end: ArrayPoint::new([3.2f64, -4f64]),
        };
        let at = 0.5;
        let (left, right) = bezier.split(at);
        // Each subcurve must match the corresponding parameter range of the original.
        let nsteps: usize = 1000;
        for t in 0..=nsteps {
            let t = t as f64 * 1f64 / (nsteps as f64);
            assert_point_approx_equal(
                bezier.eval(t / 2.0),
                left.eval(t),
                ABSOLUTE_TOLERANCE,
                RELATIVE_TOLERANCE,
            );
            assert_point_approx_equal(
                bezier.eval((t * 0.5) + 0.5),
                right.eval(t),
                ABSOLUTE_TOLERANCE,
                RELATIVE_TOLERANCE,
            );
        }
    }

    #[test]
    fn bounding_box_contains() {
        // Allow for floating-point error when comparing samples with the box.
        let bezier = QuadraticBezier {
            start: ArrayPoint::new([0f64, 1.77f64]),
            ctrl: ArrayPoint::new([4.3f64, -3f64]),
            end: ArrayPoint::new([3.2f64, 4f64]),
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
        // A cubic with four controls requires eight knots.
        let curve = QuadraticBezier {
            start: ArrayPoint::new([0f64, 1.77f64]),
            ctrl: ArrayPoint::new([4.3f64, 3f64]),
            end: ArrayPoint::new([3.2f64, -4f64]),
        };
        assert!(
            curve.distance_to_point(ArrayPoint::new([-5.1, -5.6]))
                > curve.distance_to_point(ArrayPoint::new([5.1, 5.6]))
        );
    }

    #[test]
    fn quadratic_api_parity() {
        let bezier = QuadraticBezier::new(
            ArrayPoint::new([0f64, 0f64]),
            ArrayPoint::new([1f64, 0f64]),
            ArrayPoint::new([2f64, 0f64]),
        )
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
        assert!((t - 0.5).abs() < 1e-6);

        let p = bezier.point_at_length(1.0, 1e-12).unwrap();
        assert_point_approx_equal(p, ArrayPoint::new([1.0, 0.0]), 1e-9, 1e-9);
    }

    #[test]
    fn quadratic_curvature_nonzero() {
        let bezier = QuadraticBezier::new(
            ArrayPoint::new([0f64, 0f64]),
            ArrayPoint::new([1f64, 1f64]),
            ArrayPoint::new([2f64, 0f64]),
        )
        .unwrap();

        let curvature = bezier.curvature(0.5);
        assert!(curvature > 0.0);

        let normal = bezier.normal(0.5).unwrap();
        assert!((normal.squared_norm() - 1.0).abs() < 1e-6);
    }

    #[test]
    fn small_scale_quadratic_remains_linear_with_unit_tangent() {
        let scale = 1e-20_f64;
        let bezier = QuadraticBezier::new(
            ArrayPoint::new([0.0, 0.0]),
            ArrayPoint::new([scale * 0.5, 0.0]),
            ArrayPoint::new([scale, 0.0]),
        )
        .unwrap();

        assert!(bezier.is_linear(scale * 1e-6));
        assert_point_approx_equal(
            bezier.tangent(0.5),
            ArrayPoint::new([1.0, 0.0]),
            ABSOLUTE_TOLERANCE,
            RELATIVE_TOLERANCE,
        );
    }

    #[test]
    fn quadratic_roots_are_scale_invariant() {
        let bezier = QuadraticBezier::new(
            ArrayPoint::new([0.0_f64]),
            ArrayPoint::new([0.0]),
            ArrayPoint::new([0.0]),
        )
        .unwrap();
        let roots = bezier.real_roots(1e-30, -3e-30, 2e-30);

        assert_eq!(roots.len(), 2);
        assert_approx_equal(roots[0], 1.0, ABSOLUTE_TOLERANCE, RELATIVE_TOLERANCE);
        assert_approx_equal(roots[1], 2.0, ABSOLUTE_TOLERANCE, RELATIVE_TOLERANCE);
    }
}
