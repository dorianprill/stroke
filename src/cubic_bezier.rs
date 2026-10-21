//! Cubic Bezier curve specialization.

use super::{
    ArcLengthError, ArrayVec, BoundingBox, LineSegment, Point, PointIndex, PointNorm,
    QuadraticBezier,
};
use crate::error::{CurveConstructionError, validate_control_points};
use crate::numeric::is_negligible;
use num_traits::{Float, NumCast};

const DEFAULT_DISTANCE_STEPS: usize = 64;

/// A cubic Bezier curve with endpoint controls `start` and `end`
/// and interior controls `ctrl1` and `ctrl2`.
///
/// For `t` in `[0, 1]`, the curve is
/// ```text
/// P(t) = (1 - t)³ * start
///            + 3 * (1 - t)² * t * ctrl1
///            + 3 * t² * (1 - t) * ctrl2
///            + t³ * end
/// ```
///
/// Bounding boxes require `PointIndex`. Distance and differential geometry
/// methods require `PointNorm`.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct CubicBezier<P> {
    pub(crate) start: P,
    pub(crate) ctrl1: P,
    pub(crate) ctrl2: P,
    pub(crate) end: P,
}

impl<P> CubicBezier<P>
where
    P: Point,
{
    /// Creates a new cubic Bezier curve from four control points.
    ///
    /// # Errors
    ///
    /// Returns [`CurveConstructionError::NonFiniteControlPoint`] with the index
    /// of the first control point containing a non-finite coordinate.
    pub fn new(start: P, ctrl1: P, ctrl2: P, end: P) -> Result<Self, CurveConstructionError> {
        validate_control_points(&[start, ctrl1, ctrl2, end])?;
        Ok(CubicBezier {
            start,
            ctrl1,
            ctrl2,
            end,
        })
    }

    /// Evaluate at `t` using endpoint-preserving de Casteljau interpolation.
    pub fn eval(&self, t: P::Scalar) -> P {
        crate::Bezier::from(*self).eval(t)
    }

    /// Evaluate at `t` using de Casteljau interpolation, as in [`Self::eval`].
    pub fn eval_casteljau(&self, t: P::Scalar) -> P {
        self.eval(t)
    }

    /// Return the 4 Bernstein basis weights in control-point order.
    ///
    /// Finite parameters outside `[0, 1]` extrapolate the polynomial.
    /// See [`Bezier::basis_functions`](crate::Bezier::basis_functions) for
    /// numerical limits and [`BezierBasisError`](crate::BezierBasisError).
    pub fn basis_functions(&self, t: P::Scalar) -> Result<[P::Scalar; 4], crate::BezierBasisError> {
        Ok(crate::bernstein::basis::<_, 4, 0>(t)?[0])
    }

    /// Return basis values and parameter derivatives through order `R`.
    ///
    /// Entry `[r][i]` corresponds to derivative order `r` and control `i`.
    /// Orders above 3 are zero. Finite parameters may extrapolate.
    /// See [`Bezier::basis_functions_with_derivatives`](crate::Bezier::basis_functions_with_derivatives)
    /// for numerical limits, errors and workspace requirements.
    pub fn basis_functions_with_derivatives<const R: usize>(
        &self,
        t: P::Scalar,
    ) -> Result<[[P::Scalar; 4]; R + 1], crate::BezierBasisError>
    where
        [(); R + 1]: Sized,
    {
        crate::bernstein::basis::<_, 4, R>(t)
    }

    /// Borrow the control points in curve order.
    pub fn control_points(&self) -> core::array::IntoIter<&P, 4> {
        [&self.start, &self.ctrl1, &self.ctrl2, &self.end].into_iter()
    }

    /// Return the control-point array, consuming the curve.
    pub fn into_control_points(self) -> [P; 4] {
        [self.start, self.ctrl1, self.ctrl2, self.end]
    }

    /// Replace all control points. Leave the curve unchanged if validation fails.
    pub fn set_control_points(&mut self, points: [P; 4]) -> Result<(), CurveConstructionError> {
        *self = Self::new(points[0], points[1], points[2], points[3])?;
        Ok(())
    }

    /// Return a curve with mapped control points. Reject non-finite results.
    /// Affine maps, such as translations, also map every point of the curve.
    pub fn try_map_control_points(
        &self,
        map: impl FnMut(P) -> P,
    ) -> Result<Self, CurveConstructionError> {
        let points = (*self).into_control_points().map(map);
        Self::new(points[0], points[1], points[2], points[3])
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
        CubicBezier {
            start: self.end,
            ctrl1: self.ctrl2,
            ctrl2: self.ctrl1,
            end: self.start,
        }
    }

    /// Approximate arc length with `segments` straight lines at equal parameter steps.
    pub fn arclen_polyline(&self, segments: usize) -> P::Scalar
    where
        P: PointNorm,
    {
        crate::Bezier {
            control_points: [self.start, self.ctrl1, self.ctrl2, self.end],
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
            control_points: [self.start, self.ctrl1, self.ctrl2, self.end],
        }
        .arclen_between(start, end, accuracy)
    }

    /// Split the curve at `t` into two sub-curves.
    pub fn split(&self, t: P::Scalar) -> (Self, Self) {
        let (left, right) = crate::Bezier::from(*self).split(t);
        (Self::from(left), Self::from(right))
    }

    /// Return the derivative as a quadratic Bezier curve.
    /// Evaluate it with `derivative().eval(t)` or use [`Self::derivative_at`].
    pub fn derivative(&self) -> QuadraticBezier<P> {
        let three = <P::Scalar as NumCast>::from(3.0).unwrap();
        QuadraticBezier {
            start: (self.ctrl1 - self.start) * three,
            ctrl: (self.ctrl2 - self.ctrl1) * three,
            end: (self.end - self.ctrl2) * three,
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

    fn are_points_colinear(&self, tolerance: P::Scalar) -> bool
    where
        P: PointNorm,
    {
        let line = self.baseline();
        line.distance_to_point(self.ctrl1) <= tolerance
            && line.distance_to_point(self.ctrl2) <= tolerance
    }

    /// Return whether every control point is within `abs(tolerance)` of the start.
    /// Zero requires equal coordinates. Non-finite tolerances return `false`.
    pub fn is_a_point(&self, tolerance: P::Scalar) -> bool
    where
        P: PointNorm,
    {
        crate::distance::is_a_point(&[self.start, self.ctrl1, self.ctrl2, self.end], tolerance)
    }

    /// Return up to three real roots of `a*t^3 + b*t^2 + c*t + d`.
    /// Uses Cardano's method, adapted from [Lyon](https://github.com/nical/lyon).
    /// Call once for each coordinate.
    #[allow(clippy::many_single_char_names)] // Match the polynomial's coefficient names.
    pub(crate) fn real_roots(
        &self,
        a: P::Scalar,
        b: P::Scalar,
        c: P::Scalar,
        d: P::Scalar,
    ) -> ArrayVec<[P::Scalar; 3]>
    where
        [P::Scalar; 3]: tinyvec::Array<Item = P::Scalar>,
    {
        let mut result = ArrayVec::new();
        let pi = <P::Scalar as NumCast>::from(core::f64::consts::PI).unwrap();
        let zero = <P::Scalar as NumCast>::from(0.0).unwrap();
        let two = <P::Scalar as NumCast>::from(2.0).unwrap();
        let three = <P::Scalar as NumCast>::from(3.0).unwrap();
        let four = <P::Scalar as NumCast>::from(4.0).unwrap();
        let nine = <P::Scalar as NumCast>::from(9.0).unwrap();
        let twenty_seven = <P::Scalar as NumCast>::from(27.0).unwrap();
        let fifty_four = <P::Scalar as NumCast>::from(54.0).unwrap();
        let coefficient_scale = a.abs().max(b.abs()).max(c.abs()).max(d.abs());

        // Treat negligible leading coefficients as lower-degree polynomials.
        if is_negligible(a, coefficient_scale) {
            if is_negligible(b, coefficient_scale) {
                if is_negligible(c, coefficient_scale) {
                    // The polynomial is constant and has no isolated roots.
                    return result;
                }
                result.push(-d / c);
                return result;
            }
            let c_squared = c * c;
            let four_bd = b * d * four;
            let delta = c_squared - four_bd;
            if is_negligible(delta, c_squared.abs() + four_bd.abs()) {
                result.push(-c / (b * two));
            } else if delta > zero {
                let sqrt_delta = delta.sqrt();
                result.push((-c - sqrt_delta) / (b * two));
                result.push((-c + sqrt_delta) / (b * two));
            }
            return result;
        }

        // Solve the cubic with Cardano's method.
        let frac_1_3 = <P::Scalar as NumCast>::from(1.0 / 3.0).unwrap();

        let bn = b / a;
        let cn = c / a;
        let dn = d / a;

        let delta0: P::Scalar = (cn * three - bn * bn) / nine;
        let delta1: P::Scalar =
            (bn * cn * nine - dn * twenty_seven - bn * bn * bn * two) / fifty_four;
        let delta_01: P::Scalar = delta0 * delta0 * delta0 + delta1 * delta1;

        if delta_01 >= zero {
            let delta_p_sqrt: P::Scalar = delta1 + delta_01.sqrt();
            let delta_m_sqrt: P::Scalar = delta1 - delta_01.sqrt();

            let s = delta_p_sqrt.signum() * delta_p_sqrt.abs().powf(frac_1_3);
            let t = delta_m_sqrt.signum() * delta_m_sqrt.abs().powf(frac_1_3);

            result.push(-bn * frac_1_3 + (s + t));

            // Don't add the repeated root when s + t == 0.
            let root_scale = s.abs().max(t.abs());
            if is_negligible(s - t, root_scale) && !is_negligible(s + t, root_scale) {
                result.push(-bn * frac_1_3 - (s + t) / two);
            }
        } else {
            let theta = (delta1 / (-delta0 * delta0 * delta0).sqrt()).acos();
            let two_sqrt_delta0 = (-delta0).sqrt() * two;
            result.push(two_sqrt_delta0 * Float::cos(theta * frac_1_3) - bn * frac_1_3);
            result
                .push(two_sqrt_delta0 * Float::cos((theta + two * pi) * frac_1_3) - bn * frac_1_3);
            result
                .push(two_sqrt_delta0 * Float::cos((theta + four * pi) * frac_1_3) - bn * frac_1_3);
        }

        result
    }

    /// Return parameters in `[0, 1]` where coordinate `axis` equals `value`.
    /// Solve the coordinate polynomial analytically.
    #[allow(dead_code)]
    fn solve_t_for_axis(&self, value: P::Scalar, axis: usize) -> ArrayVec<[P::Scalar; 3]>
    where
        P: PointIndex + PointNorm,
        [P::Scalar; 3]: tinyvec::Array<Item = P::Scalar>,
    {
        let mut result = ArrayVec::new();
        // Handle coincident controls and linear curves before the cubic solver.
        let zero = <P::Scalar as NumCast>::from(0.0).unwrap();
        if self.is_a_point(zero)
            || (self.are_points_colinear(zero) && (self.start - self.end).squared_norm() == zero)
        {
            return result;
        }
        let three = <P::Scalar as NumCast>::from(3.0).unwrap();
        let six = <P::Scalar as NumCast>::from(6.0).unwrap();
        let a = -self.start[axis] + self.ctrl1[axis] * three - self.ctrl2[axis] * three
            + self.end[axis];
        let b = self.start[axis] * three - self.ctrl1[axis] * six + self.ctrl2[axis] * three;
        let c = -self.start[axis] * three + self.ctrl1[axis] * three;
        let d = self.start[axis] - value;

        let roots = self.real_roots(a, b, c, d);
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
        [P::Scalar; 2]: tinyvec::Array<Item = P::Scalar>,
        [P::Scalar; 3]: tinyvec::Array<Item = P::Scalar>,
        [P::Scalar; 4]: tinyvec::Array<Item = P::Scalar>,
    {
        // Write the quadratic derivative in the power basis, a*t^2 + b*t + c.
        let zero = <P::Scalar as NumCast>::from(0.0).unwrap();
        let one = <P::Scalar as NumCast>::from(1.0).unwrap();
        let mut bounds = [(zero, zero); P::DIM];
        let derivative = self.derivative();
        let two = <P::Scalar as NumCast>::from(2.0).unwrap();
        let a: P = derivative.start + derivative.ctrl * -two + derivative.end;
        let b: P = derivative.start * -two + derivative.ctrl * two;
        let c: P = derivative.start;
        if !a.is_finite() || !b.is_finite() || !c.is_finite() {
            return crate::Bezier {
                control_points: [self.start, self.ctrl1, self.ctrl2, self.end],
            }
            .bounding_box();
        }

        // Coordinate extrema occur at derivative roots or curve endpoints.
        for dim in 0..P::DIM {
            let mut extrema: ArrayVec<[P::Scalar; 4]> = ArrayVec::new();
            extrema.extend(
                derivative
                    .real_roots(a[dim], b[dim], c[dim])
                    .iter()
                    .copied(),
            );
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

    /// Return a conservative box around the control points.
    pub fn control_box(&self) -> BoundingBox<P::Scalar, { P::DIM }>
    where
        P: PointIndex,
    {
        BoundingBox::from_points(&[self.start, self.ctrl1, self.ctrl2, self.end])
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ArrayPoint;
    use crate::CurveConstructionError;
    use crate::test_support::{assert_approx_equal, assert_point_approx_equal};
    use core::f64::consts::PI;

    const ABSOLUTE_TOLERANCE: f64 = 1e-12;
    const RELATIVE_TOLERANCE: f64 = 1e-12;

    #[test]
    fn rejects_non_finite_control_points() {
        let result = CubicBezier::new(
            ArrayPoint::new([0.0_f64]),
            ArrayPoint::new([1.0]),
            ArrayPoint::new([f64::NEG_INFINITY]),
            ArrayPoint::new([f64::NAN]),
        );

        assert_eq!(
            result,
            Err(CurveConstructionError::NonFiniteControlPoint { index: 2 })
        );
    }
    #[test]
    fn circle_approximation_error() {
        let circle = |p: ArrayPoint<f64, 2>| -> f64 {
            p.into_iter().map(|x| x * x).sum::<f64>().sqrt() - 1f64
        };

        // Use Mortensen's quarter-circle controls to minimize radial error.
        // Reference: http://spencermortensen.com/articles/bezier-circle/
        let c = 0.551915024494;
        let max_drift_perc = 0.019608; // radial drift percent
        let max_error = max_drift_perc * 0.01; // absolute max radial error

        let bezier_quadrant_1 = CubicBezier {
            start: ArrayPoint::new([0f64, 1f64]),
            ctrl1: ArrayPoint::new([c, 1f64]),
            ctrl2: ArrayPoint::new([1f64, c]),
            end: ArrayPoint::new([1f64, 0f64]),
        };
        let bezier_quadrant_2 = CubicBezier {
            start: ArrayPoint::new([1f64, 0f64]),
            ctrl1: ArrayPoint::new([1f64, -c]),
            ctrl2: ArrayPoint::new([c, -1f64]),
            end: ArrayPoint::new([0f64, -1f64]),
        };
        let bezier_quadrant_3 = CubicBezier {
            start: ArrayPoint::new([0f64, -1f64]),
            ctrl1: ArrayPoint::new([-c, -1f64]),
            ctrl2: ArrayPoint::new([-1f64, -c]),
            end: ArrayPoint::new([-1f64, 0f64]),
        };
        let bezier_quadrant_4 = CubicBezier {
            start: ArrayPoint::new([-1f64, 0f64]),
            ctrl1: ArrayPoint::new([-1f64, c]),
            ctrl2: ArrayPoint::new([-c, 1f64]),
            end: ArrayPoint::new([0f64, 1f64]),
        };
        let nsteps = 1000;
        for t in 0..=nsteps {
            let t = t as f64 * 1f64 / (nsteps as f64);

            let point = bezier_quadrant_1.eval(t);
            let contour = circle(point);
            assert!(contour.abs() <= max_error);

            let point = bezier_quadrant_2.eval(t);
            let contour = circle(point);
            assert!(contour.abs() <= max_error);

            let point = bezier_quadrant_3.eval(t);
            let contour = circle(point);
            assert!(contour.abs() <= max_error);

            let point = bezier_quadrant_4.eval(t);
            let contour = circle(point);
            assert!(contour.abs() <= max_error);
        }
    }

    #[test]
    fn circle_circumference_approximation() {
        // Approximate a unit circle with four cubic curves. Measure each with
        // 1,000 polyline segments. The comparison includes both curve-fitting
        // error and the polyline length error.
        let c = 0.551915024494;
        let max_error = 1e-2;
        let nsteps = 1e3 as usize;
        let pi = PI;
        let tau = 2. * pi;

        let bezier_quadrant_1 = CubicBezier {
            start: ArrayPoint::new([0f64, 1f64]),
            ctrl1: ArrayPoint::new([c, 1f64]),
            ctrl2: ArrayPoint::new([1f64, c]),
            end: ArrayPoint::new([1f64, 0f64]),
        };
        let bezier_quadrant_2 = CubicBezier {
            start: ArrayPoint::new([1f64, 0f64]),
            ctrl1: ArrayPoint::new([1f64, -c]),
            ctrl2: ArrayPoint::new([c, -1f64]),
            end: ArrayPoint::new([0f64, -1f64]),
        };
        let bezier_quadrant_3 = CubicBezier {
            start: ArrayPoint::new([0f64, -1f64]),
            ctrl1: ArrayPoint::new([-c, -1f64]),
            ctrl2: ArrayPoint::new([-1f64, -c]),
            end: ArrayPoint::new([-1f64, 0f64]),
        };
        let bezier_quadrant_4 = CubicBezier {
            start: ArrayPoint::new([-1f64, 0f64]),
            ctrl1: ArrayPoint::new([-1f64, c]),
            ctrl2: ArrayPoint::new([-c, 1f64]),
            end: ArrayPoint::new([0f64, 1f64]),
        };
        let circumference = bezier_quadrant_1.arclen_polyline(nsteps)
            + bezier_quadrant_2.arclen_polyline(nsteps)
            + bezier_quadrant_3.arclen_polyline(nsteps)
            + bezier_quadrant_4.arclen_polyline(nsteps);
        assert!(((tau + max_error) > circumference) && ((tau - max_error) < circumference));
    }

    #[test]
    fn eval_equivalence_casteljau() {
        // Both evaluation entry points use de Casteljau interpolation.
        let bezier = CubicBezier::new(
            ArrayPoint::new([0f64, 1.77f64]),
            ArrayPoint::new([1.1f64, -1f64]),
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
    fn split_equivalence() {
        let bezier = CubicBezier {
            start: ArrayPoint::new([0f64, 1.77f64]),
            ctrl1: ArrayPoint::new([2.9f64, 0f64]),
            ctrl2: ArrayPoint::new([4.3f64, 3f64]),
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
    fn cubic_api_parity() {
        let bezier = CubicBezier::new(
            ArrayPoint::new([0f64, 0f64]),
            ArrayPoint::new([1f64, 0f64]),
            ArrayPoint::new([2f64, 0f64]),
            ArrayPoint::new([3f64, 0f64]),
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

        let t = bezier.t_at_length(1.5, 1e-12).unwrap();
        assert!((t - 0.5).abs() < 1e-6);

        let p = bezier.point_at_length(1.5, 1e-12).unwrap();
        assert_point_approx_equal(p, ArrayPoint::new([1.5, 0.0]), 1e-9, 1e-9);
    }

    #[test]
    fn bounding_box_contains() {
        // Allow for floating-point error when comparing samples with the box.
        let bezier = CubicBezier {
            start: ArrayPoint::new([0f64, 1.77f64]),
            ctrl1: ArrayPoint::new([2.9f64, 0f64]),
            ctrl2: ArrayPoint::new([4.3f64, -3f64]),
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
    fn small_scale_cubic_remains_linear_with_unit_tangent() {
        let scale = 1e-20_f64;
        let bezier = CubicBezier::new(
            ArrayPoint::new([0.0, 0.0]),
            ArrayPoint::new([scale / 3.0, 0.0]),
            ArrayPoint::new([scale * 2.0 / 3.0, 0.0]),
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
    fn cubic_roots_are_scale_invariant() {
        let bezier = CubicBezier::new(
            ArrayPoint::new([0.0_f64]),
            ArrayPoint::new([0.0]),
            ArrayPoint::new([0.0]),
            ArrayPoint::new([0.0]),
        )
        .unwrap();
        let mut roots = bezier.real_roots(1e-30, -6e-30, 11e-30, -6e-30);
        roots.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap());

        assert_eq!(roots.len(), 3);
        for (actual, expected) in roots.into_iter().zip([1.0, 2.0, 3.0]) {
            assert_approx_equal(actual, expected, ABSOLUTE_TOLERANCE, RELATIVE_TOLERANCE);
        }
    }

    #[test]
    fn distance_to_point() {
        // A cubic with four controls requires eight knots.
        let curve = CubicBezier {
            start: ArrayPoint::new([0f64, 1.77f64]),
            ctrl1: ArrayPoint::new([1.1f64, -1f64]),
            ctrl2: ArrayPoint::new([4.3f64, 3f64]),
            end: ArrayPoint::new([3.2f64, -4f64]),
        };
        assert!(
            curve.distance_to_point(ArrayPoint::new([-5.1, -5.6]))
                > curve.distance_to_point(ArrayPoint::new([5.1, 5.6]))
        );
    }
}
