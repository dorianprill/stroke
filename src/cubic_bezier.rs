//! Cubic Bezier curve specialization.

use super::{
    ArcLengthError, ArrayVec, BoundingBox, LineSegment, Point, PointDot, PointIndex, PointNorm,
    QuadraticBezier,
};
use crate::error::{CurveConstructionError, validate_control_points};
use crate::numeric::is_negligible;
use num_traits::{Float, NumCast};

const DEFAULT_DISTANCE_STEPS: usize = 64;
const LOCAL_SEARCH_ITERS: usize = 16;

/// Cubic Bezier curve defined by four points: start, two control points, and end.
///
/// The curve is defined by:
/// ```text
/// ∀ t ∈ [0..1], P(t) = (1 - t)³ * start
///            + 3 * (1 - t)² * t * ctrl1
///            + 3 * t² * (1 - t) * ctrl2
///            + t³ * end
/// ```
///
/// Methods that need component access or norms add `PointIndex`/`PointNorm`
/// bounds as required.
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

    /// Evaluate a `CubicBezier` curve at t by direct evaluation of the polynomial (not numerically stable)
    pub fn eval(&self, t: P::Scalar) -> P {
        let one = <P::Scalar as NumCast>::from(1.0).unwrap();
        let three = <P::Scalar as NumCast>::from(3.0).unwrap();
        let one_t = one - t;
        let one_t2 = one_t * one_t;
        let one_t3 = one_t2 * one_t;
        let t2 = t * t;
        let t3 = t2 * t;

        self.start * one_t3
            + self.ctrl1 * (t * one_t2 * three)
            + self.ctrl2 * (t2 * one_t * three)
            + self.end * t3
    }

    /// Evaluate a `CubicBezier` curve at t using the numerically stable De Casteljau algorithm
    pub fn eval_casteljau(&self, t: P::Scalar) -> P {
        // unrolled de casteljau algorithm
        // _1ab is the first iteration from first (a) to second (b) control point and so on
        let ctrl_1ab = self.start + (self.ctrl1 - self.start) * t;
        let ctrl_1bc = self.ctrl1 + (self.ctrl2 - self.ctrl1) * t;
        let ctrl_1cd = self.ctrl2 + (self.end - self.ctrl2) * t;
        // second iteration
        let ctrl_2ab = ctrl_1ab + (ctrl_1bc - ctrl_1ab) * t;
        let ctrl_2bc = ctrl_1bc + (ctrl_1cd - ctrl_1bc) * t;
        // third iteration, return final point on the curve ctrl_3ab
        ctrl_2ab + (ctrl_2bc - ctrl_2ab) * t
    }

    /// Borrow the control points in curve order.
    pub fn control_points(&self) -> core::array::IntoIter<&P, 4> {
        [&self.start, &self.ctrl1, &self.ctrl2, &self.end].into_iter()
    }

    /// Consume the curve, returning its control-point array.
    pub fn into_control_points(self) -> [P; 4] {
        [self.start, self.ctrl1, self.ctrl2, self.end]
    }

    /// Replace all control points atomically, leaving this curve unchanged on error.
    pub fn set_control_points(&mut self, points: [P; 4]) -> Result<(), CurveConstructionError> {
        *self = Self::new(points[0], points[1], points[2], points[3])?;
        Ok(())
    }

    /// Map control points and validate the result without modifying this curve.
    /// This maps the entire curve geometrically only when `map` is affine.
    pub fn try_map_control_points(
        &self,
        map: impl FnMut(P) -> P,
    ) -> Result<Self, CurveConstructionError> {
        let points = (*self).into_control_points().map(map);
        Self::new(points[0], points[1], points[2], points[3])
    }

    /// Evaluate the derivative with respect to the unit-domain parameter.
    pub fn derivative_at(&self, t: P::Scalar) -> P {
        crate::Bezier::from(*self).derivative_at(t)
    }

    /// Return the derivative magnitude with respect to the unit-domain parameter.
    pub fn speed(&self, t: P::Scalar) -> P::Scalar
    where
        P: PointNorm,
    {
        self.derivative_at(t).norm()
    }

    /// Extract a finite ordered range in `[0, 1]`, reparameterized to `[0, 1]`.
    /// Equal range endpoints yield a constant curve.
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

    /// Approximate the arc length by flattening the curve into `segments` lines.
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
        crate::Bezier {
            control_points: [self.start, self.ctrl1, self.ctrl2, self.end],
        }
        .arclen_between(start, end, accuracy)
    }

    /// Split the curve at `t` into two sub-curves.
    pub fn split(&self, t: P::Scalar) -> (Self, Self) {
        // unrolled de casteljau algorithm
        // _1ab is the first iteration from first (a) to second (b) control point and so on
        let ctrl_1ab = self.start + (self.ctrl1 - self.start) * t;
        let ctrl_1bc = self.ctrl1 + (self.ctrl2 - self.ctrl1) * t;
        let ctrl_1cd = self.ctrl2 + (self.end - self.ctrl2) * t;
        // second iteration
        let ctrl_2ab = ctrl_1ab + (ctrl_1bc - ctrl_1ab) * t;
        let ctrl_2bc = ctrl_1bc + (ctrl_1cd - ctrl_1bc) * t;
        // third iteration, final point on the curve
        let ctrl_3ab = ctrl_2ab + (ctrl_2bc - ctrl_2ab) * t;

        (
            CubicBezier {
                start: self.start,
                ctrl1: ctrl_1ab,
                ctrl2: ctrl_2ab,
                end: ctrl_3ab,
            },
            CubicBezier {
                start: ctrl_3ab,
                ctrl1: ctrl_2bc,
                ctrl2: ctrl_1cd,
                end: self.end,
            },
        )
    }

    /// Return the derivative curve.
    /// The derivative is also a bezier curve but of degree n-1 (cubic->quadratic)
    /// Since it returns the derivative function, `eval()` needs to be called separately
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
        let one = <P::Scalar as NumCast>::from(1.0).unwrap();
        let zero = <P::Scalar as NumCast>::from(0.0).unwrap();
        let dir = self.derivative().eval(t);
        let len = dir.norm();
        if len == zero { dir } else { dir * (one / len) }
    }

    /// Return the curvature magnitude at `t`.
    pub fn curvature(&self, t: P::Scalar) -> P::Scalar
    where
        P: PointNorm + PointDot,
    {
        let zero = <P::Scalar as NumCast>::from(0.0).unwrap();
        let v = self.derivative().eval(t);
        let a = self.derivative().derivative().eval(t);
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
    {
        let zero = <P::Scalar as NumCast>::from(0.0).unwrap();
        let v = self.derivative().eval(t);
        let v2 = v.squared_norm();
        if v2 == zero {
            return None;
        }
        let a = self.derivative().derivative().eval(t);
        let dot = PointDot::dot(&v, &a);
        let a_perp = a - v * (dot / v2);
        let n2 = a_perp.squared_norm();
        if n2 == zero {
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
        let t2 = t * t;
        let three = <P::Scalar as NumCast>::from(3.0).unwrap();
        let six = <P::Scalar as NumCast>::from(6.0).unwrap();
        let nine = <P::Scalar as NumCast>::from(9.0).unwrap();
        let twelve = <P::Scalar as NumCast>::from(12.0).unwrap();
        let c0 = t * -three + t * six - three;
        let c1 = t2 * nine - t * twelve + three;
        let c2 = t2 * -nine + t * six;
        let c3 = t2 * three;

        self.start[axis] * c0 + self.ctrl1[axis] * c1 + self.ctrl2[axis] * c2 + self.end[axis] * c3
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
        Ok(self.eval_casteljau(self.t_at_length(s, accuracy)?))
    }

    /// Returns the line segment formed by the curve's start and end points.
    pub fn baseline(&self) -> LineSegment<P> {
        LineSegment {
            start: self.start,
            end: self.end,
        }
    }

    /// Checks whether, given some tolerance, the curve can be considered linear.
    pub fn is_linear(&self, tolerance: P::Scalar) -> bool
    where
        P: PointNorm,
    {
        let tolerance = tolerance.abs();
        if (self.start - self.end).squared_norm() <= tolerance * tolerance {
            return false;
        }
        // else check if ctrl points lie on baseline
        self.are_points_colinear(tolerance)
    }

    fn are_points_colinear(&self, tolerance: P::Scalar) -> bool
    where
        P: PointNorm,
    {
        let line = self.baseline();
        line.distance_to_point(self.ctrl1) <= tolerance
            && line.distance_to_point(self.ctrl2) <= tolerance
    }

    // Returs if the whole set of control points can be considered one singular point
    // given some tolerance.
    // TODO use machine epsilon vs squared norm OK?
    /// Checks whether, given some tolerance, all control points are coincident.
    pub fn is_a_point(&self, tolerance: P::Scalar) -> bool
    where
        P: PointNorm,
    {
        let tolerance = tolerance.abs();
        let tolerance_squared = tolerance * tolerance;
        // Use <= so that tolerance can be zero.
        (self.start - self.end).squared_norm() <= tolerance_squared
            && (self.start - self.ctrl1).squared_norm() <= tolerance_squared
            && (self.end - self.ctrl2).squared_norm() <= tolerance_squared
    }

    /// Compute the real roots of the cubic bezier function with
    /// parameters of the form a*t^3 + b*t^2 + c*t + d for each dimension
    /// using cardano's algorithm (code adapted from github.com/nical/lyon)
    /// returns an `ArrayVec` of the present roots (max 3)
    #[allow(clippy::many_single_char_names)] // this is math, get over it
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

        // check if can be handled below cubic order
        if is_negligible(a, coefficient_scale) {
            if is_negligible(b, coefficient_scale) {
                if is_negligible(c, coefficient_scale) {
                    // no solutions
                    return result;
                }
                // is linear equation
                result.push(-d / c);
                return result;
            }
            // is quadratic equation
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

        // is cubic equation -> use cardano's algorithm
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

    /// Solves the cubic bezier function given a particular coordinate axis value
    /// by solving the roots for the axis functions
    /// Parameters:
    /// value: the coordinate value on the particular axis
    /// axis: the index of the axis
    /// Returns those roots of the function that are in the interval [0.0, 1.0].
    #[allow(dead_code)]
    fn solve_t_for_axis(&self, value: P::Scalar, axis: usize) -> ArrayVec<[P::Scalar; 3]>
    where
        P: PointIndex + PointNorm,
        [P::Scalar; 3]: tinyvec::Array<Item = P::Scalar>,
    {
        let mut result = ArrayVec::new();
        // check if all points are the same or if the curve is really just a line
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
        // calculate coefficients for the derivative: at^2 + bt + c
        // from the expansion of the cubic bezier curve: sum_i=0_to_3( binomial(3, i) * t^i * (1-t)^(n-i) )
        // yields coeffcients
        // po: [1, -2,  1]
        // p1: [0,  2, -2]
        // p2: [0,  0,  1]
        //      c   b   a
        let zero = <P::Scalar as NumCast>::from(0.0).unwrap();
        let one = <P::Scalar as NumCast>::from(1.0).unwrap();
        let mut bounds = [(zero, zero); P::DIM];
        let derivative = self.derivative();
        // calculate coefficients for derivative
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

        // calculate roots for t over x axis and plug them into the bezier function
        //  to get x,y values (make vec 2 bigger for t=0,t=1 values)
        // loop over any of the points dimensions (they're all the same)
        for dim in 0..P::DIM {
            let mut extrema: ArrayVec<[P::Scalar; 4]> = ArrayVec::new();
            extrema.extend(
                derivative
                    .real_roots(a[dim], b[dim], c[dim])
                    .iter()
                    .copied(),
            );
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
            // .unwrap() is ok as it can never be empty as it always at least contains the endpoints
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
    use crate::CurveConstructionError;
    use crate::PointN;
    use crate::test_support::{assert_approx_equal, assert_point_approx_equal};
    use core::f64::consts::PI;

    const ABSOLUTE_TOLERANCE: f64 = 1e-12;
    const RELATIVE_TOLERANCE: f64 = 1e-12;

    #[test]
    fn rejects_non_finite_control_points() {
        let result = CubicBezier::new(
            PointN::new([0.0_f64]),
            PointN::new([1.0]),
            PointN::new([f64::NEG_INFINITY]),
            PointN::new([f64::NAN]),
        );

        assert_eq!(
            result,
            Err(CurveConstructionError::NonFiniteControlPoint { index: 2 })
        );
    }
    #[test]
    fn circle_approximation_error() {
        // define closure for unit circle
        let circle =
            |p: PointN<f64, 2>| -> f64 { p.into_iter().map(|x| x * x).sum::<f64>().sqrt() - 1f64 };

        // define control points for 4 bezier segments
        // control points are chosen for minimum radial distance error
        // according to: http://spencermortensen.com/articles/bezier-circle/
        // TODO don't hardcode values
        let c = 0.551915024494;
        let max_drift_perc = 0.019608; // radial drift percent
        let max_error = max_drift_perc * 0.01; // absolute max radial error

        let bezier_quadrant_1 = CubicBezier {
            start: PointN::new([0f64, 1f64]),
            ctrl1: PointN::new([c, 1f64]),
            ctrl2: PointN::new([1f64, c]),
            end: PointN::new([1f64, 0f64]),
        };
        let bezier_quadrant_2 = CubicBezier {
            start: PointN::new([1f64, 0f64]),
            ctrl1: PointN::new([1f64, -c]),
            ctrl2: PointN::new([c, -1f64]),
            end: PointN::new([0f64, -1f64]),
        };
        let bezier_quadrant_3 = CubicBezier {
            start: PointN::new([0f64, -1f64]),
            ctrl1: PointN::new([-c, -1f64]),
            ctrl2: PointN::new([-1f64, -c]),
            end: PointN::new([-1f64, 0f64]),
        };
        let bezier_quadrant_4 = CubicBezier {
            start: PointN::new([-1f64, 0f64]),
            ctrl1: PointN::new([-1f64, c]),
            ctrl2: PointN::new([-c, 1f64]),
            end: PointN::new([0f64, 1f64]),
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
        // define control points for 4 cubic bezier segments to best approximate a unit circle
        // control points are chosen for minimum radial distance error, see circle_approximation_error() in this file
        // given this, the circumference will also be close to 2*pi
        // (remember arclen also works by linear approximation, not the true integral, so we have to accept error)!
        // This approximation is unfeasable if desired accuracy is greater than ~2 decimal places (at 1000 steps)
        // TODO don't hardcode values, solve for them
        let c = 0.551915024494;
        let max_error = 1e-2;
        let nsteps = 1e3 as usize;
        let pi = PI;
        let tau = 2. * pi;

        let bezier_quadrant_1 = CubicBezier {
            start: PointN::new([0f64, 1f64]),
            ctrl1: PointN::new([c, 1f64]),
            ctrl2: PointN::new([1f64, c]),
            end: PointN::new([1f64, 0f64]),
        };
        let bezier_quadrant_2 = CubicBezier {
            start: PointN::new([1f64, 0f64]),
            ctrl1: PointN::new([1f64, -c]),
            ctrl2: PointN::new([c, -1f64]),
            end: PointN::new([0f64, -1f64]),
        };
        let bezier_quadrant_3 = CubicBezier {
            start: PointN::new([0f64, -1f64]),
            ctrl1: PointN::new([-c, -1f64]),
            ctrl2: PointN::new([-1f64, -c]),
            end: PointN::new([-1f64, 0f64]),
        };
        let bezier_quadrant_4 = CubicBezier {
            start: PointN::new([-1f64, 0f64]),
            ctrl1: PointN::new([-1f64, c]),
            ctrl2: PointN::new([-c, 1f64]),
            end: PointN::new([0f64, 1f64]),
        };
        let circumference = bezier_quadrant_1.arclen_polyline(nsteps)
            + bezier_quadrant_2.arclen_polyline(nsteps)
            + bezier_quadrant_3.arclen_polyline(nsteps)
            + bezier_quadrant_4.arclen_polyline(nsteps);
        //dbg!(circumference);
        //dbg!(tau);
        assert!(((tau + max_error) > circumference) && ((tau - max_error) < circumference));
    }

    #[test]
    fn eval_equivalence_casteljau() {
        // all eval methods should be approximately equivalent for well defined test cases
        // and not equivalent where numerical stability becomes an issue for normal eval
        let bezier = CubicBezier::new(
            PointN::new([0f64, 1.77f64]),
            PointN::new([1.1f64, -1f64]),
            PointN::new([4.3f64, 3f64]),
            PointN::new([3.2f64, -4f64]),
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
        // chose some arbitrary control points and construct a cubic bezier
        let bezier = CubicBezier {
            start: PointN::new([0f64, 1.77f64]),
            ctrl1: PointN::new([2.9f64, 0f64]),
            ctrl2: PointN::new([4.3f64, 3f64]),
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
            assert_point_approx_equal(
                bezier.eval(t / 2.0),
                left.eval(t),
                ABSOLUTE_TOLERANCE,
                RELATIVE_TOLERANCE,
            );
            // right
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
            PointN::new([0f64, 0f64]),
            PointN::new([1f64, 0f64]),
            PointN::new([2f64, 0f64]),
            PointN::new([3f64, 0f64]),
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
        assert_point_approx_equal(p, PointN::new([1.5, 0.0]), 1e-9, 1e-9);
    }

    #[test]
    fn bounding_box_contains() {
        // check if bounding box for a curve contains all points (with some approximation error)
        let bezier = CubicBezier {
            start: PointN::new([0f64, 1.77f64]),
            ctrl1: PointN::new([2.9f64, 0f64]),
            ctrl2: PointN::new([4.3f64, -3f64]),
            end: PointN::new([3.2f64, 4f64]),
        };

        let bounds = bezier.bounding_box();

        let max_err = 1e-2;

        let nsteps: usize = 100;
        for t in 0..=nsteps {
            let t = t as f64 * 1f64 / (nsteps as f64);
            let p = bezier.eval_casteljau(t);
            //dbg!(t);
            //dbg!(p);
            //dbg!(xmin-max_err, ymin-max_err, xmax+max_err, ymax+max_err);
            for (idx, axis) in p.into_iter().enumerate() {
                assert!((axis >= (bounds[idx].0 - max_err)) && (axis <= (bounds[idx].1 + max_err)))
            }
        }
    }

    #[test]
    fn small_scale_cubic_remains_linear_with_unit_tangent() {
        let scale = 1e-20_f64;
        let bezier = CubicBezier::new(
            PointN::new([0.0, 0.0]),
            PointN::new([scale / 3.0, 0.0]),
            PointN::new([scale * 2.0 / 3.0, 0.0]),
            PointN::new([scale, 0.0]),
        )
        .unwrap();

        assert!(bezier.is_linear(scale * 1e-6));
        assert_point_approx_equal(
            bezier.tangent(0.5),
            PointN::new([1.0, 0.0]),
            ABSOLUTE_TOLERANCE,
            RELATIVE_TOLERANCE,
        );
    }

    #[test]
    fn cubic_roots_are_scale_invariant() {
        let bezier = CubicBezier::new(
            PointN::new([0.0_f64]),
            PointN::new([0.0]),
            PointN::new([0.0]),
            PointN::new([0.0]),
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
        // degree 3, 4 control points => 4+3+1=8 knots
        let curve = CubicBezier {
            start: PointN::new([0f64, 1.77f64]),
            ctrl1: PointN::new([1.1f64, -1f64]),
            ctrl2: PointN::new([4.3f64, 3f64]),
            end: PointN::new([3.2f64, -4f64]),
        };
        assert!(
            curve.distance_to_point(PointN::new([-5.1, -5.6]))
                > curve.distance_to_point(PointN::new([5.1, 5.6]))
        );
    }
}
