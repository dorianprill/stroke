//! Const-generic B-spline curves.

use core::slice::*;

use super::{Point, PointDot, PointIndex, PointNorm};
use crate::find_root::FindRoot;
use crate::roots::{RootFindingError, root_newton_raphson};
use num_traits::{Float, NumCast};

const ROOT_TOLERANCE: f64 = 1e-6;
const MAX_ROOT_ITER: usize = 64;
const DEFAULT_DISTANCE_STEPS: usize = 64;
const DEFAULT_LENGTH_STEPS: usize = 64;
const MAX_LENGTH_ITERS: usize = 32;
const LOCAL_SEARCH_ITERS: usize = 16;

/// Errors returned by B-spline construction or evaluation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum BSplineError {
    /// The curve has fewer control points than required by its degree.
    TooFewControlPoints,
    /// The knot vector has fewer entries than required by the degree and control points.
    TooFewKnots,
    /// The knot vector has more entries than required by the degree and control points.
    TooManyKnots,
    /// The knot vector is not non-decreasing.
    DecreasingKnots,
    /// The evaluation parameter is outside the knot domain.
    KnotDomainViolation,
    /// The de Boor evaluation encountered a NaN alpha value.
    EvaluationAlphaIsNaN,
}

#[derive(Debug, Clone, Copy, PartialEq)]
enum KnotVectorKind {
    Unclamped,
    ClampedStart,
    ClampedEnd,
    Clamped,
}

impl Default for KnotVectorKind {
    fn default() -> Self {
        KnotVectorKind::Unclamped
    }
}

impl core::fmt::Display for BSplineError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            BSplineError::TooFewControlPoints => write!(f, "Too few control points"),
            BSplineError::TooFewKnots => write!(f, "Insufficient number of knots"),
            BSplineError::TooManyKnots => write!(f, "Excess number of knots"),
            BSplineError::DecreasingKnots => write!(f, "Knots must be non-decreasing in value"),
            BSplineError::KnotDomainViolation => write!(f, "Parameter outside of knot span"),
            BSplineError::EvaluationAlphaIsNaN => {
                write!(f, "Cox-deBoor evaluation yields alpha NaN value")
            }
        }
    }
}

/// Const-generic B-spline curve with fixed control-point and knot counts.
///
/// The scalar type used for the knots and interpolation is `P::Scalar`.
///
/// # Parameters
/// - `P`: point type implementing [`Point`](crate::point::Point).
/// - `C`: number of control points.
/// - `K`: number of knots (`K = C + D + 1`).
/// - `D`: degree (order - 1).
#[derive(Clone, Copy)]
pub struct BSpline<P, const K: usize, const C: usize, const D: usize>
where
    P: Point,
{
    /// Knot vector
    knots: [P::Scalar; K],
    /// Knot Vector Kind (Clamped or Unclamped)
    knot_kind: KnotVectorKind,
    /// Control points
    control_points: [P; C],
}

impl<P, const K: usize, const C: usize, const D: usize> Default for BSpline<P, K, C, D>
where
    P: Point + Default,
    P::Scalar: Default,
{
    fn default() -> Self {
        BSpline {
            knots: core::array::from_fn(|_| P::Scalar::default()),
            knot_kind: KnotVectorKind::default(),
            control_points: core::array::from_fn(|_| P::default()),
        }
    }
}

impl<P, const K: usize, const C: usize, const D: usize> BSpline<P, { K }, { C }, { D }>
where
    P: Point,
{
    /// Create a new B-spline curve that interpolates
    /// the `control_points` using a piecewise polynomial of `degree` within intervals specified by the `knots`.
    /// The knots _must_ be sorted in non-decreasing order; the constructor rejects decreasing knots.
    /// The degree is defined as `curve_order - 1`.
    /// Desired curve must have a valid number of control points and knots in relation to its degree or the constructor will return an error.
    /// A B-Spline curve requires at least one more control point than the degree (`control_points.len() >
    /// degree`) and the number of knots should be equal to `control_points.len() + degree + 1`.
    /// K = C + D + 1
    pub fn new(
        knots: [P::Scalar; K],
        control_points: [P; C],
    ) -> Result<BSpline<P, { K }, { C }, { D }>, BSplineError> {
        if control_points.len() <= { D } {
            Err(BSplineError::TooFewControlPoints)
        // n_knots = n_control_points + 1 + degree
        } else if knots.len() < control_points.len() + 1 + D {
            Err(BSplineError::TooFewKnots)
        } else if knots.len() > control_points.len() + 1 + D {
            Err(BSplineError::TooManyKnots)
        // if knots array is not non-decreasing
        } else if knots.windows(2).any(|w| w[0] > w[1]) {
            Err(BSplineError::DecreasingKnots)
        } else {
            let is_clamped_front = knots[0..=D]
                .iter()
                .all(|&knot| (knot - knots[0]).abs() < P::Scalar::epsilon());
            let is_clamped_back = knots[knots.len() - D - 1..]
                .iter()
                .all(|&knot| (knot - knots[knots.len() - 1]).abs() < P::Scalar::epsilon());

            let knot_kind = match (is_clamped_front, is_clamped_back) {
                (false, false) => KnotVectorKind::Unclamped,
                (true, false) => KnotVectorKind::ClampedStart,
                (false, true) => KnotVectorKind::ClampedEnd,
                (true, true) => KnotVectorKind::Clamped,
            };

            Ok(BSpline {
                knots,
                knot_kind,
                control_points,
            })
        }
    }

    /// Compute a point on the curve at `t` using iterative de boor algorithm.
    /// The parameter **must** be in the inclusive range of values returned
    /// by `knot_domain()`. If `t` is out of bounds a KnotDomainViolation error is returned.
    pub fn eval(&self, t: P::Scalar) -> Result<P, BSplineError>
    where
        [(); D + 1]: Sized,
    {
        let (knot_start, knot_end) = self.knot_domain();

        // Ensure t is within the knot domain
        if t < knot_start || t > knot_end {
            return Err(BSplineError::KnotDomainViolation);
        }

        // Calculate the point on the B-Spline at `t` using de Boor's algorithm
        let p = self.de_boor_iterative(t)?;
        Ok(p)
    }

    /// Return the start point of the curve.
    pub fn start(&self) -> Result<P, BSplineError>
    where
        [(); D + 1]: Sized,
    {
        let (kmin, _) = self.knot_domain();
        self.eval(kmin)
    }

    /// Return the end point of the curve.
    pub fn end(&self) -> Result<P, BSplineError>
    where
        [(); D + 1]: Sized,
    {
        let (_, kmax) = self.knot_domain();
        self.eval(kmax)
    }

    /// Return a curve with reversed direction.
    pub fn reverse(&self) -> Self {
        let control_points = core::array::from_fn(|i| self.control_points[C - 1 - i]);
        let knot_min = self.knots[0];
        let knot_max = self.knots[K - 1];
        let knots = core::array::from_fn(|i| knot_min + knot_max - self.knots[K - 1 - i]);
        let knot_kind = match self.knot_kind {
            KnotVectorKind::ClampedStart => KnotVectorKind::ClampedEnd,
            KnotVectorKind::ClampedEnd => KnotVectorKind::ClampedStart,
            other => other,
        };

        BSpline {
            knots,
            knot_kind,
            control_points,
        }
    }

    /// Returns an iterator over the control points.
    pub fn control_points(&self) -> Iter<'_, P> {
        self.control_points.iter()
    }

    /// Returns an iterator over the knots.
    pub fn knots(&self) -> Iter<'_, P::Scalar> {
        self.knots.iter()
    }

    /// Returns the knot domain of the B-Spline.
    /// The knot domain is the range of values over which the B-Spline is defined.
    /// The knot domain is defined as the minimum and maximum knot values.
    /// For a clamped B-Spline the first and last knot has multiplicity D+1.
    /// For an unclamped B-Spline the first and last knot has multiplicity 1.
    pub fn knot_domain(&self) -> (P::Scalar, P::Scalar) {
        match self.knot_kind {
            KnotVectorKind::Clamped => self.knot_domain_clamped(),
            KnotVectorKind::ClampedStart => self.knot_domain_clamped_start(),
            KnotVectorKind::ClampedEnd => self.knot_domain_clamped_end(),
            KnotVectorKind::Unclamped => self.knot_domain_unclamped(),
        }
    }

    /// Return true if `t` lies within the inclusive knot domain.
    ///
    /// # Examples
    /// ```rust
    /// use stroke::{BSpline, PointN};
    ///
    /// let knots: [f64; 4] = [0.0, 0.0, 1.0, 1.0];
    /// let points = [PointN::new([0.0]), PointN::new([2.0])];
    /// let curve: BSpline<PointN<f64, 1>, 4, 2, 1> = BSpline::new(knots, points).unwrap();
    ///
    /// assert!(curve.domain_contains(0.5));
    /// assert!(!curve.domain_contains(-0.1));
    /// ```
    pub fn domain_contains(&self, t: P::Scalar) -> bool {
        let (kmin, kmax) = self.knot_domain();
        t >= kmin && t <= kmax
    }

    /// Return the knot span index for `t`.
    ///
    /// The span index `i` satisfies `knots[i] <= t < knots[i + 1]`, with the
    /// conventional special-case that `t == domain_max` maps to the last span.
    /// Use this value when evaluating basis functions for `t`.
    ///
    /// # Examples
    /// ```rust
    /// use stroke::{BSpline, PointN};
    ///
    /// let knots: [f64; 4] = [0.0, 0.0, 1.0, 1.0];
    /// let points = [PointN::new([0.0]), PointN::new([2.0])];
    /// let curve: BSpline<PointN<f64, 1>, 4, 2, 1> = BSpline::new(knots, points).unwrap();
    ///
    /// let span = curve.knot_span(0.25).unwrap();
    /// assert_eq!(span, 1);
    /// ```
    pub fn knot_span(&self, t: P::Scalar) -> Result<usize, BSplineError> {
        if !self.domain_contains(t) {
            return Err(BSplineError::KnotDomainViolation);
        }
        self.knot_span_start_for_t(t)
            .ok_or(BSplineError::KnotDomainViolation)
    }

    /// Evaluate the non-zero basis functions for a given `span` and `t`.
    ///
    /// Returns the `D + 1` basis values `N_{i,p}(t)` for the active span,
    /// ordered from `i = span - D` through `i = span`.
    /// The caller should pass a `span` obtained from `knot_span(t)`.
    ///
    /// # Examples
    /// ```rust
    /// use stroke::{BSpline, PointN};
    ///
    /// let knots: [f64; 4] = [0.0, 0.0, 1.0, 1.0];
    /// let points = [PointN::new([0.0]), PointN::new([2.0])];
    /// let curve: BSpline<PointN<f64, 1>, 4, 2, 1> = BSpline::new(knots, points).unwrap();
    ///
    /// let t = 0.25;
    /// let span = curve.knot_span(t).unwrap();
    /// let basis = curve.basis_functions(span, t).unwrap();
    /// assert!((basis[0] - 0.75).abs() < 1e-6);
    /// assert!((basis[1] - 0.25).abs() < 1e-6);
    /// ```
    pub fn basis_functions(
        &self,
        span: usize,
        t: P::Scalar,
    ) -> Result<[P::Scalar; D + 1], BSplineError>
    where
        [(); D + 1]: Sized,
    {
        if !self.domain_contains(t) {
            return Err(BSplineError::KnotDomainViolation);
        }
        if span < D || span + D >= K {
            return Err(BSplineError::KnotDomainViolation);
        }

        let zero = <P::Scalar as NumCast>::from(0.0).unwrap();
        let one = <P::Scalar as NumCast>::from(1.0).unwrap();
        let mut left = [zero; D + 1];
        let mut right = [zero; D + 1];
        let mut basis = [zero; D + 1];
        basis[0] = one;

        for j in 1..=D {
            left[j] = t - self.knots[span + 1 - j];
            right[j] = self.knots[span + j] - t;
            let mut saved = zero;
            for r in 0..j {
                let denom = right[r + 1] + left[j - r];
                let temp = if denom.abs() <= P::Scalar::epsilon() {
                    zero
                } else {
                    basis[r] / denom
                };
                basis[r] = saved + right[r + 1] * temp;
                saved = left[j - r] * temp;
            }
            basis[j] = saved;
        }

        Ok(basis)
    }

    /// Evaluate the basis functions and their derivatives up to order `R`.
    ///
    /// Returns a `(R + 1) x (D + 1)` table `ders[k][j]` where `k` is the
    /// derivative order and `j` indexes the basis functions for the active span
    /// (`j = 0` corresponds to `i = span - D`). The zeroth row `ders[0]` matches
    /// `basis_functions(span, t)`. Requires `R <= D`.
    ///
    /// # Examples
    /// ```rust
    /// use stroke::{BSpline, PointN};
    ///
    /// let knots: [f64; 4] = [0.0, 0.0, 1.0, 1.0];
    /// let points = [PointN::new([0.0]), PointN::new([2.0])];
    /// let curve: BSpline<PointN<f64, 1>, 4, 2, 1> = BSpline::new(knots, points).unwrap();
    ///
    /// let t = 0.25;
    /// let span = curve.knot_span(t).unwrap();
    /// let ders = curve.basis_functions_with_derivatives::<1>(span, t).unwrap();
    /// assert!((ders[0][0] - 0.75).abs() < 1e-6);
    /// assert!((ders[0][1] - 0.25).abs() < 1e-6);
    /// assert!((ders[1][0] + 1.0).abs() < 1e-6);
    /// assert!((ders[1][1] - 1.0).abs() < 1e-6);
    /// ```
    pub fn basis_functions_with_derivatives<const R: usize>(
        &self,
        span: usize,
        t: P::Scalar,
    ) -> Result<[[P::Scalar; D + 1]; R + 1], BSplineError>
    where
        [(); D + 1]: Sized,
        [(); R + 1]: Sized,
        [(); D - R]: Sized,
    {
        if !self.domain_contains(t) {
            return Err(BSplineError::KnotDomainViolation);
        }
        if span < D || span + D >= K {
            return Err(BSplineError::KnotDomainViolation);
        }

        let zero = <P::Scalar as NumCast>::from(0.0).unwrap();
        let one = <P::Scalar as NumCast>::from(1.0).unwrap();
        let mut left = [zero; D + 1];
        let mut right = [zero; D + 1];
        let mut ndu = [[zero; D + 1]; D + 1];
        ndu[0][0] = one;

        for j in 1..=D {
            left[j] = t - self.knots[span + 1 - j];
            right[j] = self.knots[span + j] - t;
            let mut saved = zero;
            for r in 0..j {
                let denom = right[r + 1] + left[j - r];
                ndu[j][r] = denom;
                let temp = if denom.abs() <= P::Scalar::epsilon() {
                    zero
                } else {
                    ndu[r][j - 1] / denom
                };
                ndu[r][j] = saved + right[r + 1] * temp;
                saved = left[j - r] * temp;
            }
            ndu[j][j] = saved;
        }

        let mut ders = [[zero; D + 1]; R + 1];
        for j in 0..=D {
            ders[0][j] = ndu[j][D];
        }

        for r in 0..=D {
            let mut a = [[zero; D + 1]; 2];
            a[0][0] = one;
            let mut s1 = 0usize;
            let mut s2 = 1usize;

            for k in 1..=R {
                let mut d = zero;
                let rk = r as isize - k as isize;
                let pk = D - k;

                if r >= k {
                    let denom = ndu[pk + 1][rk as usize];
                    a[s2][0] = if denom.abs() <= P::Scalar::epsilon() {
                        zero
                    } else {
                        a[s1][0] / denom
                    };
                    d = a[s2][0] * ndu[rk as usize][pk];
                }

                let j1 = if rk >= -1 { 1 } else { (-rk) as usize };
                let j2 = if r <= pk { k - 1 } else { D - r };

                if j1 <= j2 {
                    for j in j1..=j2 {
                        let denom = ndu[pk + 1][(rk + j as isize) as usize];
                        a[s2][j] = if denom.abs() <= P::Scalar::epsilon() {
                            zero
                        } else {
                            (a[s1][j] - a[s1][j - 1]) / denom
                        };
                        d = d + a[s2][j] * ndu[(rk + j as isize) as usize][pk];
                    }
                }

                if r <= pk {
                    let denom = ndu[pk + 1][r];
                    a[s2][k] = if denom.abs() <= P::Scalar::epsilon() {
                        zero
                    } else {
                        -a[s1][k - 1] / denom
                    };
                    d = d + a[s2][k] * ndu[r][pk];
                }

                ders[k][r] = d;
                core::mem::swap(&mut s1, &mut s2);
            }
        }

        let mut factor = one;
        for k in 1..=R {
            let mult = <P::Scalar as NumCast>::from((D - k + 1) as f64).unwrap();
            factor = factor * mult;
            for j in 0..=D {
                ders[k][j] = ders[k][j] * factor;
            }
        }

        Ok(ders)
    }

    fn knot_domain_clamped_start(&self) -> (P::Scalar, P::Scalar) {
        // Start from knots[D], end at knots[knots.len() - 1]
        (self.knots[D], self.knots[self.knots.len() - 1])
    }

    fn knot_domain_clamped_end(&self) -> (P::Scalar, P::Scalar) {
        // Start from knots[0], end at knots[knots.len() - 1 - D] (C = K - D - 1)
        (self.knots[0], self.knots[self.knots.len() - 1 - D])
    }

    // Returns the knot domain for a clamped B-Spline
    // where the first and last knot has multiplicity D+1
    fn knot_domain_clamped(&self) -> (P::Scalar, P::Scalar) {
        // The valid domain is from knots[D] to knots[C] (C = K - D - 1)
        (self.knots[D], self.knots[self.knots.len() - 1 - D])
    }

    // Returns the knot domain for an unclamped B-Spline
    // where the first and last knot has multiplicity 1
    // fn knot_domain_unclamped(&self) -> (P::Scalar, P::Scalar) {
    //     (self.knots[0], self.knots[self.knots.len() - 1])
    // }
    fn knot_domain_unclamped(&self) -> (P::Scalar, P::Scalar) {
        // The valid domain is from knots[D] to knots[C] (C = K - D - 1)
        (self.knots[D], self.knots[C])
    }

    /// Return the bounding box of the curve as an array of (min, max) tuples for each dimension.
    pub fn bounding_box(&self) -> [(P::Scalar, P::Scalar); P::DIM]
    where
        P: PointIndex,
        [(); D + 1]: Sized,
        [(); D - 1]: Sized,
        [(); (D - 1) + 1]: Sized, // keep (D - 1) + 1; compiler doesn't normalize to D
        [(); C - 1]: Sized,
        [(); K - 2]: Sized,
    {
        let (kmin, kmax) = self.knot_domain();
        let tolerance = <P::Scalar as NumCast>::from(ROOT_TOLERANCE).unwrap();
        let half = <P::Scalar as NumCast>::from(0.5).unwrap();
        let zero = <P::Scalar as NumCast>::from(0.0).unwrap();

        let start = self.eval(kmin).unwrap_or(self.control_points[0]);
        let end = self.eval(kmax).unwrap_or(self.control_points[C - 1]);

        let derivative = self.derivative();
        let mut bounds = [(zero, zero); P::DIM];

        for dim in 0..P::DIM {
            let mut min = start[dim];
            let mut max = min;
            let end_value = end[dim];
            if end_value < min {
                min = end_value;
            }
            if end_value > max {
                max = end_value;
            }

            let mut update = |t: P::Scalar| {
                if let Ok(p) = self.eval(t) {
                    let value = p[dim];
                    if value < min {
                        min = value;
                    }
                    if value > max {
                        max = value;
                    }
                }
            };

            for i in 0..K - 1 {
                let mut t0 = self.knots[i];
                let mut t1 = self.knots[i + 1];
                if t1 <= t0 {
                    continue;
                }
                if t1 < kmin || t0 > kmax {
                    continue;
                }
                if t0 < kmin {
                    t0 = kmin;
                }
                if t1 > kmax {
                    t1 = kmax;
                }

                let steps = (D + 1) * 2;
                let steps_scalar = <P::Scalar as NumCast>::from(steps as f64).unwrap();
                let span = t1 - t0;
                let mut prev_t = t0;
                let mut prev_f = match derivative.eval(prev_t) {
                    Ok(p) => p[dim],
                    Err(_) => continue,
                };
                if prev_f.abs() <= tolerance {
                    update(prev_t);
                }

                for j in 1..=steps {
                    let alpha = <P::Scalar as NumCast>::from(j as f64).unwrap() / steps_scalar;
                    let t = t0 + span * alpha;
                    let f = match derivative.eval(t) {
                        Ok(p) => p[dim],
                        Err(_) => break,
                    };
                    if f.abs() <= tolerance {
                        update(t);
                    }
                    if prev_f * f < zero {
                        if let Some(root) =
                            self.bisect_axis_root(&derivative, dim, prev_t, t, tolerance, half)
                        {
                            update(root);
                        }
                    }
                    prev_t = t;
                    prev_f = f;
                }
            }

            bounds[dim] = (min, max);
        }

        bounds
    }

    fn bisect_axis_root(
        &self,
        derivative: &BSpline<P, { K - 2 }, { C - 1 }, { D - 1 }>,
        axis: usize,
        mut a: P::Scalar,
        mut b: P::Scalar,
        tolerance: P::Scalar,
        half: P::Scalar,
    ) -> Option<P::Scalar>
    where
        P: PointIndex,
        [(); D - 1]: Sized,
        [(); (D - 1) + 1]: Sized, // keep (D - 1) + 1; compiler doesn't normalize to D
        [(); C - 1]: Sized,
        [(); K - 2]: Sized,
    {
        let zero = <P::Scalar as NumCast>::from(0.0).unwrap();
        let mut fa = derivative.eval(a).ok().map(|p| p[axis])?;
        let fb = derivative.eval(b).ok().map(|p| p[axis])?;
        if fa.abs() <= tolerance {
            return Some(a);
        }
        if fb.abs() <= tolerance {
            return Some(b);
        }
        if fa * fb > zero {
            return None;
        }

        for _ in 0..MAX_ROOT_ITER {
            let mid = (a + b) * half;
            let fm = derivative.eval(mid).ok().map(|p| p[axis])?;
            if fm.abs() <= tolerance || (b - a).abs() <= tolerance {
                return Some(mid);
            }
            if fa * fm <= zero {
                b = mid;
            } else {
                a = mid;
                fa = fm;
            }
        }

        Some((a + b) * half)
    }

    /// Approximate the minimum distance between given `point` and the curve.
    /// Uses a coarse sampling pass over the knot domain and a local search around
    /// the best sample. `nsteps` is the number of coarse samples.
    pub fn distance_to_point_approx(
        &self,
        point: P,
        nsteps: usize,
    ) -> Result<P::Scalar, BSplineError>
    where
        P: PointNorm,
        [(); D + 1]: Sized,
    {
        let nsteps = nsteps.max(1);
        let (kmin, kmax) = self.knot_domain();
        let span = kmax - kmin;
        let nsteps_scalar = <P::Scalar as NumCast>::from(nsteps as f64).unwrap();
        let half = <P::Scalar as NumCast>::from(0.5).unwrap();
        let three = <P::Scalar as NumCast>::from(3.0).unwrap();

        let mut best_i = 0usize;
        let mut best_d = (self.eval(kmin)? - point).squared_norm();
        for i in 1..=nsteps {
            let i_scalar = <P::Scalar as NumCast>::from(i as f64).unwrap();
            let t = kmin + span * (i_scalar / nsteps_scalar);
            let candidate = self.eval(t)?;
            let distance = (candidate - point).squared_norm();
            if distance < best_d {
                best_d = distance;
                best_i = i;
            }
        }

        let left_i = if best_i == 0 { 0 } else { best_i - 1 };
        let right_i = if best_i == nsteps { nsteps } else { best_i + 1 };
        let mut left =
            kmin + span * (<P::Scalar as NumCast>::from(left_i as f64).unwrap() / nsteps_scalar);
        let mut right =
            kmin + span * (<P::Scalar as NumCast>::from(right_i as f64).unwrap() / nsteps_scalar);

        for _ in 0..LOCAL_SEARCH_ITERS {
            let third = (right - left) / three;
            let t1 = left + third;
            let t2 = right - third;
            let d1 = (self.eval(t1)? - point).squared_norm();
            let d2 = (self.eval(t2)? - point).squared_norm();
            if d1 < d2 {
                right = t2;
            } else {
                left = t1;
            }
        }

        let t = (left + right) * half;
        Ok((self.eval(t)? - point).squared_norm().sqrt())
    }

    /// Approximate the minimum distance between given `point` and the curve using
    /// a default sampling resolution.
    pub fn distance_to_point(&self, point: P) -> Result<P::Scalar, BSplineError>
    where
        P: PointNorm,
        [(); D + 1]: Sized,
    {
        self.distance_to_point_approx(point, DEFAULT_DISTANCE_STEPS)
    }

    /// Iteratively compute de Boor's B-spline algorithm, this computes the recursive
    /// de Boor algorithm tree from the bottom up. At each level we use the results
    /// from the previous one to compute this level and store the results in the
    /// array indices we no longer need to compute the current level (the left one
    /// used computing node j).
    fn de_boor_iterative(&self, t: P::Scalar) -> Result<P, BSplineError>
    where
        [(); D + 1]: Sized,
    {
        // Special case when t equals the last knot value
        if t == self.knots[self.control_points.len()] {
            return Ok(self.control_points[self.control_points.len() - 1]);
        }

        let k = self.knot_span_start_for_t(t);

        let k = match k {
            Some(r) => r,
            None => return Err(BSplineError::KnotDomainViolation),
        };

        // Create a second array for storing intermediate results
        let mut d: [P; D + 1] = core::array::from_fn(|_| self.control_points[0]);

        // Copy the control points needed for the first iteration
        // For a B-Spline of degree D, D+1 control points in the knot span contribute to the result
        for j in 0..=D {
            d[j] = self.control_points[j + k - D];
        }

        // // Adjust the start index to prevent negative values
        // let start_index = if k >= D { k - D } else { 0 };
        // // Calculate the number of control points to copy
        // let num_points = if k >= D { D + 1 } else { k + 1 };
        // // Copy the control points needed for the first iteration
        // for j in 0..num_points {
        //     d[j] = self.control_points[start_index + j];
        // }

        // de Boor's algorithm
        // The result overwrites the right of the two point operands used in the current iteration
        // The first iteration computes the result for the first knot span,
        // the second iteration for the next knot span, and so on.
        for r in 1..=D {
            for j in (r..=D).rev() {
                // calculate alpha, avoid division by zero
                let numerator = t - self.knots[j + k - D];
                let denominator = self.knots[j + 1 + k - r] - self.knots[j + k - D];
                let zero = <P::Scalar as NumCast>::from(0.0).unwrap();
                let one = <P::Scalar as NumCast>::from(1.0).unwrap();
                let alpha: P::Scalar;
                if numerator < P::Scalar::epsilon() {
                    // we are at the start of a knot span (interpolation d[j-1] -> d[j])
                    alpha = zero;
                } else if denominator.abs() < P::Scalar::epsilon() {
                    // we are at then end of a knot span (interpolation d[j-1] -> d[j])
                    // this will lead to d[j] simply being copied to d[j]
                    alpha = one;
                } else {
                    // inside a knot span (interpolation d[j-1] -> d[j])
                    alpha = numerator / denominator;
                }
                if alpha.is_nan() {
                    return Err(BSplineError::EvaluationAlphaIsNaN);
                }
                // The right of two points is overwritten by the result
                // to be used in the next iteration
                d[j] = d[j - 1] * (one - alpha) + d[j] * alpha;
            }
        }
        Ok(d[D])
    }

    fn knot_span_start_for_t(&self, t: P::Scalar) -> Option<usize> {
        match self.knot_kind {
            KnotVectorKind::Unclamped => self.knot_span_start_for_t_unclamped(t),
            KnotVectorKind::ClampedStart => self.knot_span_start_for_t_clamped_start(t),
            KnotVectorKind::Clamped => self.knot_span_start_for_t_clamped(t),
            KnotVectorKind::ClampedEnd => self.knot_span_start_for_t_clamped_end(t),
        }
    }

    /// Return the index of the start of the knot span for given `t` for an unclamped B-Spline.
    /// This is the index of the first element that satisfies:
    ///     knots[i] <= t < knots[i+1]
    /// Because the knot vector is non-decreasing, this function uses binary search.
    /// If no element greater than the value passed is found, the function returns None.
    fn knot_span_start_for_t_unclamped(&self, t: P::Scalar) -> Option<usize> {
        if t == self.knots[C] {
            return Some(C - 1);
        }

        let mut low = D;
        let mut high = C - 1;
        while low <= high {
            let mid = (low + high) / 2;
            if t >= self.knots[mid] && t < self.knots[mid + 1] {
                return Some(mid);
            }
            if t < self.knots[mid] {
                if mid == 0 {
                    return None;
                }
                high = mid - 1;
            } else {
                low = mid + 1;
            }
        }
        None
    }

    /// Return the index of the start of the knot span for given `t` for a clamped B-Spline.
    /// This is the index of the first element that satisfies:
    ///     knots[i] <= t < knots[i+1] where i >= degree and i < knots.len() - degree - 1
    /// Because the knot vector is non-decreasing, this function uses binary search.
    /// If no element greater than the value passed is found, the function returns None.
    fn knot_span_start_for_t_clamped(&self, t: P::Scalar) -> Option<usize> {
        let n = self.control_points.len() - 1; // n = number of control points - 1
        let p = D; // degree
        if t == self.knots[n + 1] {
            return Some(n);
        }
        let mut low = p;
        let mut high = n;
        while low <= high {
            let mid = (low + high) / 2;
            if t >= self.knots[mid] && t < self.knots[mid + 1] {
                return Some(mid);
            } else if t < self.knots[mid] {
                high = mid - 1;
            } else {
                low = mid + 1;
            }
        }
        None
    }

    fn knot_span_start_for_t_clamped_start(&self, t: P::Scalar) -> Option<usize> {
        let n = self.control_points.len() - 1; // n = number of control points - 1
        let p = D; // Degree of the spline
        let m = self.knots.len() - 1; // Last index of the knot vector

        // Special case when t equals the last knot value
        if t == self.knots[m] {
            return Some(n);
        }

        // Valid parameter domain is from knots[p] to knots[m]
        if t < self.knots[p] || t > self.knots[m] {
            return None;
        }

        // Initialize binary search indices
        let mut low = p;
        let mut high = m - 1; // Since we compare t with knots[mid + 1], high should be m - 1

        // Binary search to find the knot span
        while low <= high {
            let mid = (low + high) / 2;
            if t >= self.knots[mid] && t < self.knots[mid + 1] {
                return Some(mid);
            } else if t < self.knots[mid] {
                high = mid - 1;
            } else {
                low = mid + 1;
            }
        }

        None // t is not within any knot span
    }

    fn knot_span_start_for_t_clamped_end(&self, t: P::Scalar) -> Option<usize> {
        let n = self.control_points.len() - 1; // n = number of control points - 1

        // Special case when t equals knots[n + 1]
        if t == self.knots[n + 1] {
            return Some(n);
        }

        // Valid parameter domain is from knots[0] to knots[n + 1]
        if t < self.knots[0] || t > self.knots[n + 1] {
            return None;
        }

        // Initialize binary search indices
        let mut low = 0;
        let mut high = n; // Since we compare t with knots[mid + 1], high should be n

        // Binary search to find the knot span
        while low <= high {
            let mid = (low + high) / 2;
            if t >= self.knots[mid] && t < self.knots[mid + 1] {
                return Some(mid);
            } else if t < self.knots[mid] {
                high = mid - 1;
            } else {
                low = mid + 1;
            }
        }

        None // t is not within any knot span
    }

    /// Approximates the arc length of the curve by flattening it with straight line segments.
    /// This approximation is unfeasable if desired accuracy is greater than ~2 decimal places.
    /// Returns an error if evaluation fails inside the knot domain.
    pub fn arclen(&self, nsteps: usize) -> Result<P::Scalar, BSplineError>
    where
        P: PointNorm,
        [(); D + 1]: Sized,
    {
        let nsteps = nsteps.max(1);
        let nsteps_scalar = <P::Scalar as NumCast>::from(nsteps as f64).unwrap();
        let zero = <P::Scalar as NumCast>::from(0.0).unwrap();
        let mut arclen: P::Scalar = zero;
        let (kmin, kmax) = self.knot_domain();
        let span = kmax - kmin;
        let dt = span / nsteps_scalar;

        for i in 0..nsteps {
            let i_scalar = <P::Scalar as NumCast>::from(i as f64).unwrap();
            let mut t0 = kmin + dt * i_scalar;
            let mut t1 = if i + 1 == nsteps { kmax } else { t0 + dt };
            if t0 < kmin {
                t0 = kmin;
            } else if t0 > kmax {
                t0 = kmax;
            }
            if t1 < kmin {
                t1 = kmin;
            } else if t1 > kmax {
                t1 = kmax;
            }

            let p1 = self.eval(t0)?;
            let p2 = self.eval(t1)?;
            arclen = arclen + (p1 - p2).squared_norm().sqrt();
        }

        Ok(arclen)
    }

    fn arclen_partial(&self, t: P::Scalar, nsteps: usize) -> Result<P::Scalar, BSplineError>
    where
        P: PointNorm,
        [(); D + 1]: Sized,
    {
        let zero = <P::Scalar as NumCast>::from(0.0).unwrap();
        let nsteps = nsteps.max(1);
        let (kmin, kmax) = self.knot_domain();
        let t = t.clamp(kmin, kmax);
        if t <= kmin {
            return Ok(zero);
        }

        let nsteps_scalar = <P::Scalar as NumCast>::from(nsteps as f64).unwrap();
        let span = t - kmin;
        let dt = span / nsteps_scalar;
        let mut arclen: P::Scalar = zero;
        let mut prev = self.eval(kmin)?;

        for i in 1..=nsteps {
            let i_scalar = <P::Scalar as NumCast>::from(i as f64).unwrap();
            let ti = if i == nsteps { t } else { kmin + dt * i_scalar };
            let p = self.eval(ti)?;
            arclen = arclen + (p - prev).squared_norm().sqrt();
            prev = p;
        }

        Ok(arclen)
    }

    /// Approximate parameter `t` at arc length `s`.
    pub fn t_at_length_approx(
        &self,
        s: P::Scalar,
        nsteps: usize,
    ) -> Result<P::Scalar, BSplineError>
    where
        P: PointNorm,
        [(); D + 1]: Sized,
    {
        let zero = <P::Scalar as NumCast>::from(0.0).unwrap();
        let half = <P::Scalar as NumCast>::from(0.5).unwrap();
        let (kmin, kmax) = self.knot_domain();

        let total = self.arclen(nsteps)?;
        if total <= P::Scalar::epsilon() {
            return Ok(kmin);
        }
        let target = s.clamp(zero, total);

        let mut lo = kmin;
        let mut hi = kmax;
        for _ in 0..MAX_LENGTH_ITERS {
            let mid = (lo + hi) * half;
            let len = self.arclen_partial(mid, nsteps)?;
            if len < target {
                lo = mid;
            } else {
                hi = mid;
            }
        }

        Ok((lo + hi) * half)
    }

    /// Approximate parameter `t` at arc length `s` using a default resolution.
    pub fn t_at_length(&self, s: P::Scalar) -> Result<P::Scalar, BSplineError>
    where
        P: PointNorm,
        [(); D + 1]: Sized,
    {
        self.t_at_length_approx(s, DEFAULT_LENGTH_STEPS)
    }

    /// Evaluate the point at arc length `s`.
    pub fn point_at_length_approx(&self, s: P::Scalar, nsteps: usize) -> Result<P, BSplineError>
    where
        P: PointNorm,
        [(); D + 1]: Sized,
    {
        let t = self.t_at_length_approx(s, nsteps)?;
        self.eval(t)
    }

    /// Evaluate the point at arc length `s` using a default resolution.
    pub fn point_at_length(&self, s: P::Scalar) -> Result<P, BSplineError>
    where
        P: PointNorm,
        [(); D + 1]: Sized,
    {
        self.point_at_length_approx(s, DEFAULT_LENGTH_STEPS)
    }

    /// Returns the derivative curve of self which has C-1 control points and K-2 knots.
    /// The derivative of an nth degree B-Spline curve is an (n-1)th degree (d) B-Spline curve,
    /// with the same knot vector, and new control points Q0...Qn-1 derived from the
    /// original control points Pi as:
    ///                 d
    /// Qi =    ----------------- (P[i+1]-P[i])
    ///         k[i+d+1] - k[i+1].
    /// with degree = curve_order - 1
    pub fn derivative(&self) -> BSpline<P, { K - 2 }, { C - 1 }, { D - 1 }>
    where
        [(); D - 1]: Sized,
        [(); C - 1]: Sized,
        [(); K - 2]: Sized,
    {
        let derivative_knots: [P::Scalar; K - 2] = core::array::from_fn(|i| self.knots[i + 1]);
        let derivative_points: [P; C - 1] = core::array::from_fn(|i| {
            let zero = <P::Scalar as NumCast>::from(0.0).unwrap();
            let denom = self.knots[i + D + 1] - self.knots[i + 1];
            let scale = if denom.abs() < P::Scalar::epsilon() {
                zero
            } else {
                <P::Scalar as NumCast>::from(D as f64).unwrap() / denom
            };
            (self.control_points[i + 1] - self.control_points[i]) * scale
        });

        let degree = D - 1;
        let is_clamped_front = derivative_knots[0..=degree]
            .iter()
            .all(|&knot| (knot - derivative_knots[0]).abs() < P::Scalar::epsilon());
        let is_clamped_back = derivative_knots[derivative_knots.len() - degree - 1..]
            .iter()
            .all(|&knot| {
                (knot - derivative_knots[derivative_knots.len() - 1]).abs() < P::Scalar::epsilon()
            });

        let knot_kind = match (is_clamped_front, is_clamped_back) {
            (false, false) => KnotVectorKind::Unclamped,
            (true, false) => KnotVectorKind::ClampedStart,
            (false, true) => KnotVectorKind::ClampedEnd,
            (true, true) => KnotVectorKind::Clamped,
        };

        BSpline {
            knots: derivative_knots,
            knot_kind,
            control_points: derivative_points,
        }
    }

    /// Return the unit tangent direction at `t`.
    pub fn tangent(&self, t: P::Scalar) -> Result<P, BSplineError>
    where
        P: PointNorm,
        [(); D - 1]: Sized,
        [(); C - 1]: Sized,
        [(); K - 2]: Sized,
        [(); D]: Sized,
        [(); D + 1]: Sized,
        [(); (D - 1) + 1]: Sized, // keep (D - 1) + 1; compiler doesn't normalize to D
    {
        let one = <P::Scalar as NumCast>::from(1.0).unwrap();
        let dir = self.derivative().eval(t)?;
        let len = dir.squared_norm().sqrt();
        if len <= P::Scalar::epsilon() {
            Ok(dir)
        } else {
            Ok(dir * (one / len))
        }
    }

    /// Return the curvature magnitude at `t`.
    ///
    /// Requires `D >= 2` so that the second derivative exists.
    pub fn curvature(&self, t: P::Scalar) -> Result<P::Scalar, BSplineError>
    where
        P: PointNorm + PointDot,
        [(); D + 1]: Sized,
        [(); D - 2]: Sized,
    {
        let span = self.knot_span(t)?;
        let ders = self.basis_functions_with_derivatives::<2>(span, t)?;
        let zero = <P::Scalar as NumCast>::from(0.0).unwrap();
        let mut v = self.control_points[0] * zero;
        let mut a = v;
        let start = span - D;
        for j in 0..=D {
            let cp = self.control_points[start + j];
            v = v + cp * ders[1][j];
            a = a + cp * ders[2][j];
        }

        let v2 = v.squared_norm();
        if v2 <= P::Scalar::epsilon() {
            return Ok(zero);
        }
        let a2 = a.squared_norm();
        let dot = PointDot::dot(&v, &a);
        let mut num = v2 * a2 - dot * dot;
        if num < zero {
            num = zero;
        }
        let denom = v2 * v2.sqrt();
        if denom <= P::Scalar::epsilon() {
            Ok(zero)
        } else {
            Ok(num.sqrt() / denom)
        }
    }

    /// Return the principal normal direction at `t`.
    ///
    /// Requires `D >= 2` so that the second derivative exists. Returns `None`
    /// if the velocity is zero or curvature is undefined.
    pub fn normal(&self, t: P::Scalar) -> Result<Option<P>, BSplineError>
    where
        P: PointNorm + PointDot,
        [(); D + 1]: Sized,
        [(); D - 2]: Sized,
    {
        let span = self.knot_span(t)?;
        let ders = self.basis_functions_with_derivatives::<2>(span, t)?;
        let zero = <P::Scalar as NumCast>::from(0.0).unwrap();
        let mut v = self.control_points[0] * zero;
        let mut a = v;
        let start = span - D;
        for j in 0..=D {
            let cp = self.control_points[start + j];
            v = v + cp * ders[1][j];
            a = a + cp * ders[2][j];
        }

        let v2 = v.squared_norm();
        if v2 <= P::Scalar::epsilon() {
            return Ok(None);
        }
        let dot = PointDot::dot(&v, &a);
        let a_perp = a - v * (dot / v2);
        let n2 = a_perp.squared_norm();
        if n2 <= P::Scalar::epsilon() {
            return Ok(None);
        }
        let one = <P::Scalar as NumCast>::from(1.0).unwrap();
        Ok(Some(a_perp * (one / n2.sqrt())))
    }

    fn derivative_curve(
        &self,
    ) -> Result<BSpline<P, { K - 2 }, { C - 1 }, { D - 1 }>, RootFindingError>
    where
        [(); D - 1]: Sized,
        [(); C - 1]: Sized,
        [(); K - 2]: Sized,
    {
        Ok(self.derivative())
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
        [(); D + 1]: Sized,
        [(); D - 1]: Sized,
        [(); D]: Sized,
        [(); C - 1]: Sized,
        [(); (D - 1) + 1]: Sized, // keep (D - 1) + 1; compiler doesn't normalize to D
        [(); K - 2]: Sized,
        P: PointIndex,
    {
        FindRoot::root_newton_axis(self, value, axis, start, eps, max_iter)
    }
}

impl<P, const K: usize, const C: usize, const D: usize> FindRoot<P> for BSpline<P, K, C, D>
where
    P: PointIndex,
    [(); D + 1]: Sized,
    [(); D - 1]: Sized,
    [(); D]: Sized,
    [(); C - 1]: Sized,
    [(); (D - 1) + 1]: Sized, // keep (D - 1) + 1; compiler doesn't normalize to D
    [(); K - 2]: Sized,
{
    fn parameter_domain(&self) -> (P::Scalar, P::Scalar) {
        self.knot_domain()
    }

    fn axis_value(&self, t: P::Scalar, axis: usize) -> Result<P::Scalar, RootFindingError> {
        self.eval(t)
            .map(|p| p[axis])
            .map_err(|_| RootFindingError::FailedToConverge)
    }

    fn axis_derivative(&self, t: P::Scalar, axis: usize) -> Result<P::Scalar, RootFindingError> {
        let derivative = self.derivative_curve()?;
        derivative
            .eval(t)
            .map(|p| p[axis])
            .map_err(|_| RootFindingError::FailedToConverge)
    }

    fn root_newton_axis(
        &self,
        value: P::Scalar,
        axis: usize,
        start: P::Scalar,
        eps: Option<P::Scalar>,
        max_iter: Option<usize>,
    ) -> Result<P::Scalar, RootFindingError> {
        let eps = eps.unwrap_or_else(|| <P::Scalar as NumCast>::from(1e-6).unwrap());
        let max_iter = max_iter.unwrap_or(64);
        let (kmin, kmax) = self.knot_domain();
        let derivative = self.derivative_curve()?;

        let clamp_t = |mut t: P::Scalar| {
            if t < kmin {
                t = kmin;
            } else if t > kmax {
                t = kmax;
            }
            t
        };

        let start = clamp_t(start);
        let fx = |x: P::Scalar| {
            let t = clamp_t(x);
            self.axis_value(t, axis).map(|v| v - value)
        };
        let dx = |x: P::Scalar| {
            let t = clamp_t(x);
            derivative
                .eval(t)
                .map(|p| p[axis])
                .map_err(|_| RootFindingError::FailedToConverge)
        };

        let root = root_newton_raphson(start, fx, dx, eps, max_iter)?;
        Ok(clamp_t(root))
    }
}

#[cfg(test)]
mod tests {
    //use std;
    use super::*;
    //use crate::num_traits::{Pow};
    use crate::{EPSILON, PointN};

    #[test]
    fn degree_1_clamped_construct_and_eval_endpoints() {
        // Define the control points for a degree 1 B-Spline (essentially a line)
        let points = [
            PointN::new([0f64, 0f64]),   // Start of the line at origin
            PointN::new([10f64, 10f64]), // End of the line at (10, 10)
        ];

        // Define a uniform clamped knot vector for a degree 1 B-spline
        let knots: [f64; 4] = [0., 0., 1., 1.];

        // Create the B-Spline
        let curve: BSpline<PointN<f64, 2>, 4, 2, 1> = BSpline::new(knots, points).unwrap();

        assert!(curve.knot_kind == KnotVectorKind::Clamped);

        // evaluate along the curve to find out of bounds errors
        let (kmin, kmax) = curve.knot_domain();
        let nsteps: usize = 100;
        for t in 0..=nsteps {
            let t = kmin + (t as f64 / nsteps as f64) * (kmax - kmin);
            //dbg!(kmin, t, kmax);
            // this curve shouldn't fail to evaluate numerically unless knot index is out of bounds
            let _ = curve.eval(t).unwrap_or_else(|e| {
                panic!("Unexpected error in eval: {:?}", e);
            });
        }

        // Evaluate the B-Spline at the start (t=0) and end (t=1) of the knot vector
        let start = curve.eval(0.).unwrap();
        let end = curve.eval(1.).unwrap();

        // Check that the start and end points match the control points
        assert_eq!(start, points[0]);
        assert_eq!(end, points[1]);

        let mid = curve.eval(0.5).unwrap();
        assert!((mid - PointN::new([5.0, 5.0])).squared_norm() < EPSILON);
    }

    #[test]
    fn degree_1_unclamped_construct_and_eval() {
        // Define the control points for a degree 1 B-Spline (essentially a line)
        let points = [
            PointN::new([0f64, 0f64]),   // Start of the line at origin
            PointN::new([10f64, 10f64]), // End of the line at (10, 10)
        ];

        // Define a uniform unclamped knot vector for a degree 1 B-spline
        let knots: [f64; 4] = [0., 1., 2., 3.];

        // Create the B-Spline
        let curve: BSpline<PointN<f64, 2>, 4, 2, 1> = BSpline::new(knots, points).unwrap();

        assert!(curve.knot_kind == KnotVectorKind::Unclamped);

        let (kmin, kmax) = curve.knot_domain();

        // evaluate along the curve to find out of bounds errors
        let nsteps: usize = 100;
        for t in 0..=nsteps {
            let t = kmin + (t as f64 / nsteps as f64) * (kmax - kmin);
            //dbg!(kmin, t, kmax);
            // this curve shouldn't fail to evaluate numerically unless knot index is out of bounds
            let _ = curve.eval(t).unwrap_or_else(|e| {
                panic!("Unexpected error in eval: {:?}", e);
            });
        }

        // Evaluate at the midpoint of the knot domain.
        let mid = curve.eval(1.5).unwrap();
        assert!((mid - PointN::new([5.0, 5.0])).squared_norm() < EPSILON);
    }

    #[test]
    fn arclen_line_approx() {
        let points = [PointN::new([0f64, 0f64]), PointN::new([10f64, 0f64])];
        let knots: [f64; 4] = [0.0, 0.0, 1.0, 1.0];
        let curve: BSpline<PointN<f64, 2>, 4, 2, 1> = BSpline::new(knots, points).unwrap();

        let expected = (points[1] - points[0]).squared_norm().sqrt();
        let length = curve.arclen(32).unwrap();
        let tolerance = 1e-6;
        assert!((length - expected).abs() <= tolerance);
    }

    #[test]
    fn degree_2_clamped_construct_and_eval_endpoints() {
        // Define the control points for a degree 2 B-Spline
        let points = [
            PointN::new([0f64, 0f64]),   // Start of the curve at origin
            PointN::new([10f64, 10f64]), // Control point off the baseline
            PointN::new([20f64, 0f64]),  // End of the curve at (20, 0)
        ];

        // Define a uniform clamped knot vector for a degree 2 B-spline
        let knots: [f64; 6] = [0., 0., 0., 1., 1., 1.];

        // Create the B-Spline
        let curve: BSpline<PointN<f64, 2>, 6, 3, 2> = BSpline::new(knots, points).unwrap();

        assert!(curve.knot_kind == KnotVectorKind::Clamped);

        let (kmin, kmax) = curve.knot_domain();

        // evaluate along the curve to find out of bounds errors
        let nsteps: usize = 100;
        for t in 0..=nsteps {
            let t = kmin + (t as f64 / nsteps as f64) * (kmax - kmin);
            //dbg!(kmin, t, kmax);
            // this curve shouldn't fail to evaluate numerically unless knot index is out of bounds
            let _ = curve.eval(t).unwrap_or_else(|e| {
                panic!("Unexpected error in eval: {:?}", e);
            });
        }

        // Evaluate the B-Spline at the start (t=0) and end (t=1) of the knot vector
        let start = curve.eval(0.).unwrap();
        let end = curve.eval(1.).unwrap();

        // Check that the start and end points match the control points
        assert_eq!(start, points[0]);
        assert_eq!(end, points[2]);

        let mid = curve.eval(0.5).unwrap();
        assert!((mid - PointN::new([10.0, 5.0])).squared_norm() < EPSILON);
    }

    #[test]
    fn degree_2_unclamped_construct_and_eval() {
        // Define the control points for a degree 2 B-Spline
        let points = [
            PointN::new([0f64, 0f64]),   // Start of the line at origin
            PointN::new([10f64, 10f64]), // End of the line at (10, 10)
            PointN::new([20f64, 0f64]),  // End of the line at (20, 0)
        ];

        // Define a uniform unclamped knot vector for a degree 2 B-spline
        let knots: [f64; 6] = [0., 1., 2., 3., 4., 5.];

        // Create the B-Spline
        let curve: BSpline<PointN<f64, 2>, 6, 3, 2> = BSpline::new(knots, points).unwrap();

        assert!(curve.knot_kind == KnotVectorKind::Unclamped);

        let (kmin, kmax) = curve.knot_domain();
        // evaluate along the curve to find out of bounds errors
        let nsteps: usize = 100;
        for t in 0..=nsteps {
            let t = kmin + (t as f64 / nsteps as f64) * (kmax - kmin);
            //dbg!(kmin, t, kmax);
            // this curve shouldn't fail to evaluate numerically unless knot index is out of bounds
            let _ = curve.eval(t).unwrap_or_else(|e| {
                panic!("Unexpected error in eval: {:?}", e);
            });
        }

        // Evaluate at the midpoint of the knot domain (t=2.5).
        let mid = curve.eval(2.5).unwrap();
        // For uniform unclamped degree-2 with 3 control points, weights are
        // [0.125, 0.75, 0.125] at t=2.5.
        let expected = points[0] * 0.125 + points[1] * 0.75 + points[2] * 0.125;
        assert!((mid - expected).squared_norm() < EPSILON);
    }

    #[test]
    fn degree_3_clamped_construct_and_eval_endpoints() {
        // degree 3, 4 control points => 4+3+1=8 knots
        let points = [
            PointN::new([0f64, 1.77f64]),
            PointN::new([1.1f64, -1f64]),
            PointN::new([4.3f64, 3f64]),
            PointN::new([3.2f64, -4f64]),
            PointN::new([2.3f64, 1.77f64]),
        ];
        // clamped b-spline to it goes through the first and last control point
        let knots: [f64; 9] = [0., 0., 0., 0., 1., 2., 2., 2., 2.];

        // try to make a b-spline with the given parameters
        let b: Result<BSpline<PointN<f64, 2>, 9, 5, 3>, BSplineError> = BSpline::new(knots, points);

        let curve = match b {
            Err(e) => panic!("Error: {:?}", e),
            Ok(b) => b,
        };

        assert!(curve.knot_kind == KnotVectorKind::Clamped);
        // evaluate the curve, t needs to be inside the knot domain!
        // we need to map [0...1] to kmin..kmax
        let (kmin, kmax) = curve.knot_domain();
        let nsteps: usize = 100;
        for t in 0..=nsteps {
            let t = kmin + (t as f64 / nsteps as f64) * (kmax - kmin);
            //dbg!(kmin, t, kmax);
            // this curve shouldn't fail to evaluate numerically unless knot index is out of bounds
            let _ = curve.eval(t).unwrap_or_else(|e| {
                panic!("Unexpected error in eval: {:?}", e);
            });
        }

        // check if eval returns the correct error if t is out of bounds
        assert_eq!(
            curve.eval(-0.1).err(),
            Some(BSplineError::KnotDomainViolation)
        );
        assert_eq!(
            curve.eval(2.1).err(),
            Some(BSplineError::KnotDomainViolation)
        );
        // check if eval works correctly for endpoints
        assert!(curve.eval(kmin).is_ok());
        assert!(curve.eval(kmax).is_ok());
        assert_eq!(curve.eval(kmin).unwrap(), points[0]);
        assert_eq!(curve.eval(kmax).unwrap(), points[4]);
    }

    #[test]
    fn degree_3_unclamped_construct_and_eval() {
        // degree 3, 4 control points => 4+3+1=8 knots
        let points = [
            PointN::new([0f64, 1.77f64]),
            PointN::new([1.1f64, -1f64]),
            PointN::new([4.3f64, 3f64]),
            PointN::new([3.2f64, -4f64]),
            PointN::new([2.3f64, 1.77f64]),
        ];
        // unclamped b-spline (not necessarily going through the first and last control point)
        let knots: [f64; 9] = [0., 1., 2., 3., 4., 5., 6., 7., 8.];

        // try to make a b-spline with the given parameters
        let b: Result<BSpline<PointN<f64, 2>, 9, 5, 3>, BSplineError> = BSpline::new(knots, points);

        let curve = match b {
            Err(e) => panic!("Error: {:?}", e),
            Ok(b) => b,
        };

        assert!(curve.knot_kind == KnotVectorKind::Unclamped);
        // evaluate the curve, t needs to be inside the knot domain!
        // we need to map [0...1] to kmin..kmax
        let (kmin, kmax) = curve.knot_domain();
        assert!(kmin == 3.0 && kmax == 5.0);
        let nsteps: usize = 100;
        for t in 0..=nsteps {
            let t = kmin + (t as f64 / nsteps as f64) * (kmax - kmin);
            //dbg!(kmin, t, kmax);
            // this curve shouldn't fail to evaluate numerically unless knot index is out of bounds
            let _ = curve.eval(t).unwrap_or_else(|e| {
                panic!("Unexpected error in eval: {:?}", e);
            });
        }

        // check if eval returns the correct error if t is out of bounds
        assert_eq!(
            curve.eval(-0.1).err(),
            Some(BSplineError::KnotDomainViolation)
        );
        assert_eq!(
            curve.eval(8.1).err(),
            Some(BSplineError::KnotDomainViolation)
        );
        // check if eval works correctly for endpoints
        // since the curve is unclampend, the curves endpoints are not necessarily on the control points
        assert!(curve.eval(kmin).is_ok());
        assert!(curve.eval(kmax).is_ok());
    }

    #[test]
    fn degree_4_clamped_construct_and_eval() {
        // degree 4, 5 control points => 5+4+1=10 knots
        let points = [
            PointN::new([0f64, 1.77f64]),
            PointN::new([1.1f64, -1f64]),
            PointN::new([4.3f64, 3f64]),
            PointN::new([3.2f64, -4f64]),
            PointN::new([2.3f64, 1.77f64]),
        ];
        // clamped b-spline (goes through the first and last control point)
        let knots: [f64; 10] = [0., 0., 0., 0., 0., 3., 3., 3., 3., 3.];

        // try to make a b-spline with the given parameters
        let b: Result<BSpline<PointN<f64, 2>, 10, 5, 4>, BSplineError> =
            BSpline::new(knots, points);

        let curve = match b {
            Err(e) => panic!("Error: {:?}", e),
            Ok(b) => b,
        };

        assert!(curve.knot_kind == KnotVectorKind::Clamped);
        // evaluate the curve, t needs to be inside the knot domain!
        // we need to map [0...1] to kmin..kmax
        let (kmin, kmax) = curve.knot_domain();
        assert!(kmin == 0.0 && kmax == 3.0);
        let nsteps: usize = 100;
        for t in 0..=nsteps {
            let t = kmin + (t as f64 / nsteps as f64) * (kmax - kmin);
            //dbg!(kmin, t, kmax);
            // this curve shouldn't fail to evaluate numerically unless knot index is out of bounds
            let _ = curve.eval(t).unwrap_or_else(|e| {
                panic!("Unexpected error in eval: {:?}", e);
            });
        }

        // check if eval returns the correct error if t is out of bounds
        assert_eq!(
            curve.eval(-0.1).err(),
            Some(BSplineError::KnotDomainViolation)
        );
        assert_eq!(
            curve.eval(3.1).err(),
            Some(BSplineError::KnotDomainViolation)
        );
        // check if eval works correctly for endpoints
        assert!(curve.eval(kmin).is_ok());
        assert!(curve.eval(kmax).is_ok());
        assert_eq!(curve.eval(kmin).unwrap(), points[0]);
        assert_eq!(curve.eval(kmax).unwrap(), points[4]);
    }

    #[test]
    fn degree_4_unclamped_construct_and_eval() {
        // degree 4, 5 control points => 5+4+1=10 knots
        let points = [
            PointN::new([0f64, 1.77f64]),
            PointN::new([1.1f64, -1f64]),
            PointN::new([4.3f64, 3f64]),
            PointN::new([3.2f64, -4f64]),
            PointN::new([2.3f64, 1.77f64]),
            PointN::new([0.0, 0.0]),
        ];
        // unclamped b-spline (not necessarily going through the first and last control point)
        let knots: [f64; 11] = [0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10.];

        // try to make a b-spline with the given parameters
        let b: Result<BSpline<PointN<f64, 2>, 11, 6, 4>, BSplineError> =
            BSpline::new(knots, points);

        let curve = match b {
            Err(e) => panic!("Error: {:?}", e),
            Ok(b) => b,
        };

        assert!(curve.knot_kind == KnotVectorKind::Unclamped);
        // evaluate the curve, t needs to be inside the knot domain!
        // we need to map [0...1] to kmin..kmax
        let (kmin, kmax) = curve.knot_domain();
        assert!(kmin == 4.0 && kmax == 6.0);
        let nsteps: usize = 100;
        for t in 0..=nsteps {
            let t = kmin + (t as f64 / nsteps as f64) * (kmax - kmin);
            //dbg!(kmin, t, kmax);
            // this curve shouldn't fail to evaluate numerically unless knot index is out of bounds
            let _ = curve.eval(t).unwrap_or_else(|e| {
                panic!("Unexpected error in eval: {:?}", e);
            });
        }

        // check if eval returns the correct error if t is out of bounds
        assert_eq!(
            curve.eval(-0.1).err(),
            Some(BSplineError::KnotDomainViolation)
        );
        assert_eq!(
            curve.eval(9.1).err(),
            Some(BSplineError::KnotDomainViolation)
        );
    }

    #[test]
    fn bspline_api_parity() {
        let points = [PointN::new([0f64, 0f64]), PointN::new([2f64, 0f64])];
        let knots: [f64; 4] = [0.0, 0.0, 1.0, 1.0];
        let curve: BSpline<PointN<f64, 2>, 4, 2, 1> = BSpline::new(knots, points).unwrap();

        let (kmin, kmax) = curve.knot_domain();
        assert_eq!(curve.start().unwrap(), curve.eval(kmin).unwrap());
        assert_eq!(curve.end().unwrap(), curve.eval(kmax).unwrap());

        let reversed = curve.reverse();
        let t = kmin + (kmax - kmin) * 0.25;
        let p0 = curve.eval(t).unwrap();
        let p1 = reversed.eval(kmin + kmax - t).unwrap();
        assert!((p0 - p1).squared_norm() < EPSILON);

        let tangent = curve.tangent(kmin + (kmax - kmin) * 0.5).unwrap();
        assert!((tangent[0] - 1.0).abs() < EPSILON);
        assert!(tangent[1].abs() < EPSILON);

        let expected_len = (points[1] - points[0]).squared_norm().sqrt();
        let t_mid = curve.t_at_length(expected_len * 0.5).unwrap();
        assert!((t_mid - (kmin + kmax) * 0.5).abs() < 1e-3);

        let p_mid = curve.point_at_length(expected_len * 0.5).unwrap();
        assert!((p_mid - PointN::new([1.0, 0.0])).squared_norm() < EPSILON);
    }

    #[test]
    fn basis_functions_linear() {
        let points = [PointN::new([0f64]), PointN::new([2f64])];
        let knots: [f64; 4] = [0.0, 0.0, 1.0, 1.0];
        let curve: BSpline<PointN<f64, 1>, 4, 2, 1> = BSpline::new(knots, points).unwrap();

        let t = 0.25;
        let span = curve.knot_span(t).unwrap();
        let basis = curve.basis_functions(span, t).unwrap();
        assert!((basis[0] - 0.75).abs() < 1e-6);
        assert!((basis[1] - 0.25).abs() < 1e-6);

        let ders = curve.basis_functions_with_derivatives::<1>(span, t).unwrap();
        assert!((ders[0][0] - 0.75).abs() < 1e-6);
        assert!((ders[0][1] - 0.25).abs() < 1e-6);
        assert!((ders[1][0] + 1.0).abs() < 1e-6);
        assert!((ders[1][1] - 1.0).abs() < 1e-6);
    }

    #[test]
    fn bspline_curvature_nonzero() {
        let knots: [f64; 6] = [0.0, 0.0, 0.0, 1.0, 1.0, 1.0];
        let points = [
            PointN::new([0.0, 0.0]),
            PointN::new([1.0, 1.0]),
            PointN::new([2.0, 0.0]),
        ];
        let curve: BSpline<PointN<f64, 2>, 6, 3, 2> = BSpline::new(knots, points).unwrap();

        let curvature = curve.curvature(0.5).unwrap();
        assert!(curvature > 0.0);

        let normal = curve.normal(0.5).unwrap().unwrap();
        assert!((normal.squared_norm() - 1.0).abs() < 1e-6);
    }

    #[test]
    fn root_newton_axis_linear() {
        let points = [PointN::new([0f64]), PointN::new([2f64])];
        let knots: [f64; 4] = [0.0, 0.0, 1.0, 1.0];
        let curve: BSpline<PointN<f64, 1>, 4, 2, 1> = BSpline::new(knots, points).unwrap();

        let root = curve
            .root_newton_axis(1.0, 0, 0.25, None, Some(64))
            .unwrap();
        assert!((root - 0.5).abs() < 1e-6);
    }

    #[test]
    fn root_newton_axis_zero_derivative() {
        let points = [PointN::new([1f64]), PointN::new([1f64])];
        let knots: [f64; 4] = [0.0, 0.0, 1.0, 1.0];
        let curve: BSpline<PointN<f64, 1>, 4, 2, 1> = BSpline::new(knots, points).unwrap();

        let result = curve.root_newton_axis(0.0, 0, 0.5, None, Some(8));
        assert!(matches!(result, Err(RootFindingError::ZeroDerivative)));
    }

    #[test]
    fn derivative_reduces_knots_and_control_points() {
        let points = [PointN::new([0f64]), PointN::new([2f64])];
        let knots: [f64; 4] = [0.0, 0.0, 1.0, 1.0];
        let curve: BSpline<PointN<f64, 1>, 4, 2, 1> = BSpline::new(knots, points).unwrap();

        let derivative = curve.derivative();
        assert_eq!(derivative.knots, [0.0, 1.0]);
        assert_eq!(derivative.control_points, [PointN::new([2.0])]);
    }

    #[test]
    fn bspline_bounding_box_contains_samples() {
        let knots: [f64; 8] = [0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 3.0, 3.0];
        let points = [
            PointN::new([0.0, 0.0]),
            PointN::new([1.0, 2.0]),
            PointN::new([2.0, -1.0]),
            PointN::new([3.0, 1.5]),
            PointN::new([4.0, 0.0]),
        ];
        let curve: BSpline<PointN<f64, 2>, 8, 5, 2> = BSpline::new(knots, points).unwrap();
        let bounds = curve.bounding_box();
        let (kmin, kmax) = curve.knot_domain();
        let steps = 200;

        for i in 0..=steps {
            let t = kmin + (kmax - kmin) * (i as f64 / steps as f64);
            let p = curve.eval(t).unwrap();
            assert!(p[0] >= bounds[0].0 - EPSILON && p[0] <= bounds[0].1 + EPSILON);
            assert!(p[1] >= bounds[1].0 - EPSILON && p[1] <= bounds[1].1 + EPSILON);
        }
    }

    #[test]
    fn domain_contains_matches_knot_domain() {
        let points = [PointN::new([0f64, 0f64]), PointN::new([2f64, 0f64])];
        let knots: [f64; 4] = [0.0, 0.0, 1.0, 1.0];
        let curve: BSpline<PointN<f64, 2>, 4, 2, 1> = BSpline::new(knots, points).unwrap();

        let (kmin, kmax) = curve.knot_domain();
        assert!(curve.domain_contains(kmin));
        assert!(curve.domain_contains(kmax));
        assert!(!curve.domain_contains(kmin - 0.1));
        assert!(!curve.domain_contains(kmax + 0.1));
    }
}
