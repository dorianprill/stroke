//! Const-generic B-spline curves.

use core::slice::*;

use super::{ArcLengthError, Bezier, BoundingBox, Point, PointIndex, PointNorm};
use crate::arc_length::validate_accuracy;
use crate::find_root::FindRoot;
use crate::roots::RootFindingError;
use num_traits::{Float, NumCast};

const DEFAULT_DISTANCE_STEPS: usize = 64;

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
    /// The knot vector contains NaN or an infinite value.
    NonFiniteKnot,
    /// A control point contains NaN or an infinite coordinate.
    NonFiniteControlPoint {
        /// Index of the first non-finite control point.
        index: usize,
    },
    /// The active knot domain has zero width.
    EmptyKnotDomain,
    /// A knot occurs more than degree + 1 times.
    ExcessKnotMultiplicity,
    /// The evaluation parameter is outside the knot domain.
    KnotDomainViolation,
    /// The supplied span does not own the evaluation parameter.
    InvalidSpan,
    /// The de Boor evaluation encountered a NaN alpha value.
    EvaluationAlphaIsNaN,
    /// Differential evaluation produced a non-finite result.
    NonFiniteComputation,
}

#[derive(Debug, Clone, Copy, PartialEq, Default)]
enum KnotVectorKind {
    #[default]
    Unclamped,
    ClampedStart,
    ClampedEnd,
    Clamped,
}

fn classify_knot_vector<T: PartialEq>(knots: &[T], degree: usize) -> KnotVectorKind {
    let is_clamped_front = knots[..=degree].windows(2).all(|pair| pair[0] == pair[1]);
    let is_clamped_back = knots[knots.len() - degree - 1..]
        .windows(2)
        .all(|pair| pair[0] == pair[1]);

    match (is_clamped_front, is_clamped_back) {
        (false, false) => KnotVectorKind::Unclamped,
        (true, false) => KnotVectorKind::ClampedStart,
        (false, true) => KnotVectorKind::ClampedEnd,
        (true, true) => KnotVectorKind::Clamped,
    }
}

impl core::fmt::Display for BSplineError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            BSplineError::TooFewControlPoints => write!(f, "Too few control points"),
            BSplineError::TooFewKnots => write!(f, "Insufficient number of knots"),
            BSplineError::TooManyKnots => write!(f, "Excess number of knots"),
            BSplineError::DecreasingKnots => write!(f, "Knots must be non-decreasing in value"),
            BSplineError::NonFiniteKnot => write!(f, "Knots must be finite"),
            BSplineError::NonFiniteControlPoint { index } => {
                write!(f, "Control point at index {index} must be finite")
            }
            BSplineError::EmptyKnotDomain => {
                write!(f, "The active knot domain must have positive width")
            }
            BSplineError::ExcessKnotMultiplicity => {
                write!(f, "Knot multiplicity exceeds degree + 1")
            }
            BSplineError::KnotDomainViolation => write!(f, "Parameter outside of knot span"),
            BSplineError::InvalidSpan => f.write_str("B-spline span does not own the parameter"),
            BSplineError::EvaluationAlphaIsNaN => {
                write!(f, "Cox-deBoor evaluation yields alpha NaN value")
            }
            BSplineError::NonFiniteComputation => {
                f.write_str("B-spline arithmetic produced a non-finite value")
            }
        }
    }
}

impl core::error::Error for BSplineError {}

/// Const-generic B-spline curve with fixed control-point and knot counts.
///
/// The scalar type used for the knots and interpolation is `P::Scalar`.
///
/// # Parameters
/// - `P`: point type implementing [`Point`].
/// - `C`: number of control points.
/// - `K`: number of knots (`K = C + D + 1`).
/// - `D`: degree (order - 1).
#[derive(Clone, Copy, PartialEq)]
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

impl<P, const K: usize, const C: usize, const D: usize> core::fmt::Debug for BSpline<P, K, C, D>
where
    P: Point + core::fmt::Debug,
    P::Scalar: core::fmt::Debug,
{
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("BSpline")
            .field("knots", &self.knots)
            .field("control_points", &self.control_points)
            .finish()
    }
}

impl<P, const K: usize, const C: usize, const D: usize> Default for BSpline<P, K, C, D>
where
    P: Point + Default,
    P::Scalar: Default,
{
    fn default() -> Self {
        let zero = <P::Scalar as NumCast>::from(0.0).unwrap();
        let one = <P::Scalar as NumCast>::from(1.0).unwrap();
        let knots = core::array::from_fn(|i| {
            if i <= D {
                zero
            } else if i >= C {
                one
            } else {
                <P::Scalar as NumCast>::from(i - D).unwrap()
                    / <P::Scalar as NumCast>::from(C - D).unwrap()
            }
        });
        Self::new(knots, core::array::from_fn(|_| P::default()))
            .expect("default B-splines require a valid layout and finite default point")
    }
}

impl<P, const K: usize, const C: usize, const D: usize> BSpline<P, { K }, { C }, { D }>
where
    P: Point,
{
    /// Create a B-spline defined by the control points and polynomial knot spans.
    /// The curve generally approximates its control polygon rather than interpolating every point.
    /// The knots and control points _must_ be finite, and knots must be sorted in
    /// non-decreasing order.
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
        } else if knots.iter().any(|knot| !knot.is_finite()) {
            Err(BSplineError::NonFiniteKnot)
        } else if let Some(index) = control_points
            .iter()
            .position(|control_point| !control_point.is_finite())
        {
            Err(BSplineError::NonFiniteControlPoint { index })
        // if knots array is not non-decreasing
        } else if knots.windows(2).any(|w| w[0] > w[1]) {
            Err(BSplineError::DecreasingKnots)
        } else if knots[D] >= knots[C] {
            Err(BSplineError::EmptyKnotDomain)
        } else if knots
            .windows(D + 2)
            .any(|window| window[0] == window[D + 1])
        {
            Err(BSplineError::ExcessKnotMultiplicity)
        } else {
            let knot_kind = classify_knot_vector(&knots, D);

            Ok(BSpline {
                knots,
                knot_kind,
                control_points,
            })
        }
    }

    /// Construct an open uniform, fully clamped B-spline on `[0, 1]`.
    /// Requires `C > D` and `K = C + D + 1`; control points must be finite.
    pub fn new_clamped_uniform(control_points: [P; C]) -> Result<Self, BSplineError> {
        if C <= D {
            return Err(BSplineError::TooFewControlPoints);
        }
        let expected = C
            .checked_add(D)
            .and_then(|n| n.checked_add(1))
            .ok_or(BSplineError::TooFewKnots)?;
        if K < expected {
            return Err(BSplineError::TooFewKnots);
        }
        if K > expected {
            return Err(BSplineError::TooManyKnots);
        }
        let zero = <P::Scalar as NumCast>::from(0.0).unwrap();
        let one = <P::Scalar as NumCast>::from(1.0).unwrap();
        let spans = <P::Scalar as NumCast>::from(C - D).unwrap();
        let knots = core::array::from_fn(|i| {
            if i <= D {
                zero
            } else if i >= C {
                one
            } else {
                <P::Scalar as NumCast>::from(i - D).unwrap() / spans
            }
        });
        Self::new(knots, control_points)
    }

    /// Consume the curve, returning its control-point array.
    pub fn into_control_points(self) -> [P; C] {
        self.control_points
    }

    /// Replace control points atomically, retaining the knot vector.
    pub fn set_control_points(&mut self, points: [P; C]) -> Result<(), BSplineError> {
        *self = Self::new(self.knots, points)?;
        Ok(())
    }

    /// Map control points and validate the result without modifying this curve.
    /// Only affine maps geometrically transform the whole curve in the same way.
    pub fn try_map_control_points(&self, map: impl FnMut(P) -> P) -> Result<Self, BSplineError> {
        Self::new(self.knots, self.control_points.map(map))
    }

    /// Evaluate the derivative with respect to the native knot parameter.
    /// At interior knots this returns the right-hand piece derivative; at the final
    /// endpoint it uses the left-hand piece. Degree-zero pieces return zero.
    /// This does not assert differentiability across a knot or jump.
    pub fn derivative_at(&self, t: P::Scalar) -> Result<P, BSplineError>
    where
        [(); D + 1]: Sized,
    {
        let span = self.knot_span(t)?;
        let zero = <P::Scalar as NumCast>::from(0.0).unwrap();
        if D == 0 {
            return Ok(self.control_points[span] * zero);
        }
        let degree = <P::Scalar as NumCast>::from(D).unwrap();
        // Evaluate the degree D-1 derivative with one D+1 workspace. This avoids
        // exposing D-1, C-1 and K-2 const-expression bounds to callers.
        let mut work = [self.control_points[0] * zero; D + 1];
        for (j, value) in work.iter_mut().enumerate().take(D) {
            let i = span - D + j;
            let start = self.knots[i + 1];
            let end = self.knots[i + D + 1];
            *value = if start == end {
                self.control_points[i] * zero
            } else {
                crate::numeric::difference_quotient(
                    self.control_points[i],
                    self.control_points[i + 1],
                    start,
                    end,
                ) * degree
            };
        }
        for r in 1..D {
            for j in (r..D).rev() {
                let i = span - D + j;
                let start = self.knots[i + 1];
                let end = self.knots[i + D - r + 1];
                let alpha = if start == end {
                    zero
                } else {
                    crate::numeric::unit_parameter(t, start, end)
                };
                work[j] = crate::numeric::lerp(work[j - 1], work[j], alpha);
            }
        }
        let result = work[D - 1];
        if result.is_finite() {
            Ok(result)
        } else {
            Err(BSplineError::NonFiniteComputation)
        }
    }

    /// Return derivative magnitude per unit of the native knot parameter.
    pub fn speed(&self, t: P::Scalar) -> Result<P::Scalar, BSplineError>
    where
        P: PointNorm,
        [(); D + 1]: Sized,
    {
        let speed = self.derivative_at(t)?.norm();
        if speed.is_finite() {
            Ok(speed)
        } else {
            Err(BSplineError::NonFiniteComputation)
        }
    }

    /// Compute a point on the curve at `t` using iterative de boor algorithm.
    /// The parameter **must** be in the inclusive range of values returned
    /// by `knot_domain()`. If `t` is out of bounds a `KnotDomainViolation` error is returned.
    pub fn eval(&self, t: P::Scalar) -> Result<P, BSplineError>
    where
        [(); D + 1]: Sized,
    {
        let (knot_start, knot_end) = self.knot_domain();

        // Ensure t is within the knot domain
        if !t.is_finite() || t < knot_start || t > knot_end {
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
        let knots = core::array::from_fn(|i| {
            let knot = self.knots[K - 1 - i];
            let offset = knot - knot_min;
            if offset.is_finite() {
                knot_max - offset
            } else {
                knot_min + (knot_max - knot)
            }
        });
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
    /// The standard active domain is always `[knots[D], knots[C]]`,
    /// independent of endpoint clamping. Knots outside it provide basis support.
    pub fn knot_domain(&self) -> (P::Scalar, P::Scalar) {
        (self.knots[D], self.knots[C])
    }

    /// Return true if `t` lies within the inclusive knot domain.
    ///
    /// # Examples
    /// ```rust
    /// use stroke::{BSpline, ArrayPoint};
    ///
    /// let knots: [f64; 4] = [0.0, 0.0, 1.0, 1.0];
    /// let points = [ArrayPoint::new([0.0]), ArrayPoint::new([2.0])];
    /// let curve: BSpline<ArrayPoint<f64, 1>, 4, 2, 1> = BSpline::new(knots, points).unwrap();
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
    /// use stroke::{BSpline, ArrayPoint};
    ///
    /// let knots: [f64; 4] = [0.0, 0.0, 1.0, 1.0];
    /// let points = [ArrayPoint::new([0.0]), ArrayPoint::new([2.0])];
    /// let curve: BSpline<ArrayPoint<f64, 1>, 4, 2, 1> = BSpline::new(knots, points).unwrap();
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
    /// Returns `InvalidSpan` if it does not own `t`, or `NonFiniteComputation`
    /// if a result cannot be represented by the scalar type.
    ///
    /// # Examples
    /// ```rust
    /// use stroke::{BSpline, ArrayPoint};
    ///
    /// let knots: [f64; 4] = [0.0, 0.0, 1.0, 1.0];
    /// let points = [ArrayPoint::new([0.0]), ArrayPoint::new([2.0])];
    /// let curve: BSpline<ArrayPoint<f64, 1>, 4, 2, 1> = BSpline::new(knots, points).unwrap();
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
        Ok(self.basis_table::<0>(span, t)?[0])
    }

    /// Evaluate the basis functions and their derivatives up to order `R`.
    ///
    /// Returns a `(R + 1) x (D + 1)` table `ders[k][j]` where `k` is the
    /// derivative order and `j` indexes the basis functions for the active span
    /// (`j = 0` corresponds to `i = span - D`). The zeroth row `ders[0]` matches
    /// `basis_functions(span, t)`. Requires `R <= D`.
    /// Span validation matches `basis_functions`; derivatives are in native
    /// knot units. Unrepresentable results return `NonFiniteComputation`.
    ///
    /// # Examples
    /// ```rust
    /// use stroke::{BSpline, ArrayPoint};
    ///
    /// let knots: [f64; 4] = [0.0, 0.0, 1.0, 1.0];
    /// let points = [ArrayPoint::new([0.0]), ArrayPoint::new([2.0])];
    /// let curve: BSpline<ArrayPoint<f64, 1>, 4, 2, 1> = BSpline::new(knots, points).unwrap();
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
        self.basis_table::<R>(span, t)
    }

    /// Validate the same one-sided span ownership as `knot_span` without
    /// adding an unchecked caller-supplied index to the degree.
    fn validate_span(&self, span: usize, t: P::Scalar) -> Result<(), BSplineError> {
        if !self.domain_contains(t) {
            return Err(BSplineError::KnotDomainViolation);
        }
        if !(D..C).contains(&span) {
            return Err(BSplineError::InvalidSpan);
        }
        let start = self.knots[span];
        let end = self.knots[span + 1];
        if start >= end || t < start || (t >= end && !(t == self.knots[C] && t == end)) {
            return Err(BSplineError::InvalidSpan);
        }
        Ok(())
    }

    /// Cox-de Boor recurrence with separate exponents. Derivatives are built
    /// from the preceding degree, delaying conversion to native scalar units
    /// until each complete row is available. Equal knots contribute zero.
    fn basis_table<const R: usize>(
        &self,
        span: usize,
        t: P::Scalar,
    ) -> Result<[[P::Scalar; D + 1]; R + 1], BSplineError>
    where
        [(); D + 1]: Sized,
        [(); R + 1]: Sized,
    {
        use crate::numeric::Scaled;

        self.validate_span(span, t)?;
        let zero = Scaled::new(<P::Scalar as NumCast>::from(0_u8).unwrap());
        let mut table = [[zero; D + 1]; R + 1];
        table[0][0] = Scaled::new(<P::Scalar as NumCast>::from(1_u8).unwrap());
        let parameter = Scaled::new(t);
        for degree in 1..=D {
            let mut next = [[zero; D + 1]; R + 1];
            let factor = Scaled::new(<P::Scalar as NumCast>::from(degree).unwrap());
            for j in 0..=degree {
                let i = span - degree + j;
                if j > 0 && self.knots[i] < self.knots[i + degree] {
                    let lower = Scaled::new(self.knots[i]);
                    let width = Scaled::new(self.knots[i + degree]).difference(lower);
                    let weight = parameter.difference(lower).quotient(width);
                    next[0][j] = next[0][j].sum(table[0][j - 1].product(weight));
                    for order in 1..=R.min(degree) {
                        next[order][j] = next[order][j]
                            .sum(table[order - 1][j - 1].product(factor).quotient(width));
                    }
                }
                if j < degree && self.knots[i + 1] < self.knots[i + degree + 1] {
                    let upper = Scaled::new(self.knots[i + degree + 1]);
                    let width = upper.difference(Scaled::new(self.knots[i + 1]));
                    let weight = upper.difference(parameter).quotient(width);
                    next[0][j] = next[0][j].sum(table[0][j].product(weight));
                    for order in 1..=R.min(degree) {
                        next[order][j] = next[order][j]
                            .difference(table[order - 1][j].product(factor).quotient(width));
                    }
                }
            }
            table = next;
        }
        let result = table.map(|row| row.map(Scaled::value));
        if result.iter().flatten().all(|value| value.is_finite()) {
            Ok(result)
        } else {
            Err(BSplineError::NonFiniteComputation)
        }
    }

    fn bezier_span(&self, span: usize) -> Bezier<P, { D + 1 }>
    where
        [(); D + 1]: Sized,
    {
        let a = self.knots[span];
        let b = self.knots[span + 1];
        let control_points = core::array::from_fn(|j| {
            let mut points: [P; D + 1] =
                core::array::from_fn(|i| self.control_points[span - D + i]);
            for r in 1..=D {
                let parameter = if r <= D - j { a } else { b };
                for i in (r..=D).rev() {
                    let index = span - D + i;
                    let alpha = crate::numeric::unit_parameter(
                        parameter,
                        self.knots[index],
                        self.knots[index + D - r + 1],
                    );
                    points[i] = crate::numeric::lerp(points[i - 1], points[i], alpha);
                }
            }
            points[D]
        });
        Bezier { control_points }
    }

    /// Iterate over the nonzero polynomial spans as Bézier curves.
    ///
    /// Each item contains its original knot-parameter range and a curve
    /// parameterized over `[0, 1]`. At a discontinuous knot, adjacent pieces
    /// expose their respective one-sided endpoints. Extraction uses the
    /// B-spline blossom (polar de Boor algorithm), without allocation.
    pub fn bezier_spans(
        &self,
    ) -> impl Iterator<Item = (core::ops::Range<P::Scalar>, Bezier<P, { D + 1 }>)> + '_
    where
        [(); D + 1]: Sized,
    {
        (D..C).filter_map(move |span| {
            let a = self.knots[span];
            let b = self.knots[span + 1];
            (a < b).then(|| (a..b, self.bezier_span(span)))
        })
    }

    /// Return the tight axis-aligned box of all nonzero polynomial spans.
    pub fn bounding_box(&self) -> BoundingBox<P::Scalar, { P::DIM }>
    where
        P: PointIndex,
        [(); D + 1]: Sized,
    {
        let mut spans = self.bezier_spans();
        let (_, first) = spans.next().expect("a valid B-spline has a nonzero span");
        let mut bounds = first.bounding_box();
        for (_, curve) in spans {
            bounds = bounds.union(&curve.bounding_box());
        }
        bounds
    }

    /// Return a conservative axis-aligned box around all control points.
    pub fn control_box(&self) -> BoundingBox<P::Scalar, { P::DIM }>
    where
        P: PointIndex,
    {
        BoundingBox::from_points(&self.control_points)
    }

    /// Approximate the minimum distance between given `point` and the curve.
    /// Uses a coarse sampling pass over the knot domain and a local search around
    /// the best sample. `nsteps` counts coarse subdivisions (at least one).
    /// Every sampled/refined candidate is retained if better; the search has no
    /// global-minimum or accuracy guarantee. Non-finite query points return NaN;
    /// distances outside the scalar range are infinite.
    pub fn distance_to_point_approx(
        &self,
        point: P,
        nsteps: usize,
    ) -> Result<P::Scalar, BSplineError>
    where
        P: PointNorm,
        [(); D + 1]: Sized,
    {
        let (kmin, kmax) = self.knot_domain();
        let zero = <P::Scalar as NumCast>::from(0.0).unwrap();
        let one = <P::Scalar as NumCast>::from(1.0).unwrap();
        crate::distance::approximate(point, nsteps, |u| {
            // Preserve endpoints and avoid forming an overflowing domain width.
            let t = if u == zero {
                kmin
            } else if u == one {
                kmax
            } else {
                kmin * (one - u) + kmax * u
            };
            self.eval(t)
        })
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
        // The span lookup includes the right endpoint. Evaluate it normally:
        // only end-clamped splines necessarily interpolate the last control.
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
                let alpha = if numerator == zero {
                    // we are at the start of a knot span (interpolation d[j-1] -> d[j])
                    zero
                } else if denominator == zero {
                    // we are at then end of a knot span (interpolation d[j-1] -> d[j])
                    // this will lead to d[j] simply being copied to d[j]
                    one
                } else {
                    // inside a knot span (interpolation d[j-1] -> d[j])
                    crate::numeric::unit_parameter(
                        t,
                        self.knots[j + k - D],
                        self.knots[j + 1 + k - r],
                    )
                };
                if alpha.is_nan() {
                    return Err(BSplineError::EvaluationAlphaIsNaN);
                }
                // The right of two points is overwritten by the result
                // to be used in the next iteration
                d[j] = crate::numeric::lerp(d[j - 1], d[j], alpha);
            }
        }
        Ok(d[D])
    }

    fn knot_span_start_for_t(&self, t: P::Scalar) -> Option<usize> {
        if !self.domain_contains(t) {
            return None;
        }
        if t == self.knots[C] {
            return (D..C).rev().find(|&i| self.knots[i] < self.knots[i + 1]);
        }
        let mut lo = D;
        let mut hi = C;
        while lo < hi {
            let mid = lo + (hi - lo) / 2;
            if t < self.knots[mid] {
                hi = mid;
            } else if t >= self.knots[mid + 1] {
                lo = mid + 1;
            } else {
                return Some(mid);
            }
        }
        None
    }

    /// Approximate length with `segments` lines per nonzero knot span.
    /// Jumps at discontinuous knots do not contribute connecting lines.
    pub fn arclen_polyline(&self, segments: usize) -> P::Scalar
    where
        P: PointNorm,
        [(); D + 1]: Sized,
    {
        self.bezier_spans().fold(
            <P::Scalar as NumCast>::from(0.0).unwrap(),
            |sum, (_, curve)| sum + curve.arclen_polyline(segments),
        )
    }

    /// Integrate length over the knot domain using adaptive Gauss-Legendre quadrature.
    /// `accuracy` is a positive finite absolute length tolerance.
    /// Discontinuous spans contribute their individual lengths, without jumps.
    pub fn arclen(&self, accuracy: P::Scalar) -> Result<P::Scalar, ArcLengthError>
    where
        P: PointNorm,
        [(); D + 1]: Sized,
    {
        let (start, end) = self.knot_domain();
        self.arclen_between(start, end, accuracy)
    }

    /// Integrate length over an ordered parameter interval inside the knot domain.
    pub fn arclen_between(
        &self,
        start: P::Scalar,
        end: P::Scalar,
        accuracy: P::Scalar,
    ) -> Result<P::Scalar, ArcLengthError>
    where
        P: PointNorm,
        [(); D + 1]: Sized,
    {
        validate_accuracy(accuracy)?;
        if !self.domain_contains(start) || !self.domain_contains(end) || start > end {
            return Err(ArcLengthError::InvalidParameter);
        }
        let zero = <P::Scalar as NumCast>::from(0.0).unwrap();
        if start == end {
            return Ok(zero);
        }
        let count = self.knots[D..=C].windows(2).filter(|w| w[0] < w[1]).count();
        let tolerance = (accuracy / <P::Scalar as NumCast>::from(count).unwrap())
            .max(P::Scalar::min_positive_value() * P::Scalar::epsilon());
        let mut total = zero;
        for (range, curve) in self.bezier_spans() {
            let a = start.max(range.start);
            let b = end.min(range.end);
            if a < b {
                let a = crate::numeric::unit_parameter(a, range.start, range.end);
                let b = crate::numeric::unit_parameter(b, range.start, range.end);
                total = total + curve.arclen_between(a, b, tolerance)?;
            }
        }
        if total.is_finite() {
            Ok(total)
        } else {
            Err(ArcLengthError::NonFiniteComputation)
        }
    }

    /// Return a knot parameter at length `s`, clamping to the domain endpoints.
    /// Zero-length pieces are skipped for positive distances.
    pub fn t_at_length(
        &self,
        s: P::Scalar,
        accuracy: P::Scalar,
    ) -> Result<P::Scalar, ArcLengthError>
    where
        P: PointNorm,
        [(); D + 1]: Sized,
    {
        validate_accuracy(accuracy)?;
        if !s.is_finite() {
            return Err(ArcLengthError::InvalidParameter);
        }
        let (start, end) = self.knot_domain();
        let zero = <P::Scalar as NumCast>::from(0.0).unwrap();
        if s <= zero {
            return Ok(start);
        }
        let count = self.knots[D..=C].windows(2).filter(|w| w[0] < w[1]).count();
        let tolerance = (accuracy / <P::Scalar as NumCast>::from(count + 1).unwrap())
            .max(P::Scalar::min_positive_value() * P::Scalar::epsilon());
        let mut remaining = s;
        let mut any_length = false;
        for (range, curve) in self.bezier_spans() {
            let length = curve.arclen(tolerance)?;
            any_length |= length > zero;
            if length > zero && remaining < length {
                let local = curve.t_at_length(remaining, tolerance)?;
                let one = <P::Scalar as NumCast>::from(1.0).unwrap();
                return Ok(range.start * (one - local) + range.end * local);
            }
            remaining = remaining - length;
        }
        Ok(if any_length { end } else { start })
    }

    /// Evaluate a point at length `s` using adaptive Gauss-Legendre quadrature.
    pub fn point_at_length(&self, s: P::Scalar, accuracy: P::Scalar) -> Result<P, ArcLengthError>
    where
        P: PointNorm,
        [(); D + 1]: Sized,
    {
        self.eval(self.t_at_length(s, accuracy)?)
            .map_err(|_| ArcLengthError::CurveEvaluationFailed)
    }

    /// Returns the derivative curve of self which has C-1 control points and K-2 knots.
    /// The derivative of an nth degree B-Spline curve is an (n-1)th degree (d) B-Spline curve,
    /// with the same knot vector, and new control points Q0...Qn-1 derived from the
    /// original control points Pi as:
    ///                 d
    /// Qi =    ----------------- (`P[i+1]` - `P[i]`)
    ///         `k[i + d + 1]` - `k[i + 1]`.
    /// with `degree = curve_order - 1`
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
            let scale = if denom == zero {
                zero
            } else {
                <P::Scalar as NumCast>::from(D as f64).unwrap() / denom
            };
            (self.control_points[i + 1] - self.control_points[i]) * scale
        });

        let knot_kind = classify_knot_vector(&derivative_knots, D - 1);

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
        [(); D + 1]: Sized,
    {
        let span = self.knot_span(t)?;
        let local = crate::numeric::unit_parameter(t, self.knots[span], self.knots[span + 1]);
        let result = self.bezier_span(span).tangent(local);
        if result.is_finite() {
            Ok(result)
        } else {
            Err(BSplineError::NonFiniteComputation)
        }
    }

    /// Return the curvature magnitude at `t`.
    ///
    /// Constant and linear pieces have zero curvature.
    pub fn curvature(&self, t: P::Scalar) -> Result<P::Scalar, BSplineError>
    where
        P: PointNorm,
        [(); D + 1]: Sized,
    {
        let span = self.knot_span(t)?;
        let local = crate::numeric::unit_parameter(t, self.knots[span], self.knots[span + 1]);
        let result = self.bezier_span(span).curvature(local);
        if result.is_finite() {
            Ok(result)
        } else {
            Err(BSplineError::NonFiniteComputation)
        }
    }

    /// Return the principal normal direction at `t`.
    ///
    /// Returns `None` for stationary or straight pieces.
    pub fn normal(&self, t: P::Scalar) -> Result<Option<P>, BSplineError>
    where
        P: PointNorm,
        [(); D + 1]: Sized,
    {
        let span = self.knot_span(t)?;
        let local = crate::numeric::unit_parameter(t, self.knots[span], self.knots[span + 1]);
        Ok(self.bezier_span(span).normal(local))
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
        [(); D + 1]: Sized,
        P: PointIndex,
    {
        FindRoot::root_newton_axis(self, value, axis, start, eps, max_iter)
    }
}

impl<P, const K: usize, const C: usize, const D: usize> FindRoot<P> for BSpline<P, K, C, D>
where
    P: PointIndex,
    [(); D + 1]: Sized,
{
    fn parameter_domain(&self) -> (P::Scalar, P::Scalar) {
        self.knot_domain()
    }

    fn axis_value(&self, t: P::Scalar, axis: usize) -> Result<P::Scalar, RootFindingError> {
        if axis >= P::DIM {
            return Err(RootFindingError::InvalidAxis);
        }
        if !t.is_finite() || !self.domain_contains(t) {
            return Err(RootFindingError::InvalidParameter);
        }
        self.eval(t)
            .map(|p| p[axis])
            .map_err(|_| RootFindingError::FailedToConverge)
    }

    fn axis_derivative(&self, t: P::Scalar, axis: usize) -> Result<P::Scalar, RootFindingError> {
        if axis >= P::DIM {
            return Err(RootFindingError::InvalidAxis);
        }
        if !t.is_finite() || !self.domain_contains(t) {
            return Err(RootFindingError::InvalidParameter);
        }
        self.derivative_at(t)
            .map(|p| p[axis])
            .map_err(|_| RootFindingError::NonFiniteComputation)
    }
}

#[cfg(test)]
mod tests {
    //use std;
    use super::*;
    //use crate::num_traits::{Pow};
    use crate::ArrayPoint;
    use crate::test_support::{assert_approx_equal, assert_point_approx_equal};

    const ABSOLUTE_TOLERANCE: f64 = 1e-12;
    const RELATIVE_TOLERANCE: f64 = 1e-12;

    #[test]
    fn rejects_non_finite_knots() {
        let points = [ArrayPoint::new([0.0_f64]), ArrayPoint::new([1.0])];

        for invalid in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let result =
                BSpline::<ArrayPoint<f64, 1>, 4, 2, 1>::new([0.0, invalid, 1.0, 1.0], points);
            assert_eq!(result.err(), Some(BSplineError::NonFiniteKnot));
        }
    }

    #[test]
    fn rejects_non_finite_control_points_with_index() {
        let points = [
            ArrayPoint::new([0.0_f64]),
            ArrayPoint::new([f64::NAN]),
            ArrayPoint::new([f64::INFINITY]),
        ];
        let result =
            BSpline::<ArrayPoint<f64, 1>, 6, 3, 2>::new([0.0, 0.0, 0.0, 1.0, 1.0, 1.0], points);

        assert_eq!(
            result.err(),
            Some(BSplineError::NonFiniteControlPoint { index: 1 })
        );
    }

    #[test]
    fn knot_clamping_is_classified_structurally() {
        let points = [ArrayPoint::new([0.0_f64]), ArrayPoint::new([1.0])];
        let curve = BSpline::<ArrayPoint<f64, 1>, 4, 2, 1>::new(
            [0.0, f64::EPSILON / 2.0, 1.0, 1.0],
            points,
        )
        .unwrap();

        assert_eq!(curve.knot_kind, KnotVectorKind::ClampedEnd);
    }

    #[test]
    fn narrow_knot_span_remains_non_degenerate() {
        let span = 1e-20_f64;
        let points = [ArrayPoint::new([0.0]), ArrayPoint::new([1.0])];
        let curve =
            BSpline::<ArrayPoint<f64, 1>, 4, 2, 1>::new([0.0, 0.0, span, span], points).unwrap();
        let t = span * 0.5;

        assert_point_approx_equal(
            curve.eval(t).unwrap(),
            ArrayPoint::new([0.5]),
            ABSOLUTE_TOLERANCE,
            RELATIVE_TOLERANCE,
        );

        let knot_span = curve.knot_span(t).unwrap();
        let basis = curve.basis_functions(knot_span, t).unwrap();
        assert_approx_equal(basis[0], 0.5, ABSOLUTE_TOLERANCE, RELATIVE_TOLERANCE);
        assert_approx_equal(basis[1], 0.5, ABSOLUTE_TOLERANCE, RELATIVE_TOLERANCE);
        assert_point_approx_equal(
            curve.tangent(t).unwrap(),
            ArrayPoint::new([1.0]),
            ABSOLUTE_TOLERANCE,
            RELATIVE_TOLERANCE,
        );
    }

    #[test]
    fn degree_1_clamped_construct_and_eval_endpoints() {
        // Define the control points for a degree 1 B-Spline (essentially a line)
        let points = [
            ArrayPoint::new([0f64, 0f64]),   // Start of the line at origin
            ArrayPoint::new([10f64, 10f64]), // End of the line at (10, 10)
        ];

        // Define a uniform clamped knot vector for a degree 1 B-spline
        let knots: [f64; 4] = [0., 0., 1., 1.];

        // Create the B-Spline
        let curve: BSpline<ArrayPoint<f64, 2>, 4, 2, 1> = BSpline::new(knots, points).unwrap();

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
        assert_point_approx_equal(
            mid,
            ArrayPoint::new([5.0, 5.0]),
            ABSOLUTE_TOLERANCE,
            RELATIVE_TOLERANCE,
        );
    }

    #[test]
    fn degree_1_unclamped_construct_and_eval() {
        // Define the control points for a degree 1 B-Spline (essentially a line)
        let points = [
            ArrayPoint::new([0f64, 0f64]),   // Start of the line at origin
            ArrayPoint::new([10f64, 10f64]), // End of the line at (10, 10)
        ];

        // Define a uniform unclamped knot vector for a degree 1 B-spline
        let knots: [f64; 4] = [0., 1., 2., 3.];

        // Create the B-Spline
        let curve: BSpline<ArrayPoint<f64, 2>, 4, 2, 1> = BSpline::new(knots, points).unwrap();

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
        assert_point_approx_equal(
            mid,
            ArrayPoint::new([5.0, 5.0]),
            ABSOLUTE_TOLERANCE,
            RELATIVE_TOLERANCE,
        );
    }

    #[test]
    fn unclamped_cubic_endpoint_matches_polynomial_segment() {
        // On [3, 4], these uniform knots and control values define (t - 3)^3.
        // In particular, the final value is 1, rather than the last control's 6.
        let curve = BSpline::<ArrayPoint<f64, 1>, 8, 4, 3>::new(
            [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0],
            [
                ArrayPoint::new([0.0]),
                ArrayPoint::new([0.0]),
                ArrayPoint::new([0.0]),
                ArrayPoint::new([6.0]),
            ],
        )
        .unwrap();
        assert_eq!(curve.knot_domain(), (3.0, 4.0));
        for t in [3.0, 3.5, 4.0 - 1e-8, 4.0] {
            let u = t - 3.0;
            assert!((curve.eval(t).unwrap()[0] - u * u * u).abs() < 1e-12);
        }
        assert!((curve.end().unwrap()[0] - 1.0).abs() < 1e-12);
        assert!((curve.reverse().start().unwrap()[0] - 1.0).abs() < 1e-12);
        assert_eq!(curve.bounding_box().into_axes(), [(0.0, 1.0)]);
    }

    #[test]
    fn arclen_line_approx() {
        let points = [
            ArrayPoint::new([0f64, 0f64]),
            ArrayPoint::new([10f64, 0f64]),
        ];
        let knots: [f64; 4] = [0.0, 0.0, 1.0, 1.0];
        let curve: BSpline<ArrayPoint<f64, 2>, 4, 2, 1> = BSpline::new(knots, points).unwrap();

        let expected = (points[1] - points[0]).squared_norm().sqrt();
        let length = curve.arclen_polyline(32);
        let tolerance = 1e-6;
        assert!((length - expected).abs() <= tolerance);
    }

    #[test]
    fn degree_2_clamped_construct_and_eval_endpoints() {
        // Define the control points for a degree 2 B-Spline
        let points = [
            ArrayPoint::new([0f64, 0f64]),   // Start of the curve at origin
            ArrayPoint::new([10f64, 10f64]), // Control point off the baseline
            ArrayPoint::new([20f64, 0f64]),  // End of the curve at (20, 0)
        ];

        // Define a uniform clamped knot vector for a degree 2 B-spline
        let knots: [f64; 6] = [0., 0., 0., 1., 1., 1.];

        // Create the B-Spline
        let curve: BSpline<ArrayPoint<f64, 2>, 6, 3, 2> = BSpline::new(knots, points).unwrap();

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
        assert_point_approx_equal(
            mid,
            ArrayPoint::new([10.0, 5.0]),
            ABSOLUTE_TOLERANCE,
            RELATIVE_TOLERANCE,
        );
    }

    #[test]
    fn degree_2_unclamped_construct_and_eval() {
        // Define the control points for a degree 2 B-Spline
        let points = [
            ArrayPoint::new([0f64, 0f64]),   // Start of the line at origin
            ArrayPoint::new([10f64, 10f64]), // End of the line at (10, 10)
            ArrayPoint::new([20f64, 0f64]),  // End of the line at (20, 0)
        ];

        // Define a uniform unclamped knot vector for a degree 2 B-spline
        let knots: [f64; 6] = [0., 1., 2., 3., 4., 5.];

        // Create the B-Spline
        let curve: BSpline<ArrayPoint<f64, 2>, 6, 3, 2> = BSpline::new(knots, points).unwrap();

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
        assert_point_approx_equal(mid, expected, ABSOLUTE_TOLERANCE, RELATIVE_TOLERANCE);
    }

    #[test]
    fn degree_3_clamped_construct_and_eval_endpoints() {
        // degree 3, 4 control points => 4+3+1=8 knots
        let points = [
            ArrayPoint::new([0f64, 1.77f64]),
            ArrayPoint::new([1.1f64, -1f64]),
            ArrayPoint::new([4.3f64, 3f64]),
            ArrayPoint::new([3.2f64, -4f64]),
            ArrayPoint::new([2.3f64, 1.77f64]),
        ];
        // clamped b-spline to it goes through the first and last control point
        let knots: [f64; 9] = [0., 0., 0., 0., 1., 2., 2., 2., 2.];

        // try to make a b-spline with the given parameters
        let b: Result<BSpline<ArrayPoint<f64, 2>, 9, 5, 3>, BSplineError> =
            BSpline::new(knots, points);

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
            ArrayPoint::new([0f64, 1.77f64]),
            ArrayPoint::new([1.1f64, -1f64]),
            ArrayPoint::new([4.3f64, 3f64]),
            ArrayPoint::new([3.2f64, -4f64]),
            ArrayPoint::new([2.3f64, 1.77f64]),
        ];
        // unclamped b-spline (not necessarily going through the first and last control point)
        let knots: [f64; 9] = [0., 1., 2., 3., 4., 5., 6., 7., 8.];

        // try to make a b-spline with the given parameters
        let b: Result<BSpline<ArrayPoint<f64, 2>, 9, 5, 3>, BSplineError> =
            BSpline::new(knots, points);

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
    fn clamped_quintic_matches_bezier() {
        let points = [
            ArrayPoint::new([0.0, 0.0]),
            ArrayPoint::new([0.5, 4.0]),
            ArrayPoint::new([2.0, -3.0]),
            ArrayPoint::new([3.0, 4.0]),
            ArrayPoint::new([4.5, -2.0]),
            ArrayPoint::new([5.0, 1.0]),
        ];
        let spline = BSpline::<ArrayPoint<f64, 2>, 12, 6, 5>::new(
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
            points,
        )
        .unwrap();
        let bezier = crate::Bezier::new(points).unwrap();

        assert_eq!(spline.knot_kind, KnotVectorKind::Clamped);
        assert_eq!(spline.knot_domain(), (0.0, 1.0));
        assert_eq!(spline.start().unwrap(), points[0]);
        assert_eq!(spline.end().unwrap(), points[5]);
        for i in 0..=1000 {
            let t = i as f64 / 1000.0;
            assert!((spline.eval(t).unwrap() - bezier.eval(t)).squared_norm() < 1e-24);
        }
    }

    #[test]
    fn degree_4_clamped_construct_and_eval() {
        // degree 4, 5 control points => 5+4+1=10 knots
        let points = [
            ArrayPoint::new([0f64, 1.77f64]),
            ArrayPoint::new([1.1f64, -1f64]),
            ArrayPoint::new([4.3f64, 3f64]),
            ArrayPoint::new([3.2f64, -4f64]),
            ArrayPoint::new([2.3f64, 1.77f64]),
        ];
        // clamped b-spline (goes through the first and last control point)
        let knots: [f64; 10] = [0., 0., 0., 0., 0., 3., 3., 3., 3., 3.];

        // try to make a b-spline with the given parameters
        let b: Result<BSpline<ArrayPoint<f64, 2>, 10, 5, 4>, BSplineError> =
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
            ArrayPoint::new([0f64, 1.77f64]),
            ArrayPoint::new([1.1f64, -1f64]),
            ArrayPoint::new([4.3f64, 3f64]),
            ArrayPoint::new([3.2f64, -4f64]),
            ArrayPoint::new([2.3f64, 1.77f64]),
            ArrayPoint::new([0.0, 0.0]),
        ];
        // unclamped b-spline (not necessarily going through the first and last control point)
        let knots: [f64; 11] = [0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10.];

        // try to make a b-spline with the given parameters
        let b: Result<BSpline<ArrayPoint<f64, 2>, 11, 6, 4>, BSplineError> =
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
        let points = [ArrayPoint::new([0f64, 0f64]), ArrayPoint::new([2f64, 0f64])];
        let knots: [f64; 4] = [0.0, 0.0, 1.0, 1.0];
        let curve: BSpline<ArrayPoint<f64, 2>, 4, 2, 1> = BSpline::new(knots, points).unwrap();

        let (kmin, kmax) = curve.knot_domain();
        assert_eq!(curve.start().unwrap(), curve.eval(kmin).unwrap());
        assert_eq!(curve.end().unwrap(), curve.eval(kmax).unwrap());

        let reversed = curve.reverse();
        let t = kmin + (kmax - kmin) * 0.25;
        let p0 = curve.eval(t).unwrap();
        let p1 = reversed.eval(kmin + kmax - t).unwrap();
        assert_point_approx_equal(p0, p1, ABSOLUTE_TOLERANCE, RELATIVE_TOLERANCE);

        let tangent = curve.tangent(kmin + (kmax - kmin) * 0.5).unwrap();
        assert_approx_equal(tangent[0], 1.0, ABSOLUTE_TOLERANCE, RELATIVE_TOLERANCE);
        assert_approx_equal(tangent[1], 0.0, ABSOLUTE_TOLERANCE, RELATIVE_TOLERANCE);

        let expected_len = (points[1] - points[0]).squared_norm().sqrt();
        let t_mid = curve.t_at_length(expected_len * 0.5, 1e-12).unwrap();
        assert!((t_mid - (kmin + kmax) * 0.5).abs() < 1e-3);

        let p_mid = curve.point_at_length(expected_len * 0.5, 1e-12).unwrap();
        assert_point_approx_equal(p_mid, ArrayPoint::new([1.0, 0.0]), 1e-9, 1e-9);
    }

    #[test]
    fn basis_functions_linear() {
        let points = [ArrayPoint::new([0f64]), ArrayPoint::new([2f64])];
        let knots: [f64; 4] = [0.0, 0.0, 1.0, 1.0];
        let curve: BSpline<ArrayPoint<f64, 1>, 4, 2, 1> = BSpline::new(knots, points).unwrap();

        let t = 0.25;
        let span = curve.knot_span(t).unwrap();
        let basis = curve.basis_functions(span, t).unwrap();
        assert!((basis[0] - 0.75).abs() < 1e-6);
        assert!((basis[1] - 0.25).abs() < 1e-6);

        let ders = curve
            .basis_functions_with_derivatives::<1>(span, t)
            .unwrap();
        assert!((ders[0][0] - 0.75).abs() < 1e-6);
        assert!((ders[0][1] - 0.25).abs() < 1e-6);
        assert!((ders[1][0] + 1.0).abs() < 1e-6);
        assert!((ders[1][1] - 1.0).abs() < 1e-6);
    }

    #[test]
    fn bspline_curvature_nonzero() {
        let knots: [f64; 6] = [0.0, 0.0, 0.0, 1.0, 1.0, 1.0];
        let points = [
            ArrayPoint::new([0.0, 0.0]),
            ArrayPoint::new([1.0, 1.0]),
            ArrayPoint::new([2.0, 0.0]),
        ];
        let curve: BSpline<ArrayPoint<f64, 2>, 6, 3, 2> = BSpline::new(knots, points).unwrap();

        let curvature = curve.curvature(0.5).unwrap();
        assert!(curvature > 0.0);

        let normal = curve.normal(0.5).unwrap().unwrap();
        assert!((normal.squared_norm() - 1.0).abs() < 1e-6);
    }

    #[test]
    fn root_newton_axis_linear() {
        let points = [ArrayPoint::new([0f64]), ArrayPoint::new([2f64])];
        let knots: [f64; 4] = [0.0, 0.0, 1.0, 1.0];
        let curve: BSpline<ArrayPoint<f64, 1>, 4, 2, 1> = BSpline::new(knots, points).unwrap();

        let root = curve
            .root_newton_axis(1.0, 0, 0.25, None, Some(64))
            .unwrap();
        assert!((root - 0.5).abs() < 1e-6);
    }

    #[test]
    fn root_newton_axis_zero_derivative() {
        let points = [ArrayPoint::new([1f64]), ArrayPoint::new([1f64])];
        let knots: [f64; 4] = [0.0, 0.0, 1.0, 1.0];
        let curve: BSpline<ArrayPoint<f64, 1>, 4, 2, 1> = BSpline::new(knots, points).unwrap();

        let result = curve.root_newton_axis(0.0, 0, 0.5, None, Some(8));
        assert!(matches!(result, Err(RootFindingError::ZeroDerivative)));
    }

    #[test]
    fn derivative_reduces_knots_and_control_points() {
        let points = [ArrayPoint::new([0f64]), ArrayPoint::new([2f64])];
        let knots: [f64; 4] = [0.0, 0.0, 1.0, 1.0];
        let curve: BSpline<ArrayPoint<f64, 1>, 4, 2, 1> = BSpline::new(knots, points).unwrap();

        let derivative = curve.derivative();
        assert_eq!(derivative.knots, [0.0, 1.0]);
        assert_eq!(derivative.control_points, [ArrayPoint::new([2.0])]);
    }

    #[test]
    fn bspline_bounding_box_contains_samples() {
        let knots: [f64; 8] = [0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 3.0, 3.0];
        let points = [
            ArrayPoint::new([0.0, 0.0]),
            ArrayPoint::new([1.0, 2.0]),
            ArrayPoint::new([2.0, -1.0]),
            ArrayPoint::new([3.0, 1.5]),
            ArrayPoint::new([4.0, 0.0]),
        ];
        let curve: BSpline<ArrayPoint<f64, 2>, 8, 5, 2> = BSpline::new(knots, points).unwrap();
        let bounds = curve.bounding_box();
        let (kmin, kmax) = curve.knot_domain();
        let steps = 200;

        for i in 0..=steps {
            let t = kmin + (kmax - kmin) * (i as f64 / steps as f64);
            let p = curve.eval(t).unwrap();
            assert!(
                p[0] >= bounds[0].0 - ABSOLUTE_TOLERANCE
                    && p[0] <= bounds[0].1 + ABSOLUTE_TOLERANCE
            );
            assert!(
                p[1] >= bounds[1].0 - ABSOLUTE_TOLERANCE
                    && p[1] <= bounds[1].1 + ABSOLUTE_TOLERANCE
            );
        }
    }

    #[test]
    fn domain_contains_matches_knot_domain() {
        let points = [ArrayPoint::new([0f64, 0f64]), ArrayPoint::new([2f64, 0f64])];
        let knots: [f64; 4] = [0.0, 0.0, 1.0, 1.0];
        let curve: BSpline<ArrayPoint<f64, 2>, 4, 2, 1> = BSpline::new(knots, points).unwrap();

        let (kmin, kmax) = curve.knot_domain();
        assert!(curve.domain_contains(kmin));
        assert!(curve.domain_contains(kmax));
        assert!(!curve.domain_contains(kmin - 0.1));
        assert!(!curve.domain_contains(kmax + 0.1));
    }
}
