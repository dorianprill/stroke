//! Const-generic B-spline curves.

use core::slice::*;

use super::{
    ArcLengthError, Bezier, BoundingBox, OrientedBoundingBox, OrientedBoundingBoxError, Point,
    PointIndex, PointNorm,
};
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
    /// Curve arithmetic produced a non-finite result.
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

/// A B-spline with `C` control points, `K` knots and degree `D`.
///
/// Knots use the point's scalar type, `P::Scalar`.
///
/// # Parameters
///
/// - `P`: point type implementing [`Point`].
/// - `C`: number of control points.
/// - `K`: number of knots. Must equal `C + D + 1`.
/// - `D`: polynomial degree, one less than the order.
#[derive(Clone, Copy, PartialEq)]
pub struct BSpline<P, const K: usize, const C: usize, const D: usize>
where
    P: Point,
{
    /// Knots in non-decreasing order.
    knots: [P::Scalar; K],
    /// Clamping at the start, end, both or neither.
    knot_kind: KnotVectorKind,
    /// Control points in parameter order.
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
    /// Create a B-spline from its knots and control points.
    ///
    /// Requires `C > D`, `K = C + D + 1`, and finite knots and controls.
    /// Knots must be in non-decreasing order, with at most `D + 1` copies of
    /// any value. The active parameter range `[knots[D], knots[C]]` must have
    /// positive width. Invalid inputs return [`BSplineError`].
    ///
    /// The controls shape the curve. The curve need not pass through them.
    pub fn new(
        knots: [P::Scalar; K],
        control_points: [P; C],
    ) -> Result<BSpline<P, { K }, { C }, { D }>, BSplineError> {
        if control_points.len() <= { D } {
            Err(BSplineError::TooFewControlPoints)
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

    /// Return the control-point array, consuming the curve.
    pub fn into_control_points(self) -> [P; C] {
        self.control_points
    }

    /// Replace the controls and keep the knots. Leave the curve unchanged on error.
    pub fn set_control_points(&mut self, points: [P; C]) -> Result<(), BSplineError> {
        *self = Self::new(self.knots, points)?;
        Ok(())
    }

    /// Return a curve with mapped controls and the same knots. Reject non-finite results.
    /// Affine maps, such as translations, also map every point of the curve.
    pub fn try_map_control_points(&self, map: impl FnMut(P) -> P) -> Result<Self, BSplineError> {
        Self::new(self.knots, self.control_points.map(map))
    }

    /// Evaluate the derivative with respect to the knot parameter `t`.
    /// At an interior knot, use the piece to the right. At the final endpoint,
    /// use the piece to the left. Derivatives across a jump are undefined.
    /// Degree-zero pieces return zero.
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
        // A D+1 workspace holds the degree D-1 derivative. Callers then need
        // no additional D-1, C-1 or K-2 const-expression bounds.
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

    /// Return speed per unit of knot parameter, the magnitude of `derivative_at(t)`.
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

    /// Evaluate at `t` using de Boor's algorithm.
    /// `t` must be within the inclusive range returned by [`Self::knot_domain`].
    /// Otherwise, return [`BSplineError::KnotDomainViolation`].
    pub fn eval(&self, t: P::Scalar) -> Result<P, BSplineError>
    where
        [(); D + 1]: Sized,
    {
        let (knot_start, knot_end) = self.knot_domain();

        if !t.is_finite() || t < knot_start || t > knot_end {
            return Err(BSplineError::KnotDomainViolation);
        }

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

    /// Return the inclusive parameter range `[knots[D], knots[C]]`.
    /// This range applies to clamped and unclamped splines. Knots outside it
    /// still contribute to the basis functions.
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
    /// The index `i` satisfies `knots[i] <= t < knots[i + 1]`. The final domain
    /// endpoint belongs to the last nonzero span. Pass this index to the basis
    /// function methods.
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
    /// Return `D + 1` values, ordered from basis index `span - D` to `span`.
    /// Obtain `span` from [`Self::knot_span`]. A span that does not contain `t`
    /// returns `InvalidSpan`. An unrepresentable result returns `NonFiniteComputation`.
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
    /// Return a `(R + 1) x (D + 1)` table `ders[k][j]`. Row `k` contains the
    /// derivatives of order `k` with respect to the knot parameter. Column `j`
    /// corresponds to basis index `span - D + j`. Row zero matches
    /// [`Self::basis_functions`]. Requires `R <= D`.
    ///
    /// Span validation matches [`Self::basis_functions`]. An unrepresentable
    /// result returns `NonFiniteComputation`.
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

    /// Check that `span` is the span selected by `knot_span(t)`.
    /// Validate the index before adding it to the degree.
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

    /// Evaluate the Cox-de Boor recurrence with separate binary exponents.
    /// Build derivatives from the preceding degree and convert to `P::Scalar`
    /// after each row is complete. Equal knots contribute zero.
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
    /// Each item contains the original knot-parameter range and a Bezier curve
    /// with parameter range `[0, 1]`. At a jump, adjacent pieces retain their
    /// own endpoints, so both one-sided limits are available. Extraction uses
    /// the polar form of de Boor's algorithm, without allocation.
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

    /// Fit a box with rotated axes around the entire active knot domain.
    ///
    /// Use the [Bezier search](Bezier::bounding_box_oriented_approx) with up to
    /// 256 full-curve trials. Include every nonzero span and both limits at jumps.
    /// Each trial combines the span bounds in one set of box axes.
    ///
    /// The search may stop before reaching even a local minimum. It has no
    /// guarantee on volume error. It uses fixed stack storage and needs no PCA
    /// or eigensolver. Floating-point margins include span-extraction work,
    /// but they are not certified interval bounds.
    /// See the [method and research references](
    /// crate::oriented_bounding_box#method-and-references) and the numerical
    /// limits in [`Bezier::bounding_box_oriented_approx`].
    ///
    /// Return [`OrientedBoundingBoxError::ZeroDimension`] for zero-dimensional
    /// points, or [`OrientedBoundingBoxError::NonFiniteComputation`] if
    /// normalization, extraction, projection or padding cannot be represented.
    ///
    /// ```
    /// #![allow(incomplete_features)]
    /// #![feature(generic_const_exprs)]
    /// use stroke::{ArrayPoint, BSpline};
    /// let curve = BSpline::<_, 7, 4, 2>::new_clamped_uniform([
    ///     ArrayPoint::from([0.0_f64, 0.0]),
    ///     ArrayPoint::from([1.0, 3.0]),
    ///     ArrayPoint::from([3.0, -1.0]),
    ///     ArrayPoint::from([4.0, 1.0]),
    /// ]).unwrap();
    /// let bounds = curve.bounding_box_oriented_approx().unwrap();
    /// assert!(bounds.contains(&curve.eval(0.5).unwrap()));
    /// ```
    pub fn bounding_box_oriented_approx(
        &self,
    ) -> Result<OrientedBoundingBox<P::Scalar, { P::DIM }>, OrientedBoundingBoxError>
    where
        P: PointIndex,
        [(); D + 1]: Sized,
        [(); C + 64]: Sized,
    {
        self.bounding_box_oriented_approx_with_budget(256)
    }

    /// Fit an oriented box with at most `max_evaluations` full-curve trials.
    ///
    /// Use [`Self::bounding_box_oriented_approx`] with a caller-chosen budget.
    /// Each trial bounds every nonzero span. The budget limits trials, not runtime
    /// or volume error. Larger budgets extend the same sequence and keep the best
    /// earlier box. This does not establish convergence.
    ///
    /// For budgets greater than one, model construction and its initial search
    /// have a fixed cost. The model retains all spans and caches at most `C + 64`
    /// control polygons on the stack, using `O((C + 64) * (D + 1) * P::DIM)` scalars.
    /// Remaining slots allow subdivision. Refinement targets a 0.1% model-volume
    /// gap at selected frames. Capacity or scalar resolution can prevent reaching
    /// the target. It does not bound the error relative to the minimum-volume box.
    /// [`OrientedBoundingBox::diagnostics`] reports model work and uncertainty.
    ///
    /// Constant controls use zero trials. One-dimensional curves and budget-one
    /// calls skip model construction. A zero budget returns
    /// [`OrientedBoundingBoxError::ZeroBudget`]. Other errors and numerical limits
    /// match [`Self::bounding_box_oriented_approx`].
    pub fn bounding_box_oriented_approx_with_budget(
        &self,
        max_evaluations: usize,
    ) -> Result<OrientedBoundingBox<P::Scalar, { P::DIM }>, OrientedBoundingBoxError>
    where
        P: PointIndex,
        [(); D + 1]: Sized,
        [(); C + 64]: Sized,
    {
        crate::oriented_bounding_box::fit_bspline(self, max_evaluations)
    }

    /// Return a conservative axis-aligned box around all control points.
    pub fn control_box(&self) -> BoundingBox<P::Scalar, { P::DIM }>
    where
        P: PointIndex,
    {
        BoundingBox::from_points(&self.control_points)
    }

    /// Estimate the minimum distance to `point`.
    /// Sample at `nsteps` equal intervals over the knot domain, then refine near
    /// the closest sample. Use at least one interval and retain the best distance.
    /// This local search can miss a closer point and has no error bound.
    /// Non-finite query points return NaN. Distances outside the scalar range are infinite.
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

    /// Estimate the minimum distance using the default number of parameter intervals.
    pub fn distance_to_point(&self, point: P) -> Result<P::Scalar, BSplineError>
    where
        P: PointNorm,
        [(); D + 1]: Sized,
    {
        self.distance_to_point_approx(point, DEFAULT_DISTANCE_STEPS)
    }

    /// Evaluate de Boor's recurrence in a fixed workspace of `D + 1` points.
    fn de_boor_iterative(&self, t: P::Scalar) -> Result<P, BSplineError>
    where
        [(); D + 1]: Sized,
    {
        // Evaluate the final parameter with the recurrence too. An unclamped
        // spline need not end at its last control point.
        let k = self.knot_span_start_for_t(t);

        let k = match k {
            Some(r) => r,
            None => return Err(BSplineError::KnotDomainViolation),
        };

        let mut d: [P; D + 1] = core::array::from_fn(|_| self.control_points[0]);

        // A degree-D span uses D+1 control points.
        for j in 0..=D {
            d[j] = self.control_points[j + k - D];
        }

        // Overwrite the right operand with each interpolation result.
        // Reverse iteration preserves the operands needed by the next pair.
        for r in 1..=D {
            for j in (r..=D).rev() {
                let numerator = t - self.knots[j + k - D];
                let denominator = self.knots[j + 1 + k - r] - self.knots[j + k - D];
                let zero = <P::Scalar as NumCast>::from(0.0).unwrap();
                let one = <P::Scalar as NumCast>::from(1.0).unwrap();
                let alpha = if numerator == zero {
                    zero
                } else if denominator == zero {
                    one
                } else {
                    crate::numeric::unit_parameter(
                        t,
                        self.knots[j + k - D],
                        self.knots[j + 1 + k - r],
                    )
                };
                if alpha.is_nan() {
                    return Err(BSplineError::EvaluationAlphaIsNaN);
                }
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
    /// `accuracy` must be finite and positive, in coordinate units.
    /// Add the lengths of the spans. Gaps at jumps add no length.
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

    /// Return a spline with degree `D - 1`, `C - 1` controls and `K - 2` knots.
    /// Remove the first and last knots. Control `i` becomes
    /// `D * (P[i + 1] - P[i]) / (knots[i + D + 1] - knots[i + 1])`.
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

    /// Find a parameter where one coordinate equals `value`, using Newton's method.
    /// `eps` bounds the coordinate error and defaults to `1e-6`. It must be finite
    /// and positive. This local search can fail even when a solution exists.
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
    use super::*;
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
        let points = [
            ArrayPoint::new([0f64, 0f64]),
            ArrayPoint::new([10f64, 10f64]),
        ];
        let knots: [f64; 4] = [0., 0., 1., 1.];
        let curve: BSpline<ArrayPoint<f64, 2>, 4, 2, 1> = BSpline::new(knots, points).unwrap();

        assert!(curve.knot_kind == KnotVectorKind::Clamped);
        let (kmin, kmax) = curve.knot_domain();
        let nsteps: usize = 100;
        for t in 0..=nsteps {
            let t = kmin + (t as f64 / nsteps as f64) * (kmax - kmin);
            let _ = curve.eval(t).unwrap_or_else(|e| {
                panic!("Unexpected error in eval: {:?}", e);
            });
        }
        let start = curve.eval(0.).unwrap();
        let end = curve.eval(1.).unwrap();
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
        let points = [
            ArrayPoint::new([0f64, 0f64]),
            ArrayPoint::new([10f64, 10f64]),
        ];
        let knots: [f64; 4] = [0., 1., 2., 3.];
        let curve: BSpline<ArrayPoint<f64, 2>, 4, 2, 1> = BSpline::new(knots, points).unwrap();

        assert!(curve.knot_kind == KnotVectorKind::Unclamped);

        let (kmin, kmax) = curve.knot_domain();
        let nsteps: usize = 100;
        for t in 0..=nsteps {
            let t = kmin + (t as f64 / nsteps as f64) * (kmax - kmin);
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
        // The final value is 1. The last control point is 6.
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
        let points = [
            ArrayPoint::new([0f64, 0f64]),
            ArrayPoint::new([10f64, 10f64]), // Control point off the baseline
            ArrayPoint::new([20f64, 0f64]),
        ];
        let knots: [f64; 6] = [0., 0., 0., 1., 1., 1.];
        let curve: BSpline<ArrayPoint<f64, 2>, 6, 3, 2> = BSpline::new(knots, points).unwrap();

        assert!(curve.knot_kind == KnotVectorKind::Clamped);

        let (kmin, kmax) = curve.knot_domain();
        let nsteps: usize = 100;
        for t in 0..=nsteps {
            let t = kmin + (t as f64 / nsteps as f64) * (kmax - kmin);
            let _ = curve.eval(t).unwrap_or_else(|e| {
                panic!("Unexpected error in eval: {:?}", e);
            });
        }
        let start = curve.eval(0.).unwrap();
        let end = curve.eval(1.).unwrap();
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
        let points = [
            ArrayPoint::new([0f64, 0f64]),
            ArrayPoint::new([10f64, 10f64]),
            ArrayPoint::new([20f64, 0f64]), // End of the line at (20, 0)
        ];
        let knots: [f64; 6] = [0., 1., 2., 3., 4., 5.];
        let curve: BSpline<ArrayPoint<f64, 2>, 6, 3, 2> = BSpline::new(knots, points).unwrap();

        assert!(curve.knot_kind == KnotVectorKind::Unclamped);

        let (kmin, kmax) = curve.knot_domain();
        let nsteps: usize = 100;
        for t in 0..=nsteps {
            let t = kmin + (t as f64 / nsteps as f64) * (kmax - kmin);
            let _ = curve.eval(t).unwrap_or_else(|e| {
                panic!("Unexpected error in eval: {:?}", e);
            });
        }

        // Evaluate at the domain midpoint, t=2.5.
        let mid = curve.eval(2.5).unwrap();
        // For uniform unclamped degree-2 with 3 control points, weights are
        // [0.125, 0.75, 0.125] at t=2.5.
        let expected = points[0] * 0.125 + points[1] * 0.75 + points[2] * 0.125;
        assert_point_approx_equal(mid, expected, ABSOLUTE_TOLERANCE, RELATIVE_TOLERANCE);
    }

    #[test]
    fn degree_3_clamped_construct_and_eval_endpoints() {
        // A cubic with four controls requires eight knots.
        let points = [
            ArrayPoint::new([0f64, 1.77f64]),
            ArrayPoint::new([1.1f64, -1f64]),
            ArrayPoint::new([4.3f64, 3f64]),
            ArrayPoint::new([3.2f64, -4f64]),
            ArrayPoint::new([2.3f64, 1.77f64]),
        ];
        // Repeated endpoint knots force interpolation of the endpoint controls.
        let knots: [f64; 9] = [0., 0., 0., 0., 1., 2., 2., 2., 2.];
        let b: Result<BSpline<ArrayPoint<f64, 2>, 9, 5, 3>, BSplineError> =
            BSpline::new(knots, points);

        let curve = match b {
            Err(e) => panic!("Error: {:?}", e),
            Ok(b) => b,
        };

        assert!(curve.knot_kind == KnotVectorKind::Clamped);
        // Sample the active knot domain, including both endpoints.
        let (kmin, kmax) = curve.knot_domain();
        let nsteps: usize = 100;
        for t in 0..=nsteps {
            let t = kmin + (t as f64 / nsteps as f64) * (kmax - kmin);
            let _ = curve.eval(t).unwrap_or_else(|e| {
                panic!("Unexpected error in eval: {:?}", e);
            });
        }
        assert_eq!(
            curve.eval(-0.1).err(),
            Some(BSplineError::KnotDomainViolation)
        );
        assert_eq!(
            curve.eval(2.1).err(),
            Some(BSplineError::KnotDomainViolation)
        );
        assert!(curve.eval(kmin).is_ok());
        assert!(curve.eval(kmax).is_ok());
        assert_eq!(curve.eval(kmin).unwrap(), points[0]);
        assert_eq!(curve.eval(kmax).unwrap(), points[4]);
    }

    #[test]
    fn degree_3_unclamped_construct_and_eval() {
        // A cubic with four controls requires eight knots.
        let points = [
            ArrayPoint::new([0f64, 1.77f64]),
            ArrayPoint::new([1.1f64, -1f64]),
            ArrayPoint::new([4.3f64, 3f64]),
            ArrayPoint::new([3.2f64, -4f64]),
            ArrayPoint::new([2.3f64, 1.77f64]),
        ];
        // Unclamped knots do not force interpolation of the endpoint controls.
        let knots: [f64; 9] = [0., 1., 2., 3., 4., 5., 6., 7., 8.];
        let b: Result<BSpline<ArrayPoint<f64, 2>, 9, 5, 3>, BSplineError> =
            BSpline::new(knots, points);

        let curve = match b {
            Err(e) => panic!("Error: {:?}", e),
            Ok(b) => b,
        };

        assert!(curve.knot_kind == KnotVectorKind::Unclamped);
        // Sample the active knot domain, including both endpoints.
        let (kmin, kmax) = curve.knot_domain();
        assert!(kmin == 3.0 && kmax == 5.0);
        let nsteps: usize = 100;
        for t in 0..=nsteps {
            let t = kmin + (t as f64 / nsteps as f64) * (kmax - kmin);
            let _ = curve.eval(t).unwrap_or_else(|e| {
                panic!("Unexpected error in eval: {:?}", e);
            });
        }
        assert_eq!(
            curve.eval(-0.1).err(),
            Some(BSplineError::KnotDomainViolation)
        );
        assert_eq!(
            curve.eval(8.1).err(),
            Some(BSplineError::KnotDomainViolation)
        );
        // An unclamped curve need not interpolate its endpoint controls.
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
        // A quartic with five controls requires ten knots.
        let points = [
            ArrayPoint::new([0f64, 1.77f64]),
            ArrayPoint::new([1.1f64, -1f64]),
            ArrayPoint::new([4.3f64, 3f64]),
            ArrayPoint::new([3.2f64, -4f64]),
            ArrayPoint::new([2.3f64, 1.77f64]),
        ];
        // Repeated endpoint knots force interpolation of the endpoint controls.
        let knots: [f64; 10] = [0., 0., 0., 0., 0., 3., 3., 3., 3., 3.];
        let b: Result<BSpline<ArrayPoint<f64, 2>, 10, 5, 4>, BSplineError> =
            BSpline::new(knots, points);

        let curve = match b {
            Err(e) => panic!("Error: {:?}", e),
            Ok(b) => b,
        };

        assert!(curve.knot_kind == KnotVectorKind::Clamped);
        // Sample the active knot domain, including both endpoints.
        let (kmin, kmax) = curve.knot_domain();
        assert!(kmin == 0.0 && kmax == 3.0);
        let nsteps: usize = 100;
        for t in 0..=nsteps {
            let t = kmin + (t as f64 / nsteps as f64) * (kmax - kmin);
            let _ = curve.eval(t).unwrap_or_else(|e| {
                panic!("Unexpected error in eval: {:?}", e);
            });
        }
        assert_eq!(
            curve.eval(-0.1).err(),
            Some(BSplineError::KnotDomainViolation)
        );
        assert_eq!(
            curve.eval(3.1).err(),
            Some(BSplineError::KnotDomainViolation)
        );
        assert!(curve.eval(kmin).is_ok());
        assert!(curve.eval(kmax).is_ok());
        assert_eq!(curve.eval(kmin).unwrap(), points[0]);
        assert_eq!(curve.eval(kmax).unwrap(), points[4]);
    }

    #[test]
    fn degree_4_unclamped_construct_and_eval() {
        // A quartic with five controls requires ten knots.
        let points = [
            ArrayPoint::new([0f64, 1.77f64]),
            ArrayPoint::new([1.1f64, -1f64]),
            ArrayPoint::new([4.3f64, 3f64]),
            ArrayPoint::new([3.2f64, -4f64]),
            ArrayPoint::new([2.3f64, 1.77f64]),
            ArrayPoint::new([0.0, 0.0]),
        ];
        // Unclamped knots do not force interpolation of the endpoint controls.
        let knots: [f64; 11] = [0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10.];
        let b: Result<BSpline<ArrayPoint<f64, 2>, 11, 6, 4>, BSplineError> =
            BSpline::new(knots, points);

        let curve = match b {
            Err(e) => panic!("Error: {:?}", e),
            Ok(b) => b,
        };

        assert!(curve.knot_kind == KnotVectorKind::Unclamped);
        // Sample the active knot domain, including both endpoints.
        let (kmin, kmax) = curve.knot_domain();
        assert!(kmin == 4.0 && kmax == 6.0);
        let nsteps: usize = 100;
        for t in 0..=nsteps {
            let t = kmin + (t as f64 / nsteps as f64) * (kmax - kmin);
            let _ = curve.eval(t).unwrap_or_else(|e| {
                panic!("Unexpected error in eval: {:?}", e);
            });
        }
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
