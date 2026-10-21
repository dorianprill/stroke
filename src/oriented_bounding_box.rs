//! Approximate oriented bounding boxes for Bezier and B-spline curves.
//!
//! A box's frame is its set of perpendicular unit axes. The fitter searches
//! for a frame with small box volume, using curve subdivision and plane rotations.
//! It uses fixed stack storage and needs no PCA or eigensolver.
//!
//! # Method and references
//!
//! Both curve types use the same search:
//!
//! 1. Translate and uniformly scale the controls. Split the curve into Bezier
//!    pieces using de Casteljau's algorithm. Cache each piece's control polygon
//!    and endpoints. The piece lies in the convex hull of its controls.
//!    Projecting the endpoints gives a lower bound on curve width. Projecting
//!    the controls gives an upper bound in exact arithmetic.
//!    Split the pieces that contribute most to the uncertain box faces.
//! 2. Build starting frames from chords between endpoints. Choose further axes
//!    from the spread perpendicular to earlier axes. Keep several starting
//!    frames, including the coordinate axes and direction-aligned frames.
//! 3. Test individual and combined Givens rotations, each of which rotates two
//!    axes within their plane. Use smaller angular steps for thinner sides.
//!    Minimize padded log volume, with total side length breaking equal scores.
//! 4. Refine the best frames against the continuous curve using Bernstein bounds.
//!    Keep the best box after each trial. Larger budgets extend the same trial
//!    sequence. Model construction is fixed for budgets greater than one.
//!
//! B-splines contribute every nonzero Bezier span and both endpoints at jumps.
//! Their model capacity grows with the control count. Each trial combines all
//! span bounds in the same frame. Samples guide the choice of axes. The final
//! box uses bounds on the continuous curve.
//!
//! See [`Bezier::bounding_box_oriented_approx_with_budget`] and
//! [`BSpline::bounding_box_oriented_approx_with_budget`] for work limits,
//! stack usage and model-refinement targets.
//!
//! ## Related methods
//!
//! This library combines the following ideas. It does not reproduce any one
//! published algorithm or inherit its convergence and approximation guarantees.
//!
//! - The [`bezier` subdivision documentation](https://bezier.readthedocs.io/en/2024.6.20/algorithms/geometric-curve-intersection.html)
//!   explains control-hull containment and recursive subdivision.
//! - [ApproxMVBB](https://gabyx.github.io/ApproxMVBB/) builds a 3D frame using
//!   approximate diameters and perpendicular projection. Our chord-based starts
//!   use a related idea in any dimension. We do not use its projected
//!   minimum-rectangle or grid search.
//! - [Chan's 2018 paper, section 4](https://tmc.web.engr.illinois.edu/low_cheby.pdf),
//!   describes epsilon-kernels that preserve relative widths in every direction.
//!   Our model refines selected frames and has no such all-direction guarantee.
//! - [Shalit and Chechik, 2014](https://proceedings.mlr.press/v32/shalit14.html),
//!   use Givens rotations to optimize orthogonal matrices. Our derivative-free
//!   trials, combined rotations, angular scales and starting-frame selection
//!   differ from their coordinate-descent algorithm.
//!
//! A small model-volume gap measures uncertainty at a chosen frame. It does
//! not measure distance from the minimum-volume orientation. Neither this target
//! nor a true epsilon-kernel would certify that the rotation search found the
//! global minimum. Floating-point margins also remain separate from the
//! exact-arithmetic hull bounds. They are not certified interval arithmetic.

use crate::{ArrayPoint, BSpline, Bezier, BoundingBox, PointIndex};
use core::array;
use num_traits::{Float, NumCast};

mod search;

pub(crate) use search::{fit, fit_bspline};

type Frame<T, const D: usize> = [[T; D]; D];

/// Failure to construct an approximate oriented curve box.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OrientedBoundingBoxError {
    /// A box requires at least one coordinate dimension.
    ZeroDimension,
    /// At least one trial evaluation must be permitted.
    ZeroBudget,
    /// Normalization, projection, or numerical padding was not representable.
    NonFiniteComputation,
}

impl core::fmt::Display for OrientedBoundingBoxError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.write_str(match self {
            Self::ZeroDimension => "an oriented box needs at least one dimension",
            Self::ZeroBudget => "the evaluation budget must be positive",
            Self::NonFiniteComputation => "oriented-box arithmetic was non-finite",
        })
    }
}

impl core::error::Error for OrientedBoundingBoxError {}

#[derive(Debug, Clone, Copy, PartialEq)]
struct Candidate<T, const D: usize> {
    axes: Frame<T, D>,
    bounds: BoundingBox<T, D>,
    // Log volume avoids overflow/underflow. Total width breaks equal scores.
    score: (T, T),
}

/// A padded box whose axes are perpendicular unit vectors to floating-point precision.
///
/// Returned by [`Bezier::bounding_box_oriented_approx`] and
/// [`BSpline::bounding_box_oriented_approx`]. The search may stop before reaching
/// even a local minimum. There is no guarantee on how close its volume is to
/// the global minimum.
///
/// Bernstein subdivision bounds the continuous curve. Floating-point margins
/// account for numerical error but do not provide certified interval bounds.
/// Translation and uniform scaling limit intermediate overflow and underflow.
///
/// Axis directions, handedness and exact fitted values may change between
/// library versions. Constructed boxes require `D > 0`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct OrientedBoundingBox<T, const D: usize> {
    origin: [T; D],
    scale: T,
    // Local intervals are normalized: world = origin + scale * axes^T * local.
    candidate: Candidate<T, D>,
    evaluations: usize,
    diagnostics: Option<OrientedBoundingBoxDiagnostics<T>>,
}

/// Work counts and model-error estimates for an oriented box.
///
/// The model consists of subdivided control polygons and their endpoints.
/// Error estimates apply at the returned orientation. They do not certify
/// floating-point error or closeness to the minimum-volume box.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct OrientedBoundingBoxDiagnostics<T> {
    /// Number of subdivision endpoints retained.
    pub representative_points: usize,
    /// Whether subdivision reached its capacity or parameter resolution.
    /// Capacity is 64 leaves for Beziers, or `C + 64` for B-splines.
    pub representative_limit_reached: bool,
    /// Number of frames scored against the subdivided control hulls.
    /// Subdivision-error checks are not included.
    pub point_set_evaluations: usize,
    /// Natural log of the padded model's upper/lower volume ratio at the returned frame.
    pub model_log_volume_gap: T,
    /// Maximum control-hull distance from a representative chord, in world units.
    /// An unrepresentable world-space error can overflow to infinity.
    pub representative_error: T,
}

fn scalar<T: Float>(value: f64) -> T {
    <T as NumCast>::from(value).unwrap()
}

fn identity<T: Float, const D: usize>() -> Frame<T, D> {
    array::from_fn(|i| array::from_fn(|j| if i == j { T::one() } else { T::zero() }))
}

fn dot<T: Float, const D: usize>(a: &[T; D], b: &[T; D]) -> T {
    a.iter()
        .zip(b)
        .fold(T::zero(), |sum, (&a, &b)| a.mul_add(b, sum))
}

fn normalize<T: Float, const D: usize>(mut v: [T; D]) -> Option<[T; D]> {
    let scale = v.iter().fold(T::zero(), |s, x| s.max(x.abs()));
    if scale == T::zero() || !scale.is_finite() {
        return None;
    }
    v = v.map(|x| x / scale);
    let norm = dot(&v, &v).sqrt();
    Some(v.map(|x| x / norm))
}

// Reorthogonalize twice so repeated rotation trials cannot accumulate drift.
fn orthogonalize<T: Float, const D: usize>(mut axes: Frame<T, D>) -> Option<Frame<T, D>> {
    for i in 0..D {
        let mut axis = axes[i];
        for _ in 0..2 {
            for previous in &axes[..i] {
                let projection = dot(&axis, previous);
                for (x, &previous) in axis.iter_mut().zip(previous) {
                    *x = (-projection).mul_add(previous, *x);
                }
            }
        }
        axes[i] = normalize(axis)?;
    }
    Some(axes)
}

// A Householder reflection supplies a complete frame from one direction.
// Either handedness is valid for a box. The sign choice avoids cancellation.
fn aligned_frame<T: Float, const D: usize>(direction: [T; D]) -> Option<Frame<T, D>> {
    let mut v = normalize(direction)?;
    v[0] = v[0] + T::one().copysign(v[0]);
    let factor = scalar::<T>(2.0) / dot(&v, &v);
    Some(array::from_fn(|i| {
        array::from_fn(|j| {
            let entry = if i == j { T::one() } else { T::zero() };
            (-factor * v[i]).mul_add(v[j], entry)
        })
    }))
}

fn rotated<T: Float, const D: usize>(
    axes: Frame<T, D>,
    i: usize,
    j: usize,
    angle: T,
) -> Frame<T, D> {
    let (s, c) = angle.sin_cos();
    let mut result = axes;
    for k in 0..D {
        result[i][k] = c.mul_add(axes[i][k], s * axes[j][k]);
        result[j][k] = (-s).mul_add(axes[i][k], c * axes[j][k]);
    }
    result
}

// Project each coordinate onto the supplied box axes and bound the resulting
// scalar polynomials. The returned intervals are normalized and unpadded.
fn projected_bounds<T: Float, const D: usize, const N: usize>(
    points: &[[T; D]; N],
    axes: &Frame<T, D>,
) -> Result<[(T, T); D], OrientedBoundingBoxError> {
    let mut bounds = [(T::zero(), T::zero()); D];
    for (axis, bounds) in axes.iter().zip(&mut bounds) {
        let projected = points.map(|p| dot(axis, &p));
        if projected.iter().any(|x| !x.is_finite()) {
            return Err(OrientedBoundingBoxError::NonFiniteComputation);
        }
        *bounds = Bezier::<ArrayPoint<T, 1>, N>::bound_scalar_axis(projected);
    }
    Ok(bounds)
}

fn evaluate<T: Float, const D: usize>(
    roundoff: &[T; D],
    margin: T,
    axes: Frame<T, D>,
    bounds_for_frame: &mut impl FnMut(&Frame<T, D>) -> Result<[(T, T); D], OrientedBoundingBoxError>,
) -> Result<Candidate<T, D>, OrientedBoundingBoxError> {
    let axes = orthogonalize(axes).ok_or(OrientedBoundingBoxError::NonFiniteComputation)?;
    let intervals = bounds_for_frame(&axes)?;
    let mut lo = [T::zero(); D];
    let mut hi = [T::zero(); D];
    for (i, axis) in axes.iter().enumerate() {
        let padding = axis
            .iter()
            .zip(roundoff)
            .fold(margin, |sum, (&a, &error)| a.abs().mul_add(error, sum));
        lo[i] = intervals[i].0 - padding;
        hi[i] = intervals[i].1 + padding;
    }
    let bounds = BoundingBox::new(lo, hi).ok_or(OrientedBoundingBoxError::NonFiniteComputation)?;
    let score = bounds
        .extent()
        .iter()
        .fold((T::zero(), T::zero()), |(volume, sum), &w| {
            (volume + w.ln(), sum + w)
        });
    if !score.0.is_finite() || !score.1.is_finite() {
        return Err(OrientedBoundingBoxError::NonFiniteComputation);
    }
    Ok(Candidate {
        axes,
        bounds,
        score,
    })
}

fn improves<T: Float, const D: usize>(trial: &Candidate<T, D>, best: &Candidate<T, D>) -> bool {
    trial.score.0 < best.score.0 || (trial.score.0 == best.score.0 && trial.score.1 < best.score.1)
}

struct Normalized<T, const D: usize, const C: usize> {
    origin: [T; D],
    scale: T,
    points: [[T; D]; C],
    roundoff: [T; D],
    margin: T,
}

fn normalized_controls<T: Float, const D: usize, const C: usize>(
    original: [[T; D]; C],
    origin: [T; D],
    arithmetic_work: T,
) -> Result<Normalized<T, D, C>, OrientedBoundingBoxError> {
    let differences = original.map(|p| array::from_fn(|i| p[i] - origin[i]));
    let scale = differences
        .iter()
        .flatten()
        .fold(T::zero(), |s, x| s.max(x.abs()));
    if !scale.is_finite() || differences.iter().flatten().any(|x| !x.is_finite()) {
        return Err(OrientedBoundingBoxError::NonFiniteComputation);
    }
    let points = if scale == T::zero() {
        differences
    } else {
        differences.map(|p| p.map(|x| x / scale))
    };
    let roundoff = array::from_fn(|i| {
        if scale == T::zero() {
            return T::zero();
        }
        let magnitude = original.iter().fold(T::zero(), |m, p| m.max(p[i].abs()));
        let factor = scalar::<T>(8.0) * arithmetic_work;
        let absolute =
            magnitude * T::epsilon() * factor + (T::min_positive_value() * T::epsilon()) * factor;
        absolute / scale
    });
    // Include de Boor span-extraction work in the B-spline margin.
    // This allowance estimates roundoff. It is not a certified interval bound.
    let margin = scalar::<T>(64.0)
        * T::epsilon()
        * (<T as NumCast>::from(D).unwrap() + arithmetic_work + T::one());
    Ok(Normalized {
        origin,
        scale,
        points,
        roundoff,
        margin,
    })
}

fn check_request(dimension: usize, budget: usize) -> Result<(), OrientedBoundingBoxError> {
    if dimension == 0 {
        return Err(OrientedBoundingBoxError::ZeroDimension);
    }
    if budget == 0 {
        return Err(OrientedBoundingBoxError::ZeroBudget);
    }
    Ok(())
}

// Constants and budget-one/1D requests need no representative model.
fn initial_box<T: Float, const D: usize, const C: usize>(
    data: &Normalized<T, D, C>,
    bounds_for_frame: &mut impl FnMut(&Frame<T, D>) -> Result<[(T, T); D], OrientedBoundingBoxError>,
) -> Result<OrientedBoundingBox<T, D>, OrientedBoundingBoxError> {
    let candidate = if data.scale == T::zero() {
        Candidate {
            axes: identity(),
            bounds: BoundingBox::new([T::zero(); D], [T::zero(); D]).unwrap(),
            score: (T::neg_infinity(), T::zero()),
        }
    } else {
        evaluate(&data.roundoff, data.margin, identity(), bounds_for_frame)?
    };
    Ok(OrientedBoundingBox {
        origin: data.origin,
        scale: data.scale,
        candidate,
        evaluations: if data.scale == T::zero() { 0 } else { 1 },
        diagnostics: None,
    })
}

#[allow(clippy::type_complexity)]
fn normalized_bspline<
    T: Float,
    P: PointIndex<Scalar = T>,
    const K: usize,
    const C: usize,
    const D: usize,
>(
    curve: &BSpline<P, K, C, D>,
) -> Result<
    (
        Normalized<T, { P::DIM }, C>,
        BSpline<ArrayPoint<T, { P::DIM }>, K, C, D>,
    ),
    OrientedBoundingBoxError,
> {
    let original = array::from_fn::<_, C, _>(|i| {
        array::from_fn(|axis| curve.control_points().as_slice()[i][axis])
    });
    let order = <T as NumCast>::from(D + 1).unwrap();
    let data = normalized_controls(original, curve.control_box().center(), order * order)?;
    let normalized: BSpline<ArrayPoint<T, { P::DIM }>, K, C, D> = BSpline::new(
        array::from_fn(|i| curve.knots().as_slice()[i]),
        data.points.map(ArrayPoint::from),
    )
    .map_err(|_| OrientedBoundingBoxError::NonFiniteComputation)?;
    Ok((data, normalized))
}

fn projected_bspline_bounds<
    T: Float,
    const DIM: usize,
    const K: usize,
    const C: usize,
    const D: usize,
>(
    curve: &BSpline<ArrayPoint<T, DIM>, K, C, D>,
    axes: &Frame<T, DIM>,
) -> Result<[(T, T); DIM], OrientedBoundingBoxError>
where
    [(); D + 1]: Sized,
{
    let mut bounds = [(T::infinity(), T::neg_infinity()); DIM];
    for (_, span) in curve.bezier_spans() {
        let points = span.into_control_points().map(ArrayPoint::into_array);
        let projected = projected_bounds(&points, axes)?;
        for (total, (lo, hi)) in bounds.iter_mut().zip(projected) {
            total.0 = total.0.min(lo);
            total.1 = total.1.max(hi);
        }
    }
    Ok(bounds)
}

impl<T: Float, const D: usize> OrientedBoundingBox<T, D> {
    /// Return model diagnostics, or `None` when the call skipped preprocessing.
    /// Constant curves, one-dimensional curves and budget-one calls skip it.
    pub fn diagnostics(&self) -> Option<OrientedBoundingBoxDiagnostics<T>> {
        self.diagnostics
    }

    /// Orthonormal box axes in world coordinates, stored as rows.
    ///
    /// Orthogonality holds to floating-point precision; handedness is unspecified.
    pub fn axes(&self) -> &[[T; D]; D] {
        &self.candidate.axes
    }

    /// Center in world coordinates. An unrepresentable coordinate can be infinite.
    pub fn center(&self) -> [T; D] {
        self.world_from_normalized(self.candidate.bounds.center())
    }

    /// Half the side lengths along the box axes.
    ///
    /// Very large values can overflow to infinity; very small values can round
    /// to zero. The internal normalized representation retains those intervals.
    pub fn half_extents(&self) -> [T; D] {
        let half = scalar::<T>(0.5);
        self.candidate
            .bounds
            .axes()
            .map(|(lo, hi)| (hi * half - lo * half) * self.scale)
    }

    /// Return one corner. `upper[i]` chooses the upper end along box axis `i`.
    ///
    /// Returns world coordinates. A coordinate outside the scalar range can
    /// overflow to infinity. Constructing one corner requires no allocation.
    pub fn corner(&self, upper: [bool; D]) -> [T; D] {
        self.world_from_normalized(array::from_fn(|i| {
            let (lo, hi) = self.candidate.bounds.axis(i);
            if upper[i] { hi } else { lo }
        }))
    }

    /// Return the number of frames evaluated against the continuous curve.
    /// Includes the initial axis-aligned frame. Excludes model-search trials
    /// and subdivision work within a trial. Constant curves use zero trials.
    pub fn evaluations(&self) -> usize {
        self.evaluations
    }

    /// Test membership using the stored frame and padded intervals.
    ///
    /// A roundoff allowance can accept points just outside a face.
    /// Non-finite points and dimension mismatches return `false`.
    /// This is a floating-point test, not an exact predicate.
    pub fn contains<P: PointIndex<Scalar = T>>(&self, point: &P) -> bool {
        if P::DIM != D || !point.is_finite() {
            return false;
        }
        if self.scale == T::zero() {
            return self.origin.iter().enumerate().all(|(i, &x)| point[i] == x);
        }
        let local = array::from_fn(|i| (point[i] - self.origin[i]) / self.scale);
        if local.iter().any(|x| !x.is_finite()) {
            return false;
        }
        // Converting a corner to world coordinates and projecting it back
        // introduces roundoff. Include that error in containment checks while
        // keeping the volume used by the search unchanged.
        let epsilon = T::epsilon();
        let quantum = T::min_positive_value() * epsilon;
        let input_error: [T; D] = array::from_fn(|i| {
            (point[i].abs() * epsilon + self.origin[i].abs() * epsilon + quantum) / self.scale
        });
        let factor = scalar::<T>(32.0) * (<T as NumCast>::from(D).unwrap() + T::one()) * epsilon;
        self.candidate
            .axes
            .iter()
            .zip(self.candidate.bounds.axes())
            .all(|(axis, &(lo, hi))| {
                let coordinate = dot(axis, &local);
                let magnitude = axis
                    .iter()
                    .zip(local)
                    .fold(lo.abs().max(hi.abs()), |sum, (&a, x)| {
                        a.abs().mul_add(x.abs(), sum)
                    });
                let error = axis
                    .iter()
                    .zip(input_error)
                    .fold(factor * magnitude, |sum, (&a, error)| {
                        a.abs().mul_add(error, sum)
                    });
                error.is_finite() && coordinate >= lo - error && coordinate <= hi + error
            })
    }

    fn world_from_normalized(&self, local: [T; D]) -> [T; D] {
        array::from_fn(|i| {
            let coordinate = (0..D).fold(T::zero(), |sum, j| {
                self.candidate.axes[j][i].mul_add(local[j], sum)
            });
            self.scale.mul_add(coordinate, self.origin[i])
        })
    }

    /// Natural logarithm of the padded area, volume, or D-dimensional measure.
    ///
    /// Returns negative infinity for a constant curve. Computing in logarithmic
    /// units avoids intermediate products overflowing or underflowing.
    pub fn log_volume(&self) -> T {
        self.candidate.score.0 + <T as NumCast>::from(D).unwrap() * self.scale.ln()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Find cubic extrema independently through power-basis derivative roots.
    // The fixtures use modest coordinates so the reference calculation is reliable.
    pub(super) fn cubic_range(p: [f64; 4]) -> (f64, f64) {
        let mut lo = p[0].min(p[3]);
        let mut hi = p[0].max(p[3]);
        let a = 3.0 * (-p[0] + 3.0 * p[1] - 3.0 * p[2] + p[3]);
        let b = 6.0 * (p[0] - 2.0 * p[1] + p[2]);
        let c = 3.0 * (p[1] - p[0]);
        let mut visit = |t: f64| {
            if t > 0.0 && t < 1.0 {
                let s = 1.0 - t;
                let value = s * s * s * p[0]
                    + 3.0 * s * s * t * p[1]
                    + 3.0 * s * t * t * p[2]
                    + t * t * t * p[3];
                lo = lo.min(value);
                hi = hi.max(value);
            }
        };
        if a == 0.0 {
            if b != 0.0 {
                visit(-c / b);
            }
        } else {
            let discriminant = b * b - 4.0 * a * c;
            if discriminant >= 0.0 {
                let q = -0.5 * (b + discriminant.sqrt().copysign(b));
                if q == 0.0 {
                    visit(-b / (2.0 * a));
                } else {
                    visit(q / a);
                    visit(c / q);
                }
            }
        }
        (lo, hi)
    }

    fn cubic_log_volume<const D: usize>(controls: &[[f64; D]; 4], frame: Frame<f64, D>) -> f64 {
        frame
            .iter()
            .map(|axis| {
                let (lo, hi) = cubic_range(controls.map(|p| dot(axis, &p)));
                (hi - lo).ln()
            })
            .sum()
    }

    fn grid_reference_2d(controls: &[[f64; 2]; 4]) -> f64 {
        let spacing = core::f64::consts::FRAC_PI_2 / 4096.0;
        let objective = |angle| cubic_log_volume(controls, rotated(identity(), 0, 1, angle));
        let mut best_angle = 0.0;
        let mut best = objective(best_angle);
        for i in 1..4096 {
            let angle = i as f64 * spacing;
            let volume = objective(angle);
            if volume < best {
                best = volume;
                best_angle = angle;
            }
        }
        // Refine the best grid cell. Thin curves can have a narrow minimum
        // between sampled angles.
        let mut lo = best_angle - spacing;
        let mut hi = best_angle + spacing;
        let fraction = (5.0_f64.sqrt() - 1.0) / 2.0;
        let mut left = hi - fraction * (hi - lo);
        let mut right = lo + fraction * (hi - lo);
        let mut a = objective(left);
        let mut b = objective(right);
        for _ in 0..64 {
            if a < b {
                hi = right;
                right = left;
                b = a;
                left = hi - fraction * (hi - lo);
                a = objective(left);
            } else {
                lo = left;
                left = right;
                a = b;
                right = lo + fraction * (hi - lo);
                b = objective(right);
            }
            best = best.min(a).min(b);
        }
        best
    }

    fn fixtures_2d() -> [(&'static str, [[f64; 2]; 4]); 6] {
        [
            (
                "Rotated arch",
                [[0.0, 0.0], [1.0, 3.0], [3.0, 3.0], [4.0, 0.0]],
                0.63,
            ),
            (
                "S bend",
                [[0.0, 0.0], [3.0, 4.0], [2.0, -3.0], [6.0, 1.0]],
                0.4,
            ),
            (
                "Thin diagonal",
                [[-3.0, -0.02], [-1.0, 0.03], [1.0, -0.03], [3.0, 0.02]],
                0.78,
            ),
            (
                "Closed loop",
                [[0.0, 0.0], [3.0, 4.0], [-3.0, 4.0], [0.0, 0.0]],
                0.37,
            ),
            (
                "Backtracking",
                [[0.0, 0.0], [4.0, 3.0], [-3.0, -2.0], [1.0, 1.0]],
                0.8,
            ),
            (
                "Skew cubic",
                [[0.0, 1.77], [1.1, -1.0], [4.3, 3.0], [3.2, -4.0]],
                0.0,
            ),
        ]
        .map(|(name, points, angle)| {
            let frame = rotated(identity(), 0, 1, angle);
            (name, points.map(|p| frame.map(|axis| dot(&axis, &p))))
        })
    }

    fn random(state: &mut u64) -> f64 {
        *state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((*state >> 11) as f64 / ((1_u64 << 53) as f64)) * 2.0 - 1.0
    }

    fn validate<T: Float, P: PointIndex<Scalar = T>, const N: usize>(
        curve: &Bezier<P, N>,
        budget: usize,
    ) where
        [(); P::DIM]: Sized,
    {
        let fit = curve
            .bounding_box_oriented_approx_with_budget(budget)
            .unwrap();
        assert!(fit.evaluations <= budget);
        for i in 0..=1024 {
            let t = <T as NumCast>::from(i).unwrap() / scalar::<T>(1024.0);
            assert!(fit.contains(&curve.eval(t)), "sample {i}/1024 escaped");
        }
        let axes = fit.candidate.axes;
        for i in 0..P::DIM {
            for j in 0..P::DIM {
                let expected = if i == j { scalar(1.0) } else { scalar(0.0) };
                assert!(
                    (dot(&axes[i], &axes[j]) - expected).abs()
                        <= scalar::<T>(16.0 * P::DIM as f64) * T::epsilon()
                );
            }
        }
        if fit.scale > T::zero() {
            let points: [_; N] = array::from_fn(|i| {
                array::from_fn(|axis| {
                    (curve.control_points().as_slice()[i][axis] - fit.origin[axis]) / fit.scale
                })
            });
            let original = array::from_fn::<_, N, _>(|i| {
                array::from_fn(|axis| curve.control_points().as_slice()[i][axis])
            });
            let data = normalized_controls(original, fit.origin, <T as NumCast>::from(N).unwrap())
                .unwrap();
            let baseline = evaluate(&data.roundoff, data.margin, identity(), &mut |axes| {
                projected_bounds(&points, axes)
            })
            .unwrap();
            assert!(!improves(&baseline, &fit.candidate));
        }
    }

    // Check every fitted face against polynomial extrema. Modest fixture scales
    // keep this independent reference calculation reliable.
    fn check_cubic_faces<const D: usize>(controls: [[f64; D]; 4])
    where
        [(); <ArrayPoint<f64, D> as crate::Point>::DIM]: Sized,
    {
        let curve = Bezier::new(controls.map(ArrayPoint::from)).unwrap();
        let fit = curve.bounding_box_oriented_approx_with_budget(256).unwrap();
        if fit.scale == 0.0 {
            return;
        }
        for (axis, &(lo, hi)) in fit.candidate.axes.iter().zip(fit.candidate.bounds.axes()) {
            let projected = controls.map(|p| {
                (0..D)
                    .map(|j| axis[j] * ((p[j] - fit.origin[j]) / fit.scale))
                    .sum()
            });
            let (minimum, maximum) = cubic_range(projected);
            assert!(minimum >= lo, "minimum {minimum} < {lo}");
            assert!(maximum <= hi, "maximum {maximum} > {hi}");
        }
        for mask in 0..(1_usize << D) {
            let corner = fit.corner(array::from_fn(|i| mask & (1 << i) != 0));
            assert!(
                fit.contains(&ArrayPoint::from(corner)),
                "box corner escaped"
            );
        }
        validate(&curve, 256);
    }

    #[test]
    fn independent_cubic_extrema_enclose_all_faces_2d_3d_4d() {
        let mut state = 0x5eed;
        for _ in 0..64 {
            check_cubic_faces(array::from_fn::<_, 4, _>(|_| {
                array::from_fn::<_, 2, _>(|_| 4.0 * random(&mut state))
            }));
            check_cubic_faces(array::from_fn::<_, 4, _>(|_| {
                array::from_fn::<_, 3, _>(|_| 4.0 * random(&mut state))
            }));
            check_cubic_faces(array::from_fn::<_, 4, _>(|_| {
                array::from_fn::<_, 4, _>(|_| 4.0 * random(&mut state))
            }));
        }
    }

    fn check_spline_faces<const DIM: usize>(curve: &BSpline<ArrayPoint<f64, DIM>, 11, 7, 3>)
    where
        [(); <ArrayPoint<f64, DIM> as crate::Point>::DIM]: Sized,
    {
        let fit = curve.bounding_box_oriented_approx_with_budget(64).unwrap();
        for (_, span) in curve.bezier_spans() {
            for (axis, &(lo, hi)) in fit.candidate.axes.iter().zip(fit.candidate.bounds.axes()) {
                let projected = span.into_control_points().map(|p| {
                    (0..DIM)
                        .map(|i| axis[i] * ((p[i] - fit.origin[i]) / fit.scale))
                        .sum()
                });
                let (minimum, maximum) = cubic_range(projected);
                assert!(minimum >= lo, "spline minimum {minimum} < {lo}");
                assert!(maximum <= hi, "spline maximum {maximum} > {hi}");
            }
        }
    }

    #[test]
    fn independent_cubic_extrema_enclose_every_spline_span() {
        let mut state = 0xb5_11;
        for _ in 0..8 {
            let controls2 = array::from_fn::<_, 7, _>(|_| {
                ArrayPoint::from(array::from_fn::<_, 2, _>(|_| 4.0 * random(&mut state)))
            });
            let controls3 = array::from_fn::<_, 7, _>(|_| {
                ArrayPoint::from(array::from_fn::<_, 3, _>(|_| 4.0 * random(&mut state)))
            });
            check_spline_faces(&BSpline::<_, 11, 7, 3>::new_clamped_uniform(controls2).unwrap());
            check_spline_faces(&BSpline::<_, 11, 7, 3>::new_clamped_uniform(controls3).unwrap());
            let knots = array::from_fn(|i| i as f64);
            check_spline_faces(&BSpline::<_, 11, 7, 3>::new(knots, controls2).unwrap());
            check_spline_faces(&BSpline::<_, 11, 7, 3>::new(knots, controls3).unwrap());
        }
    }

    #[test]
    fn zero_dimension_and_zero_budget_are_rejected() {
        let empty = Bezier::new([ArrayPoint::<f64, 0>::from([])]).unwrap();
        assert!(matches!(
            empty.bounding_box_oriented_approx_with_budget(1),
            Err(OrientedBoundingBoxError::ZeroDimension)
        ));
        let curve = Bezier::new([ArrayPoint::from([0.0_f64])]).unwrap();
        assert!(matches!(
            curve.bounding_box_oriented_approx_with_budget(0),
            Err(OrientedBoundingBoxError::ZeroBudget)
        ));
    }

    #[test]
    fn budgets_are_hard_limits_and_results_are_deterministic() {
        let curve = Bezier::new(fixtures_2d()[1].1.map(ArrayPoint::from)).unwrap();
        for budget in [1, 2, 3, 8, 17, 64, 256, 1024] {
            validate(&curve, budget);
            let a = curve
                .bounding_box_oriented_approx_with_budget(budget)
                .unwrap();
            let b = curve
                .bounding_box_oriented_approx_with_budget(budget)
                .unwrap();
            assert_eq!(a.candidate.axes, b.candidate.axes);
            assert_eq!(a.candidate.bounds, b.candidate.bounds);
            assert_eq!(a.evaluations, b.evaluations);
        }
    }

    #[test]
    fn constants_lines_backtracking_and_closed_curves() {
        validate(
            &Bezier::new([ArrayPoint::from([2.0_f64, -3.0, 4.0])]).unwrap(),
            64,
        );
        validate(
            &Bezier::new([ArrayPoint::from([-2.0_f64]), ArrayPoint::from([3.0])]).unwrap(),
            64,
        );
        validate(
            &Bezier::new([
                ArrayPoint::from([0.0_f64, 0.0, 0.0]),
                ArrayPoint::from([4.0, 8.0, 12.0]),
            ])
            .unwrap(),
            256,
        );
        for (_, controls) in fixtures_2d() {
            check_cubic_faces(controls);
        }
        let curve = Bezier::new([ArrayPoint::from([2.0_f64, -3.0])]).unwrap();
        let fit = curve.bounding_box_oriented_approx_with_budget(1).unwrap();
        assert_eq!(fit.evaluations, 0);
        assert!(!fit.contains(&ArrayPoint::from([2.0, -2.0])));
        assert!(!fit.contains(&ArrayPoint::from([2.0, f64::NAN])));
        assert!(!fit.contains(&ArrayPoint::from([2.0])));
    }

    #[test]
    fn generic_degrees_and_dimensions() {
        validate(
            &Bezier::new([[0.0_f64, 0.0], [2.0, 4.0], [4.0, 0.0]].map(ArrayPoint::from)).unwrap(),
            256,
        );
        validate(
            &Bezier::new(
                [
                    [0.0_f64, 0.0, 0.0, 0.0],
                    [1.0, 3.0, 2.0, 4.0],
                    [4.0, -2.0, 1.0, -1.0],
                    [5.0, 1.0, 4.0, 2.0],
                    [3.0, -1.0, 3.0, 1.0],
                ]
                .map(ArrayPoint::from),
            )
            .unwrap(),
            256,
        );
        let mut state = 123;
        validate(
            &Bezier::new(array::from_fn::<_, 9, _>(|_| {
                ArrayPoint::from(array::from_fn::<_, 6, _>(|_| random(&mut state)))
            }))
            .unwrap(),
            256,
        );
    }

    #[test]
    fn scale_and_translation_do_not_break_containment() {
        for scale in [1e-200, 1e-100, 1.0, 1e100, 1e200] {
            let curve = Bezier::new(
                fixtures_2d()[0]
                    .1
                    .map(|p| ArrayPoint::from(p.map(|x| x * scale))),
            )
            .unwrap();
            validate(&curve, 256);
        }
        for offset in [1e3, 1e9, 1e12] {
            let curve = Bezier::new(
                fixtures_2d()[0]
                    .1
                    .map(|p| ArrayPoint::from(p.map(|x| x + offset))),
            )
            .unwrap();
            validate(&curve, 256);
        }
        validate(
            &Bezier::new([[-f64::MAX, 0.0], [f64::MAX, 1.0]].map(ArrayPoint::from)).unwrap(),
            64,
        );
        validate(
            &Bezier::new(
                [[0.0_f64, 0.0], [f64::from_bits(1), f64::from_bits(2)]].map(ArrayPoint::from),
            )
            .unwrap(),
            64,
        );
    }

    #[test]
    fn f32_supports_different_scales_and_dimensions() {
        for scale in [1e-20_f32, 1.0, 1e20] {
            validate(
                &Bezier::new(
                    fixtures_2d()[1]
                        .1
                        .map(|p| ArrayPoint::from(p.map(|x| x as f32 * scale))),
                )
                .unwrap(),
                256,
            );
            validate(
                &Bezier::new(
                    [
                        [0.0_f32, 0.0, 0.0],
                        [1.0, 3.0, 2.0],
                        [4.0, -2.0, 1.0],
                        [5.0, 1.0, 4.0],
                    ]
                    .map(|p| ArrayPoint::from(p.map(|x| x * scale))),
                )
                .unwrap(),
                256,
            );
        }
        let curve = Bezier::new(
            fixtures_2d()[1]
                .1
                .map(|p| ArrayPoint::from(p.map(|x| x as f32 + 1e6))),
        )
        .unwrap();
        validate(&curve, 256);
        validate(
            &Bezier::new(
                [[0.0_f32, 0.0], [f32::from_bits(1), f32::from_bits(2)]].map(ArrayPoint::from),
            )
            .unwrap(),
            64,
        );
    }

    #[test]
    fn unrepresentable_normalized_roundoff_returns_error() {
        let curve =
            Bezier::new([[f64::MAX, 0.0], [f64::MAX, f64::from_bits(1)]].map(ArrayPoint::from))
                .unwrap();
        assert!(matches!(
            curve.bounding_box_oriented_approx_with_budget(64),
            Err(OrientedBoundingBoxError::NonFiniteComputation)
        ));
    }

    #[test]
    fn two_dimensional_quality_is_close_to_dense_independent_reference() {
        for (name, controls) in fixtures_2d() {
            let curve = Bezier::new(controls.map(ArrayPoint::from)).unwrap();
            let fit = curve
                .bounding_box_oriented_approx_with_budget(1024)
                .unwrap();
            let ratio = (fit.log_volume() - grid_reference_2d(&controls)).exp();
            assert!(ratio <= 1.01, "{name}: OBB/grid = {ratio}");
        }
    }

    #[cfg(feature = "nalgebra")]
    #[test]
    fn nalgebra_uses_the_same_point_index_contract() {
        let curve = Bezier::new([
            nalgebra::SVector::<f64, 3>::new(0.0, 0.0, 0.0),
            nalgebra::SVector::<f64, 3>::new(1.0, 3.0, 2.0),
            nalgebra::SVector::<f64, 3>::new(4.0, -2.0, 1.0),
            nalgebra::SVector::<f64, 3>::new(5.0, 1.0, 4.0),
        ])
        .unwrap();
        validate(&curve, 256);
    }
}
