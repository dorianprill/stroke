//! B-spline path utilities.

use core::slice;

use super::{ArrayVec, BSpline, Point, PointDot, PointIndex, PointNorm};
use num_traits::{Float, NumCast};

const DEFAULT_LENGTH_STEPS: usize = 64;
const MAX_LENGTH_ITERS: usize = 32;

/// A path composed of multiple B-spline segments of the same degree/knots layout.
///
/// Component-aware helpers (e.g. bounding boxes) require `PointIndex`.
///
/// # Examples
/// ```rust
/// use stroke::{BSpline, BSplinePath, PointN};
///
/// let knots: [f64; 4] = [0.0, 0.0, 1.0, 1.0];
/// let seg1: BSpline<PointN<f64, 2>, 4, 2, 1> = BSpline::new(
///     knots,
///     [
///         PointN::new([0.0, 0.0]),
///         PointN::new([1.0, 0.0]),
///     ],
/// )
/// .unwrap();
/// let seg2: BSpline<PointN<f64, 2>, 4, 2, 1> = BSpline::new(
///     knots,
///     [
///         PointN::new([1.0, 0.0]),
///         PointN::new([2.0, 0.0]),
///     ],
/// )
/// .unwrap();
///
/// let mut path: BSplinePath<PointN<f64, 2>, 4, 2, 1, 2> = BSplinePath::new();
/// path.push(seg1);
/// path.push(seg2);
///
/// assert_eq!(path.start().unwrap(), PointN::new([0.0, 0.0]));
/// assert_eq!(path.end().unwrap(), PointN::new([2.0, 0.0]));
///
/// let t = path.t_at_length(1.0).unwrap();
/// let p = path.point_at_length(1.0).unwrap();
/// let tangent = path.tangent(0.25).unwrap();
/// # let _ = (t, p, tangent);
/// ```
pub struct BSplinePath<P, const K: usize, const C: usize, const D: usize, const N: usize>
where
    P: Point,
{
    segments: ArrayVec<[BSpline<P, K, C, D>; N]>,
}

impl<P, const K: usize, const C: usize, const D: usize, const N: usize> BSplinePath<P, K, C, D, N>
where
    P: Point,
    [BSpline<P, K, C, D>; N]: tinyvec::Array<Item = BSpline<P, K, C, D>>,
{
    /// Create an empty path with capacity `N`.
    pub fn new() -> Self {
        BSplinePath {
            segments: ArrayVec::new(),
        }
    }

    /// Return the number of segments currently stored.
    pub fn len(&self) -> usize {
        self.segments.len()
    }

    /// Return true if the path is empty.
    pub fn is_empty(&self) -> bool {
        self.segments.len() == 0
    }

    /// Iterate over the stored segments.
    pub fn segments(&self) -> slice::Iter<'_, BSpline<P, K, C, D>> {
        self.segments.iter()
    }

    /// Push a segment, returning false if the path is at capacity.
    pub fn push(&mut self, segment: BSpline<P, K, C, D>) -> bool {
        if self.segments.len() < self.segments.capacity() {
            self.segments.push(segment);
            true
        } else {
            false
        }
    }

    /// Evaluate a point along the path for `t` in `[0, 1]`.
    /// Values outside the range are clamped. Returns None for empty paths.
    pub fn eval(&self, t: P::Scalar) -> Option<P>
    where
        [(); D + 1]: Sized,
    {
        let (index, local_t) = self.segment_parameter(t)?;
        let segment = &self.segments[index];
        let (kmin, kmax) = segment.knot_domain();
        let t_mapped = kmin + (kmax - kmin) * local_t;
        segment.eval(t_mapped).ok()
    }

    /// Return the start point of the path.
    pub fn start(&self) -> Option<P>
    where
        [(); D + 1]: Sized,
    {
        let zero = <P::Scalar as NumCast>::from(0.0).unwrap();
        self.eval(zero)
    }

    /// Return the end point of the path.
    pub fn end(&self) -> Option<P>
    where
        [(); D + 1]: Sized,
    {
        let one = <P::Scalar as NumCast>::from(1.0).unwrap();
        self.eval(one)
    }

    /// Return a path with reversed direction.
    pub fn reverse(&self) -> Self {
        let mut reversed = BSplinePath::new();
        for segment in self.segments.iter().rev() {
            let _ = reversed.push(segment.reverse());
        }
        reversed
    }

    /// Return the unit tangent direction at `t`.
    pub fn tangent(&self, t: P::Scalar) -> Option<P>
    where
        P: PointNorm,
        [(); D - 1]: Sized,
        [(); C - 1]: Sized,
        [(); K - 2]: Sized,
        [(); D]: Sized,
        [(); D + 1]: Sized,
        [(); (D - 1) + 1]: Sized, // keep (D - 1) + 1; compiler doesn't normalize to D
    {
        let (index, local_t) = self.segment_parameter(t)?;
        let segment = &self.segments[index];
        let (kmin, kmax) = segment.knot_domain();
        let t_mapped = kmin + (kmax - kmin) * local_t;
        segment.tangent(t_mapped).ok()
    }

    /// Return the curvature magnitude at `t`.
    pub fn curvature(&self, t: P::Scalar) -> Option<P::Scalar>
    where
        P: PointNorm + PointDot,
        [(); D + 1]: Sized,
        [(); D - 2]: Sized,
    {
        let (index, local_t) = self.segment_parameter(t)?;
        let segment = &self.segments[index];
        let (kmin, kmax) = segment.knot_domain();
        let t_mapped = kmin + (kmax - kmin) * local_t;
        segment.curvature(t_mapped).ok()
    }

    /// Return the principal normal direction at `t`.
    ///
    /// Returns `None` if the path is empty or curvature is undefined.
    pub fn normal(&self, t: P::Scalar) -> Option<P>
    where
        P: PointNorm + PointDot,
        [(); D + 1]: Sized,
        [(); D - 2]: Sized,
    {
        let (index, local_t) = self.segment_parameter(t)?;
        let segment = &self.segments[index];
        let (kmin, kmax) = segment.knot_domain();
        let t_mapped = kmin + (kmax - kmin) * local_t;
        segment.normal(t_mapped).ok().flatten()
    }

    /// Approximate parameter `t` at arc length `s`.
    pub fn t_at_length_approx(&self, s: P::Scalar, nsteps: usize) -> Option<P::Scalar>
    where
        P: PointNorm,
        [(); D + 1]: Sized,
    {
        let zero = <P::Scalar as NumCast>::from(0.0).unwrap();
        let one = <P::Scalar as NumCast>::from(1.0).unwrap();
        let half = <P::Scalar as NumCast>::from(0.5).unwrap();

        let total = self.arclen_partial(one, nsteps)?;
        if total <= P::Scalar::epsilon() {
            return Some(zero);
        }
        let target = s.clamp(zero, total);

        let mut lo = zero;
        let mut hi = one;
        for _ in 0..MAX_LENGTH_ITERS {
            let mid = (lo + hi) * half;
            let len = self.arclen_partial(mid, nsteps)?;
            if len < target {
                lo = mid;
            } else {
                hi = mid;
            }
        }

        Some((lo + hi) * half)
    }

    /// Approximate parameter `t` at arc length `s` using a default resolution.
    pub fn t_at_length(&self, s: P::Scalar) -> Option<P::Scalar>
    where
        P: PointNorm,
        [(); D + 1]: Sized,
    {
        self.t_at_length_approx(s, DEFAULT_LENGTH_STEPS)
    }

    /// Evaluate the point at arc length `s`.
    pub fn point_at_length_approx(&self, s: P::Scalar, nsteps: usize) -> Option<P>
    where
        P: PointNorm,
        [(); D + 1]: Sized,
    {
        let t = self.t_at_length_approx(s, nsteps)?;
        self.eval(t)
    }

    /// Evaluate the point at arc length `s` using a default resolution.
    pub fn point_at_length(&self, s: P::Scalar) -> Option<P>
    where
        P: PointNorm,
        [(); D + 1]: Sized,
    {
        self.point_at_length_approx(s, DEFAULT_LENGTH_STEPS)
    }

    /// Return a conservative bounding box based on all control points.
    /// Returns None for empty paths.
    pub fn bounding_box(&self) -> Option<[(P::Scalar, P::Scalar); P::DIM]>
    where
        P: PointIndex,
    {
        let mut iter = self.segments.iter();
        let first = match iter.next() {
            Some(segment) => segment_control_bounds(segment),
            None => return None,
        };

        let mut bounds = first;
        for segment in iter {
            let segment_bounds = segment_control_bounds(segment);
            for dim in 0..P::DIM {
                if segment_bounds[dim].0 < bounds[dim].0 {
                    bounds[dim].0 = segment_bounds[dim].0;
                }
                if segment_bounds[dim].1 > bounds[dim].1 {
                    bounds[dim].1 = segment_bounds[dim].1;
                }
            }
        }

        Some(bounds)
    }

    fn segment_parameter(&self, t: P::Scalar) -> Option<(usize, P::Scalar)> {
        let count = self.segments.len();
        if count == 0 {
            return None;
        }

        let zero = <P::Scalar as NumCast>::from(0.0).unwrap();
        let one = <P::Scalar as NumCast>::from(1.0).unwrap();
        let t = t.clamp(zero, one);

        let count_scalar = <P::Scalar as NumCast>::from(count as f64).unwrap();
        let scaled = t * count_scalar;
        if scaled >= count_scalar {
            return Some((count - 1, one));
        }

        let index = num_traits::cast::<P::Scalar, usize>(scaled.floor()).unwrap_or(0);
        let index_scalar = <P::Scalar as NumCast>::from(index as f64).unwrap();
        let local = scaled - index_scalar;
        Some((index, local))
    }

    fn arclen_partial(&self, t: P::Scalar, nsteps: usize) -> Option<P::Scalar>
    where
        P: PointNorm,
        [(); D + 1]: Sized,
    {
        let zero = <P::Scalar as NumCast>::from(0.0).unwrap();
        let one = <P::Scalar as NumCast>::from(1.0).unwrap();
        let t = t.clamp(zero, one);
        if t <= zero {
            return Some(zero);
        }

        let nsteps = nsteps.max(1);
        let nsteps_scalar = <P::Scalar as NumCast>::from(nsteps as f64).unwrap();
        let dt = t / nsteps_scalar;
        let mut arclen = zero;
        let mut prev = self.eval(zero)?;

        for i in 1..=nsteps {
            let i_scalar = <P::Scalar as NumCast>::from(i as f64).unwrap();
            let ti = if i == nsteps { t } else { dt * i_scalar };
            let p = self.eval(ti)?;
            arclen = arclen + (p - prev).squared_norm().sqrt();
            prev = p;
        }

        Some(arclen)
    }
}

impl<P, const K: usize, const C: usize, const D: usize, const N: usize> Default
    for BSplinePath<P, K, C, D, N>
where
    P: Point,
    [BSpline<P, K, C, D>; N]: tinyvec::Array<Item = BSpline<P, K, C, D>>,
{
    fn default() -> Self {
        Self::new()
    }
}

fn segment_control_bounds<P, const K: usize, const C: usize, const D: usize>(
    segment: &BSpline<P, K, C, D>,
) -> [(P::Scalar, P::Scalar); P::DIM]
where
    P: PointIndex,
{
    let mut iter = segment.control_points();
    let first = *iter
        .next()
        .expect("BSpline must have at least one control point");
    let mut bounds = [(
        <P::Scalar as NumCast>::from(0.0).unwrap(),
        <P::Scalar as NumCast>::from(0.0).unwrap(),
    ); P::DIM];
    for dim in 0..P::DIM {
        let value = first[dim];
        bounds[dim] = (value, value);
    }

    for point in iter {
        for dim in 0..P::DIM {
            let value = point[dim];
            if value < bounds[dim].0 {
                bounds[dim].0 = value;
            }
            if value > bounds[dim].1 {
                bounds[dim].1 = value;
            }
        }
    }

    bounds
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{EPSILON, PointN, PointNorm};

    #[test]
    fn bspline_path_eval_segments() {
        let knots: [f64; 4] = [0.0, 0.0, 1.0, 1.0];
        let s1: BSpline<PointN<f64, 2>, 4, 2, 1> =
            BSpline::new(knots, [PointN::new([0.0, 0.0]), PointN::new([1.0, 0.0])]).unwrap();
        let s2: BSpline<PointN<f64, 2>, 4, 2, 1> =
            BSpline::new(knots, [PointN::new([1.0, 0.0]), PointN::new([1.0, 1.0])]).unwrap();

        let mut path: BSplinePath<PointN<f64, 2>, 4, 2, 1, 4> = BSplinePath::new();
        path.push(s1);
        path.push(s2);

        let p1 = path.eval(0.25).unwrap();
        assert!((p1 - PointN::new([0.5, 0.0])).squared_norm() < EPSILON);

        let p2 = path.eval(0.75).unwrap();
        assert!((p2 - PointN::new([1.0, 0.5])).squared_norm() < EPSILON);
    }

    #[test]
    fn bspline_path_clamps_out_of_range() {
        let knots: [f64; 4] = [0.0, 0.0, 1.0, 1.0];
        let s1: BSpline<PointN<f64, 2>, 4, 2, 1> =
            BSpline::new(knots, [PointN::new([0.0, 0.0]), PointN::new([1.0, 0.0])]).unwrap();
        let s2: BSpline<PointN<f64, 2>, 4, 2, 1> =
            BSpline::new(knots, [PointN::new([1.0, 0.0]), PointN::new([1.0, 1.0])]).unwrap();

        let mut path: BSplinePath<PointN<f64, 2>, 4, 2, 1, 4> = BSplinePath::new();
        path.push(s1);
        path.push(s2);

        let start = path.eval(-1.0).unwrap();
        assert_eq!(start, PointN::new([0.0, 0.0]));

        let end = path.eval(2.0).unwrap();
        assert_eq!(end, PointN::new([1.0, 1.0]));
    }

    #[test]
    fn bspline_path_bounds_control_points() {
        let knots: [f64; 4] = [0.0, 0.0, 1.0, 1.0];
        let s1: BSpline<PointN<f64, 2>, 4, 2, 1> =
            BSpline::new(knots, [PointN::new([0.0, 2.0]), PointN::new([1.0, -1.0])]).unwrap();
        let s2: BSpline<PointN<f64, 2>, 4, 2, 1> =
            BSpline::new(knots, [PointN::new([-1.0, 0.5]), PointN::new([2.0, -2.0])]).unwrap();

        let mut path: BSplinePath<PointN<f64, 2>, 4, 2, 1, 4> = BSplinePath::new();
        path.push(s1);
        path.push(s2);

        let bounds = path.bounding_box().unwrap();
        assert_eq!(bounds[0], (-1.0, 2.0));
        assert_eq!(bounds[1], (-2.0, 2.0));
    }

    #[test]
    fn bspline_path_api_parity() {
        let knots: [f64; 4] = [0.0, 0.0, 1.0, 1.0];
        let s1: BSpline<PointN<f64, 2>, 4, 2, 1> =
            BSpline::new(knots, [PointN::new([0.0, 0.0]), PointN::new([1.0, 0.0])]).unwrap();
        let s2: BSpline<PointN<f64, 2>, 4, 2, 1> =
            BSpline::new(knots, [PointN::new([1.0, 0.0]), PointN::new([2.0, 0.0])]).unwrap();

        let mut path: BSplinePath<PointN<f64, 2>, 4, 2, 1, 2> = BSplinePath::new();
        path.push(s1);
        path.push(s2);

        assert_eq!(path.start().unwrap(), PointN::new([0.0, 0.0]));
        assert_eq!(path.end().unwrap(), PointN::new([2.0, 0.0]));

        let reversed = path.reverse();
        let p0 = path.eval(0.25).unwrap();
        let p1 = reversed.eval(0.75).unwrap();
        assert!((p0 - p1).squared_norm() < EPSILON);

        let tangent = path.tangent(0.25).unwrap();
        assert!((tangent[0] - 1.0).abs() < EPSILON);
        assert!(tangent[1].abs() < EPSILON);

        let t = path.t_at_length(1.0).unwrap();
        assert!((t - 0.5).abs() < 1e-6);

        let p = path.point_at_length(1.0).unwrap();
        assert!((p - PointN::new([1.0, 0.0])).squared_norm() < EPSILON);
    }

    #[test]
    fn bspline_path_curvature_nonzero() {
        let knots: [f64; 6] = [0.0, 0.0, 0.0, 1.0, 1.0, 1.0];
        let seg: BSpline<PointN<f64, 2>, 6, 3, 2> = BSpline::new(
            knots,
            [
                PointN::new([0.0, 0.0]),
                PointN::new([1.0, 1.0]),
                PointN::new([2.0, 0.0]),
            ],
        )
        .unwrap();

        let mut path: BSplinePath<PointN<f64, 2>, 6, 3, 2, 1> = BSplinePath::new();
        path.push(seg);

        let curvature = path.curvature(0.5).unwrap();
        assert!(curvature > 0.0);

        let normal = path.normal(0.5).unwrap();
        assert!((normal.squared_norm() - 1.0).abs() < 1e-6);
    }
}
