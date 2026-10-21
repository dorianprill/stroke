//! Euclidean distances and nearest-point estimates from sampling and local search.

use crate::PointNorm;
use num_traits::{Float, NumCast, Zero};

const LOCAL_SEARCH_ITERS: usize = 16;

/// Coincidence relative to the first control point, without squaring distances.
pub(crate) fn is_a_point<P: PointNorm>(points: &[P], tolerance: P::Scalar) -> bool {
    if !tolerance.is_finite() {
        return false;
    }
    let tolerance = tolerance.abs();
    let first = points[0];
    points.iter().all(|&point| {
        if !first.is_finite() || !point.is_finite() {
            return false;
        }
        if tolerance == P::Scalar::zero() {
            (0..P::DIM).all(|axis| first[axis] == point[axis])
        } else {
            between(first, point) <= tolerance
        }
    })
}

/// Compute Euclidean distance with repeated `hypot` calls.
/// Return infinity for a distance outside the scalar range, or NaN if either
/// point has a non-finite coordinate.
pub(crate) fn between<P: PointNorm>(a: P, b: P) -> P::Scalar {
    if !a.is_finite() || !b.is_finite() {
        return P::Scalar::nan();
    }
    let mut distance = P::Scalar::zero();
    for axis in 0..P::DIM {
        distance = distance.hypot(a[axis] - b[axis]);
    }
    distance
}

/// Sample the unit interval, then use ternary search near the closest sample.
/// Retain the best distance found. The caller maps the unit parameter to its
/// curve's domain. This estimate has no global-minimum or error guarantee.
pub(crate) fn approximate<P: PointNorm, E>(
    point: P,
    segments: usize,
    mut eval: impl FnMut(P::Scalar) -> Result<P, E>,
) -> Result<P::Scalar, E> {
    if !point.is_finite() {
        return Ok(P::Scalar::nan());
    }
    let segments = segments.max(1);
    let count = <P::Scalar as NumCast>::from(segments).unwrap();
    let zero = P::Scalar::zero();
    let three = <P::Scalar as NumCast>::from(3.0).unwrap();
    let half = <P::Scalar as NumCast>::from(0.5).unwrap();
    let mut best_i = 0;
    let mut best = between(eval(zero)?, point);
    if best.is_nan() {
        return Ok(best);
    }
    for i in 1..=segments {
        let t = <P::Scalar as NumCast>::from(i).unwrap() / count;
        let distance = between(eval(t)?, point);
        if distance.is_nan() {
            return Ok(distance);
        }
        if distance < best {
            best = distance;
            best_i = i;
        }
    }
    let mut left = <P::Scalar as NumCast>::from(best_i.saturating_sub(1)).unwrap() / count;
    let right_i = if best_i == segments {
        segments
    } else {
        best_i + 1
    };
    let mut right = <P::Scalar as NumCast>::from(right_i).unwrap() / count;
    for _ in 0..LOCAL_SEARCH_ITERS {
        let third = (right - left) / three;
        let t1 = left + third;
        let t2 = right - third;
        let d1 = between(eval(t1)?, point);
        let d2 = between(eval(t2)?, point);
        if d1.is_nan() || d2.is_nan() {
            return Ok(P::Scalar::nan());
        }
        best = best.min(d1).min(d2);
        if d1 < d2 {
            right = t2;
        } else {
            left = t1;
        }
    }
    let distance = between(eval(left * half + right * half)?, point);
    if distance.is_nan() {
        Ok(distance)
    } else {
        Ok(best.min(distance))
    }
}
