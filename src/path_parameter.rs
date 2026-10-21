//! Shared conversion between a normalized path parameter and segment-local coordinates.

use crate::ArcLengthError;
use crate::numeric::unit_parameter;
use num_traits::{Float, NumCast};

/// Compute the same segment boundaries for evaluation and inverse length lookup.
/// Multiplication cannot always undo a rounded division.
pub(crate) fn boundary<T: Float>(index: usize, count: usize) -> T {
    <T as NumCast>::from(index).unwrap() / <T as NumCast>::from(count).unwrap()
}

pub(crate) fn decode<T: Float>(t: T, count: usize) -> Option<(usize, T)> {
    if count == 0 || t.is_nan() {
        return None;
    }
    let t = t.clamp(T::zero(), T::one());
    if t == T::one() {
        return Some((count - 1, T::one()));
    }
    // Exact boundaries select the segment on the right. Nearby representable
    // parameters must stay in their own segment.
    let mut lo = 0;
    let mut hi = count;
    while hi - lo > 1 {
        let mid = lo + (hi - lo) / 2;
        if t < boundary(mid, count) {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    let start = boundary(lo, count);
    let end = boundary(lo + 1, count);
    Some((lo, unit_parameter(t, start, end)))
}

pub(crate) fn encode<T: Float>(index: usize, local: T, count: usize) -> Result<T, ArcLengthError> {
    let start: T = boundary(index, count);
    let end: T = boundary(index + 1, count);
    if start >= end {
        // The scalar cannot distinguish this path's segment boundaries.
        return Err(ArcLengthError::FailedToConverge);
    }
    let mut t = (start * (T::one() - local) + end * local).clamp(start, end);
    if index + 1 < count && t == end {
        // Rounding to the next boundary would select the wrong segment.
        // Bisect to find the last representable parameter inside this interval.
        // This works without assumptions about the scalar's bit layout.
        let half = <T as NumCast>::from(0.5).unwrap();
        let mut inside = start;
        loop {
            let mid = inside * half + end * half;
            if mid <= inside || mid >= end {
                break;
            }
            inside = mid;
        }
        t = inside;
    }
    Ok(t)
}

/// Intersect a path interval with one segment, retaining the same boundaries
/// used by evaluation and inverse-length lookup.
pub(crate) fn local_range<T: Float>(
    start: T,
    end: T,
    index: usize,
    count: usize,
) -> Option<(T, T)> {
    let lo = boundary(index, count);
    let hi = boundary(index + 1, count);
    let start = start.max(lo);
    let end = end.min(hi);
    (start < end).then(|| (unit_parameter(start, lo, hi), unit_parameter(end, lo, hi)))
}
