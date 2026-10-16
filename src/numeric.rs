use num_traits::{Float, NumCast, One, Zero};

/// Affine interpolation that preserves exact endpoints and avoids subtracting
/// opposite finite extremes unless that subtraction is representable.
pub(crate) fn lerp<P: crate::Point>(start: P, end: P, t: P::Scalar) -> P {
    if t == P::Scalar::zero() {
        return start;
    }
    if t == P::Scalar::one() {
        return end;
    }
    let difference = end - start;
    if difference.is_finite() {
        start + difference * t
    } else {
        start * (P::Scalar::one() - t) + end * t
    }
}

/// Map a parameter inside a finite, positive-width interval to [0, 1],
/// including intervals whose width overflows the scalar representation.
pub(crate) fn unit_parameter<T: Float>(t: T, start: T, end: T) -> T {
    let width = end - start;
    if width.is_finite() {
        (t - start) / width
    } else {
        let half = <T as NumCast>::from(0.5).unwrap();
        (t * half - start * half) / (end * half - start * half)
    }
}

/// Visit Bernstein polynomial roots in parameter order. Zero coefficients are
/// ignored when counting sign variations; an identically-zero polynomial has
/// no isolated roots. Adjacent subdivisions may report the same boundary root.
pub(crate) fn bernstein_roots<T: Float, const N: usize>(
    control: [T; N],
    tolerance: T,
    mut emit: impl FnMut(T),
) {
    let scale = control
        .iter()
        .fold(T::zero(), |scale, value| scale.max(value.abs()));
    if scale == T::zero() || !scale.is_finite() {
        return;
    }
    let control = control.map(|value| value / scale);
    roots_recursive(control, T::zero(), T::one(), tolerance, 0, &mut emit);
}

fn roots_recursive<T: Float, const N: usize>(
    control: [T; N],
    start: T,
    end: T,
    tolerance: T,
    depth: usize,
    emit: &mut impl FnMut(T),
) {
    if control.iter().all(|value| *value == T::zero()) {
        return;
    }
    if control[0] == T::zero() {
        emit(start);
    }
    let mut previous = T::zero();
    let mut variations = 0;
    for value in control {
        if value == T::zero() {
            continue;
        }
        if previous != T::zero() && value.signum() != previous {
            variations += 1;
        }
        previous = value.signum();
    }
    let half = <T as NumCast>::from(0.5).unwrap();
    let mid = start * half + end * half;
    if variations > 0 {
        if depth >= 64 || end - start <= tolerance || mid == start || mid == end {
            let a = control[0];
            let b = control[N - 1];
            let root = if a != T::zero() && b != T::zero() && a.signum() != b.signum() {
                start + (end - start) * (-a / (b - a))
            } else {
                mid
            };
            emit(root);
        } else {
            let mut left = control;
            let mut right = control;
            let mut work = control;
            for level in 0..N {
                left[level] = work[0];
                right[N - 1 - level] = work[N - 1 - level];
                for i in 0..N - 1 - level {
                    work[i] = work[i] * half + work[i + 1] * half;
                }
            }
            roots_recursive(left, start, mid, tolerance, depth + 1, emit);
            roots_recursive(right, mid, end, tolerance, depth + 1, emit);
        }
    }
    if control[N - 1] == T::zero() {
        emit(end);
    }
}

/// Return whether `value` is negligible relative to quantities with magnitude
/// `scale`. This is intended for algebraic cancellation tests where both values
/// have the same units, not for geometric distance comparisons.
pub(crate) fn is_negligible<T: Float>(value: T, scale: T) -> bool {
    if value == T::zero() {
        return true;
    }
    if !value.is_finite() || !scale.is_finite() {
        return false;
    }

    let roundoff_factor = <T as NumCast>::from(8_u8).unwrap();
    let scale = scale.abs();
    scale != T::zero() && value.abs() / scale <= T::epsilon() * roundoff_factor
}
