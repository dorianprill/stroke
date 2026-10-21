use num_traits::{Float, NumCast, One, Zero};

/// A scalar with a separate binary exponent for intermediate curve arithmetic.
/// The mantissa retains the scalar's precision; the exponent extends its range.
#[derive(Clone, Copy)]
pub(crate) struct Scaled<T> {
    fraction: T,
    exponent: i32,
}

impl<T: Float> Scaled<T> {
    pub(crate) fn new(value: T) -> Self {
        Self::normalize(value, 0)
    }

    fn normalize(value: T, exponent: i32) -> Self {
        if value == T::zero() || !value.is_finite() {
            return Self {
                fraction: value,
                exponent: 0,
            };
        }
        let (mantissa, power, sign) = value.integer_decode();
        let bits = (u64::BITS - 1 - mantissa.leading_zeros()) as i32;
        let two = T::one() + T::one();
        let fraction = <T as NumCast>::from(mantissa).unwrap() * two.powi(-bits);
        Self {
            fraction: if sign < 0 { -fraction } else { fraction },
            exponent: exponent + power as i32 + bits,
        }
    }

    pub(crate) fn sum(self, other: Self) -> Self {
        if self.fraction == T::zero() {
            return other;
        }
        if other.fraction == T::zero() {
            return self;
        }
        let (large, small) = if self.exponent >= other.exponent {
            (self, other)
        } else {
            (other, self)
        };
        let two = T::one() + T::one();
        Self::normalize(
            large.fraction + small.fraction * two.powi(small.exponent - large.exponent),
            large.exponent,
        )
    }

    pub(crate) fn difference(self, other: Self) -> Self {
        self.sum(Self {
            fraction: -other.fraction,
            exponent: other.exponent,
        })
    }

    pub(crate) fn product(self, other: Self) -> Self {
        Self::normalize(
            self.fraction * other.fraction,
            self.exponent + other.exponent,
        )
    }

    pub(crate) fn quotient(self, other: Self) -> Self {
        Self::normalize(
            self.fraction / other.fraction,
            self.exponent - other.exponent,
        )
    }

    pub(crate) fn value(self) -> T {
        let two = T::one() + T::one();
        let first = self.exponent / 2;
        (self.fraction * two.powi(first)) * two.powi(self.exponent - first)
    }
}

/// Restore a short product/quotient without overflowing or underflowing its
/// scale factors first. `Float::integer_decode` specifies a binary decomposition;
/// this uses that trait contract, not an IEEE layout or a concrete scalar type.
pub(crate) fn scaled_product_ratio<T: Float>(
    numerators: &[T],
    denominators: &[T],
    power: i32,
) -> T {
    let one = T::one();
    let two = one + one;
    let half = one / two;
    let mut fraction = one;
    let mut exponent = power;
    let mut sign = one;
    for (factors, divide) in [(numerators, false), (denominators, true)] {
        for &value in factors {
            if value == T::zero() {
                return T::zero();
            }
            let (mantissa, exp, value_sign) = value.integer_decode();
            let bits = (u64::BITS - 1 - mantissa.leading_zeros()) as i32;
            let normalized = <T as NumCast>::from(mantissa).unwrap() * two.powi(-bits);
            let exp = exp as i32 + bits;
            if value_sign < 0 {
                sign = -sign;
            }
            if divide {
                fraction = fraction / normalized;
                exponent -= exp;
                if fraction < one {
                    fraction = fraction * two;
                    exponent -= 1;
                }
            } else {
                fraction = fraction * normalized;
                exponent += exp;
                if fraction >= two {
                    fraction = fraction * half;
                    exponent += 1;
                }
            }
        }
    }
    // Splitting the exponent also handles results near the subnormal and
    // overflow limits without prematurely rounding the power of two to 0/inf.
    let first = exponent / 2;
    sign * ((fraction * two.powi(first)) * two.powi(exponent - first))
}

/// Midpoint as a high/low pair. Retaining the addition error avoids shifting a
/// line by an ulp when its midpoint itself falls between representable values.
pub(crate) fn midpoint_parts<T: Float>(a: T, b: T) -> (T, T) {
    let half = T::one() / (T::one() + T::one());
    let sum = a + b;
    let (a, b, scale) = if sum.is_finite() {
        (a, b, half)
    } else {
        (a * half, b * half, T::one())
    };
    let sum = a + b;
    let tail = sum - a;
    let error = (a - (sum - tail)) + (b - tail);
    (sum * scale, error * scale)
}

/// Divide a point without first overflowing the reciprocal of a tiny divisor.
pub(crate) fn divide<P: crate::Point>(value: P, divisor: P::Scalar) -> P {
    let reciprocal = P::Scalar::one() / divisor;
    if divisor != P::Scalar::zero() && !reciprocal.is_finite() {
        let factor = P::Scalar::one() / divisor.abs().sqrt();
        (value * factor) * (factor * divisor.signum())
    } else {
        value * reciprocal
    }
}

/// Difference quotient with finite endpoints, avoiding overflowing differences.
pub(crate) fn difference_quotient<P: crate::Point>(
    a: P,
    b: P,
    start: P::Scalar,
    end: P::Scalar,
) -> P {
    let width = end - start;
    if !width.is_finite() {
        let half = <P::Scalar as NumCast>::from(0.5).unwrap();
        return divide(b * half - a * half, end * half - start * half);
    }
    let difference = b - a;
    if difference.is_finite() {
        divide(difference, width)
    } else {
        divide(b, width) - divide(a, width)
    }
}

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
