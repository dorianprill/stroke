use num_traits::{Float, NumCast};

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
