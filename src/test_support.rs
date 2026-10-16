use core::fmt::Debug;

use num_traits::Float;

use crate::PointNorm;

#[track_caller]
pub(crate) fn assert_approx_equal<T>(
    actual: T,
    expected: T,
    absolute_tolerance: T,
    relative_tolerance: T,
) where
    T: Float + Debug,
{
    let error = (actual - expected).abs();
    let scale = actual.abs().max(expected.abs());
    let allowed = absolute_tolerance.max(relative_tolerance * scale);

    assert!(
        error <= allowed,
        "actual={actual:?}, expected={expected:?}, error={error:?}, allowed={allowed:?}",
    );
}

#[track_caller]
pub(crate) fn assert_point_approx_equal<P>(
    actual: P,
    expected: P,
    absolute_tolerance: P::Scalar,
    relative_tolerance: P::Scalar,
) where
    P: PointNorm + Debug,
    P::Scalar: Debug,
{
    let error = (actual - expected).squared_norm().sqrt();
    let actual_magnitude = actual.squared_norm().sqrt();
    let expected_magnitude = expected.squared_norm().sqrt();
    let scale = actual_magnitude.max(expected_magnitude);
    let allowed = absolute_tolerance.max(relative_tolerance * scale);

    assert!(
        error <= allowed,
        "actual={actual:?}, expected={expected:?}, distance={error:?}, allowed={allowed:?}",
    );
}
