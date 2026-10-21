//! Adaptive Gauss-Legendre integration and inverse arc-length lookup.

use core::fmt;

use num_traits::{Float, NumCast};

const MAX_DEPTH: usize = 48;
const MAX_PANELS: usize = 4096;

/// Errors produced while numerically evaluating an arc length.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ArcLengthError {
    /// The requested accuracy must be finite and greater than zero.
    InvalidAccuracy,
    /// A parameter or length is non-finite, or a parameter range is invalid.
    InvalidParameter,
    /// Curve arithmetic produced a non-finite speed or integral.
    NonFiniteComputation,
    /// Subdivision reached its work limit before meeting the error estimate.
    FailedToConverge,
    /// A curve evaluation failed within the requested parameter range.
    CurveEvaluationFailed,
}

impl fmt::Display for ArcLengthError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ArcLengthError::InvalidAccuracy => {
                write!(
                    f,
                    "arc-length accuracy must be finite and greater than zero"
                )
            }
            ArcLengthError::InvalidParameter => write!(f, "invalid arc-length parameter"),
            ArcLengthError::NonFiniteComputation => {
                write!(f, "arc-length evaluation produced a non-finite value")
            }
            ArcLengthError::FailedToConverge => {
                write!(f, "arc-length integration failed to converge")
            }
            ArcLengthError::CurveEvaluationFailed => {
                write!(f, "curve evaluation failed inside its valid domain")
            }
        }
    }
}

impl core::error::Error for ArcLengthError {}

pub(crate) fn inverse<T: Float>(
    s: T,
    accuracy: T,
    total: T,
    mut length: impl FnMut(T, T, T) -> Result<T, ArcLengthError>,
) -> Result<T, ArcLengthError> {
    validate_accuracy(accuracy)?;
    if !s.is_finite() {
        return Err(ArcLengthError::InvalidParameter);
    }
    let zero = T::zero();
    let one = T::one();
    if total == zero || s <= zero {
        return Ok(zero);
    }
    if s >= total {
        return Ok(one);
    }
    let half = <T as NumCast>::from(0.5).unwrap();
    let budget = accuracy / <T as NumCast>::from(128_u8).unwrap();
    let budget = budget.max(T::min_positive_value() * T::epsilon());
    let tolerance =
        (accuracy * half).max(total * T::epsilon() * <T as NumCast>::from(8_u8).unwrap());
    let mut remaining = s;
    let mut lo = zero;
    let mut hi = one;
    for _ in 0..64 {
        let mid = lo * half + hi * half;
        if mid == lo || mid == hi {
            return Err(ArcLengthError::FailedToConverge);
        }
        let left = length(lo, mid, budget)?;
        if (left - remaining).abs() <= tolerance {
            return Ok(mid);
        }
        if left < remaining {
            lo = mid;
            remaining = remaining - left;
        } else {
            hi = mid;
        }
    }
    Err(ArcLengthError::FailedToConverge)
}

pub(crate) fn validate_accuracy<T: Float>(accuracy: T) -> Result<(), ArcLengthError> {
    if accuracy.is_finite() && accuracy > T::zero() {
        Ok(())
    } else {
        Err(ArcLengthError::InvalidAccuracy)
    }
}

/// Adaptively integrate `f` over `[a, b]` using paired 8- and 16-point
/// Gauss-Legendre estimates.
pub(crate) fn integrate<T, F>(a: T, b: T, accuracy: T, mut f: F) -> Result<T, ArcLengthError>
where
    T: Float,
    F: FnMut(T) -> Result<T, ArcLengthError>,
{
    validate_accuracy(accuracy)?;
    if a == b {
        return Ok(T::zero());
    }

    let mut panels = 0;
    integrate_recursive(a, b, accuracy, 0, &mut panels, &mut f)
}

fn integrate_recursive<T, F>(
    a: T,
    b: T,
    accuracy: T,
    depth: usize,
    panels: &mut usize,
    f: &mut F,
) -> Result<T, ArcLengthError>
where
    T: Float,
    F: FnMut(T) -> Result<T, ArcLengthError>,
{
    *panels += 1;
    if *panels > MAX_PANELS {
        return Err(ArcLengthError::FailedToConverge);
    }

    let coarse = gauss_legendre_8(a, b, f)?;
    let refined = gauss_legendre_16(a, b, f)?;
    if !coarse.is_finite() || !refined.is_finite() {
        return Err(ArcLengthError::NonFiniteComputation);
    }

    let roundoff_factor = <T as NumCast>::from(64_u8).unwrap();
    let roundoff_floor = refined.abs() * T::epsilon() * roundoff_factor;
    if (refined - coarse).abs() <= accuracy.max(roundoff_floor) {
        return Ok(refined);
    }
    if depth >= MAX_DEPTH {
        return Err(ArcLengthError::FailedToConverge);
    }

    let half = <T as NumCast>::from(0.5).unwrap();
    let mid = a * half + b * half;
    if mid == a || mid == b {
        return Err(ArcLengthError::FailedToConverge);
    }
    let child_accuracy = accuracy * half;
    let left = integrate_recursive(a, mid, child_accuracy, depth + 1, panels, f)?;
    let right = integrate_recursive(mid, b, child_accuracy, depth + 1, panels, f)?;
    let result = left + right;
    if result.is_finite() {
        Ok(result)
    } else {
        Err(ArcLengthError::NonFiniteComputation)
    }
}

fn gauss_legendre_8<T, F>(a: T, b: T, f: &mut F) -> Result<T, ArcLengthError>
where
    T: Float,
    F: FnMut(T) -> Result<T, ArcLengthError>,
{
    const NODES: [f64; 4] = [
        0.183_434_642_495_649_8,
        0.525_532_409_916_329,
        0.796_666_477_413_626_7,
        0.960_289_856_497_536_3,
    ];
    const WEIGHTS: [f64; 4] = [
        0.362_683_783_378_362,
        0.313_706_645_877_887_3,
        0.222_381_034_453_374_5,
        0.101_228_536_290_376_3,
    ];
    gauss_legendre_symmetric(a, b, &NODES, &WEIGHTS, f)
}

fn gauss_legendre_16<T, F>(a: T, b: T, f: &mut F) -> Result<T, ArcLengthError>
where
    T: Float,
    F: FnMut(T) -> Result<T, ArcLengthError>,
{
    const NODES: [f64; 8] = [
        0.095_012_509_837_637_44,
        0.281_603_550_779_258_9,
        0.458_016_777_657_227_4,
        0.617_876_244_402_643_8,
        0.755_404_408_355_003,
        0.865_631_202_387_831_8,
        0.944_575_023_073_232_6,
        0.989_400_934_991_649_9,
    ];
    const WEIGHTS: [f64; 8] = [
        0.189_450_610_455_068_5,
        0.182_603_415_044_923_6,
        0.169_156_519_395_002_54,
        0.149_595_988_816_576_74,
        0.124_628_971_255_533_88,
        0.095_158_511_682_492_78,
        0.062_253_523_938_647_89,
        0.027_152_459_411_754_096,
    ];
    gauss_legendre_symmetric(a, b, &NODES, &WEIGHTS, f)
}

fn gauss_legendre_symmetric<T, F, const N: usize>(
    a: T,
    b: T,
    nodes: &[f64; N],
    weights: &[f64; N],
    f: &mut F,
) -> Result<T, ArcLengthError>
where
    T: Float,
    F: FnMut(T) -> Result<T, ArcLengthError>,
{
    let half = <T as NumCast>::from(0.5).unwrap();
    let midpoint = a * half + b * half;
    let radius = (b - a) * half;
    let mut sum = T::zero();

    for (&node, &weight) in nodes.iter().zip(weights) {
        let node = <T as NumCast>::from(node).unwrap();
        let weight = <T as NumCast>::from(weight).unwrap();
        let offset = radius * node;
        let left = f(midpoint - offset)?;
        let right = f(midpoint + offset)?;
        if !left.is_finite() || !right.is_finite() {
            return Err(ArcLengthError::NonFiniteComputation);
        }
        sum = sum + weight * left * radius + weight * right * radius;
    }

    Ok(sum)
}
