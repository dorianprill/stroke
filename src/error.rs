//! Errors shared by curve constructors.

use core::fmt;

use crate::Point;

/// Errors returned when constructing a Bezier curve or line segment.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CurveConstructionError {
    /// A generic Bezier curve must contain at least one control point.
    NoControlPoints,
    /// A control point contains NaN or an infinite coordinate.
    NonFiniteControlPoint {
        /// Index of the first non-finite control point.
        index: usize,
    },
}

impl fmt::Display for CurveConstructionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NoControlPoints => write!(f, "a curve requires at least one control point"),
            Self::NonFiniteControlPoint { index } => {
                write!(f, "control point at index {index} must be finite")
            }
        }
    }
}

impl core::error::Error for CurveConstructionError {}

/// Errors returned when extracting a Bezier subcurve.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CurveParameterError {
    /// Endpoints must be finite and satisfy `0 <= start <= end <= 1`.
    InvalidRange,
    /// Curve arithmetic produced a non-finite control point.
    NonFiniteComputation,
}

impl fmt::Display for CurveParameterError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidRange => f.write_str("invalid subcurve parameter range"),
            Self::NonFiniteComputation => {
                f.write_str("subcurve arithmetic produced a non-finite point")
            }
        }
    }
}

impl core::error::Error for CurveParameterError {}

pub(crate) fn validate_control_points<P: Point>(
    points: &[P],
) -> Result<(), CurveConstructionError> {
    match points.iter().position(|point| !point.is_finite()) {
        Some(index) => Err(CurveConstructionError::NonFiniteControlPoint { index }),
        None => Ok(()),
    }
}
