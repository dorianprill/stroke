//! Errors shared by curve construction, parameter queries and Bezier bases.

use core::fmt;

use crate::Point;

/// Invalid control points supplied to a Bezier or line constructor.
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

/// Failure to evaluate Bernstein basis values or their parameter derivatives.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BezierBasisError {
    /// The parameter is NaN or infinite. Finite extrapolation is allowed.
    NonFiniteParameter,
    /// A requested basis value or derivative is not representable as a finite scalar.
    NonFiniteComputation,
}

impl fmt::Display for BezierBasisError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::NonFiniteParameter => "the Bezier basis parameter must be finite",
            Self::NonFiniteComputation => "Bezier basis arithmetic produced a non-finite value",
        })
    }
}

impl core::error::Error for BezierBasisError {}

pub(crate) fn validate_control_points<P: Point>(
    points: &[P],
) -> Result<(), CurveConstructionError> {
    match points.iter().position(|point| !point.is_finite()) {
        Some(index) => Err(CurveConstructionError::NonFiniteControlPoint { index }),
        None => Ok(()),
    }
}
