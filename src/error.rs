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

pub(crate) fn validate_control_points<P: Point>(
    points: &[P],
) -> Result<(), CurveConstructionError> {
    match points.iter().position(|point| !point.is_finite()) {
        Some(index) => Err(CurveConstructionError::NonFiniteControlPoint { index }),
        None => Ok(()),
    }
}
