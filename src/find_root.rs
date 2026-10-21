//! Root-finding helpers for curve types that support per-axis evaluation.

use crate::point::PointIndex;
use crate::roots::{RootFindingError, root_newton_raphson};
use num_traits::{Float, NumCast};

/// Find parameters where a curve coordinate equals a target value.
pub trait FindRoot<P: PointIndex> {
    /// Return the inclusive parameter domain for the curve.
    fn parameter_domain(&self) -> (P::Scalar, P::Scalar);

    /// Evaluate the curve along one axis.
    fn axis_value(&self, t: P::Scalar, axis: usize) -> Result<P::Scalar, RootFindingError>;

    /// Evaluate the curve's derivative along one axis.
    fn axis_derivative(&self, t: P::Scalar, axis: usize) -> Result<P::Scalar, RootFindingError>;

    /// Find a parameter where coordinate `axis` equals `value`, using Newton's method.
    /// `eps` bounds the coordinate error and defaults to `1e-6`. It must be finite
    /// and positive. Clamp the initial guess and each update to the curve's domain.
    /// Check the coordinate error at the returned parameter before reporting success.
    /// This local search can fail even when a solution exists.
    fn root_newton_axis(
        &self,
        value: P::Scalar,
        axis: usize,
        start: P::Scalar,
        eps: Option<P::Scalar>,
        max_iter: Option<usize>,
    ) -> Result<P::Scalar, RootFindingError> {
        if axis >= P::DIM {
            return Err(RootFindingError::InvalidAxis);
        }
        if !value.is_finite() {
            return Err(RootFindingError::InvalidParameter);
        }
        let eps = eps.unwrap_or_else(|| <P::Scalar as NumCast>::from(1e-6).unwrap());
        let max_iter = max_iter.unwrap_or(64);
        root_newton_raphson(
            start,
            self.parameter_domain(),
            |t| self.axis_value(t, axis).map(|v| v - value),
            |t| self.axis_derivative(t, axis),
            eps,
            max_iter,
        )
    }
}
