//! Root-finding helpers for curve types that support per-axis evaluation.

use crate::point::PointIndex;
use crate::roots::{RootFindingError, root_newton_raphson};
use num_traits::{Float, NumCast};

/// Helper trait for 1D root finding along a curve axis.
pub trait FindRoot<P: PointIndex> {
    /// Return the inclusive parameter domain for the curve.
    fn parameter_domain(&self) -> (P::Scalar, P::Scalar);

    /// Evaluate the curve along one axis.
    fn axis_value(&self, t: P::Scalar, axis: usize) -> Result<P::Scalar, RootFindingError>;

    /// Evaluate the curve's derivative along one axis.
    fn axis_derivative(&self, t: P::Scalar, axis: usize) -> Result<P::Scalar, RootFindingError>;

    /// Find a root for a particular axis using Newton-Raphson on the scalar axis function.
    /// `eps` is a positive finite absolute tolerance in coordinate units (default
    /// `1e-6`), not parameter units. Finite initial guesses and every Newton step
    /// are clamped to the domain. Success checks the actual returned residual.
    /// This is a local search: failure does not prove that no root exists.
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
