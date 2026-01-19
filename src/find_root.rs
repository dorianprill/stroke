use crate::point::Point;
use crate::roots::{root_newton_raphson, RootFindingError};
use crate::NativeFloat;

pub trait FindRoot<P: Point> {
    /// Return the inclusive parameter domain for the curve.
    fn parameter_domain(&self) -> (P::Scalar, P::Scalar);

    /// Evaluate the curve along one axis.
    fn axis_value(&self, t: P::Scalar, axis: usize) -> Result<P::Scalar, RootFindingError>;

    /// Evaluate the curve's derivative along one axis.
    fn axis_derivative(&self, t: P::Scalar, axis: usize) -> Result<P::Scalar, RootFindingError>;

    /// Find a root for a particular axis using Newton-Raphson on the scalar axis function.
    fn root_newton_axis(
        &self,
        value: P::Scalar,
        axis: usize,
        start: P::Scalar,
        eps: Option<P::Scalar>,
        max_iter: Option<usize>,
    ) -> Result<P::Scalar, RootFindingError> {
        let eps = eps.unwrap_or_else(|| P::Scalar::from(1e-6 as NativeFloat));
        let max_iter = max_iter.unwrap_or(64);
        let (kmin, kmax) = self.parameter_domain();

        let clamp_t = |mut t: P::Scalar| {
            if t < kmin {
                t = kmin;
            } else if t > kmax {
                t = kmax;
            }
            t
        };

        let start = clamp_t(start);
        let fx = |x: P::Scalar| {
            let t = clamp_t(x);
            self.axis_value(t, axis).map(|v| v - value)
        };
        let dx = |x: P::Scalar| {
            let t = clamp_t(x);
            self.axis_derivative(t, axis)
        };

        let root = root_newton_raphson(start, fx, dx, eps, max_iter)?;
        Ok(clamp_t(root))
    }
}
