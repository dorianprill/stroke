//! Common spline abstraction for evaluatable curve types.
use super::Point;

/// Trait implemented by spline types that can be evaluated at a parameter `t`.
pub trait Spline<P: Point> {
    /// Evaluate the spline at parameter `t`.
    fn eval(&self, t: P::Scalar) -> P;
}
