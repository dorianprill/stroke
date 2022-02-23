/// spline.rs
/// Trait for common abstractions over all spline types (Bezier, B-Spline)
use super::Point;

pub trait Spline<P: Point> {
    fn eval(&self, t: P::Scalar) -> P;
}