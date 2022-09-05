use nalgebra::SVector;
/// spline.rs
/// Trait for common abstractions over all spline types (Bezier, B-Spline)

pub trait Spline<T, const DIM: usize> {
    fn eval(&self, t: T) -> SVector<T, DIM>;
}
