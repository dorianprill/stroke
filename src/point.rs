use core::ops::{Add, Sub};

/// Trait defined over generic 2D points P which themselves are generic over F
/// Many libraries already provide Point-types and the mathematical operations 
/// that we need for working with curves, so that implementing methods requires mostly wrapping.
/// Keeping the trait as minimal as possible to make integration with other libraries easy
pub trait Point: Add + Sub + Copy + PartialEq + Default + IntoIterator      
{
    type Scalar;
    // TODO
    // define an iterator for the dimensions, would make generic derivatives easier
    // define a dim() method to get the number of dimensions (~ len())
    fn distance(&self, other: Self) -> Self::Scalar;
    fn abs(&self) -> Self::Scalar;
}
