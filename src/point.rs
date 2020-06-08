use core::ops::{Add, Sub};

/// Trait defined over generic 2D points P which themselves are generic over F
/// Many libraries already provide Point-types and the mathematical operations 
/// that we need for working with curves, so that implementing methods requires mostly wrapping.
/// Keeping the trait as minimal as possible to make integration with other libraries easy
pub trait Point: Add + Sub + Copy + PartialEq + Default         
{
    type Scalar;
    fn x(&self) -> Self::Scalar;
    fn y(&self) -> Self::Scalar; // return 0.0 if dim < 2 ?
    // TODO / TBD 
    // - maybe add a z() function for most common use cases, returning 0.0 if dim < 3 ?
    // - define an iterator for the dimensions, would make generic derivatives easier
    //   OR define a dim() method to depend on that ?
    fn distance(&self, other: Self) -> Self::Scalar;
    fn abs(&self) -> Self::Scalar;
}