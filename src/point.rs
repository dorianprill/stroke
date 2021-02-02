use core::ops::{Add, Sub};

/// Trait defined over generic 2D points P which themselves are generic over F
/// Many libraries already provide Point-types and the mathematical operations 
/// that we need for working with curves, so that implementing methods requires mostly wrapping.
/// Keeping the trait as minimal as possible to make integration with other libraries easy
pub trait Point: Add + Sub + Copy + PartialEq + Default + IntoIterator      
{
    type Scalar;

    // Returns the component of the Point on its axis corresponding to index e.g. [0, 1, 2] -> [x, y, z]
    // TODO naming: component/axis/at/dim ?
    fn axis(&self, index: usize) -> Self::Scalar;

    // Returns the distance between the two Points self and other
    // TODO maybe squared distance is enough (simpler computation) 
    fn distance(&self, other: Self) -> Self::Scalar;

    // Returns the L2 Norm of the Point interpreted as a Vector
    // TODO naming abs/len/l2
    fn abs(&self) -> Self::Scalar;
}
