use core::ops::{Add, Sub};
/// The Point trait is the only interface on which the library relies.
/// The associated constant DIM is necessary so that the memory layout of
/// its implementing type can be made known to the library, whenever new instances are returned.
/// TODO (TBD) specify enclosed type T (any float) of the Point into the trait instead? if so how?
pub trait Point: Add + Sub + Copy + PartialEq + Default + IntoIterator      
{
    type Scalar;
    const DIM: usize;
    // Returns the component of the Point on its axis corresponding to index e.g. [0, 1, 2] -> [x, y, z]
    // TODO remove, use mutable iterator instead (?)
    fn axis(&self, index: usize) -> Self::Scalar;

    // Returns the squared L2-Norm of the Point interpreted as a Vector
    // TODO this could be moved into the library because computability is ensured by its existing trait bounds
    fn squared_length(&self) -> Self::Scalar;
}
