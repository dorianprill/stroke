use super::NativeFloat;
use core::{
    iter::Sum,
    ops::{Add, Div, Mul, Sub},
};
use num_traits::Float;

/// The Point trait is the only interface on which the library relies.
/// The associated constant DIM is necessary so that the memory layout of
/// its implementing type can be made known to the library, whenever new instances are returned.
pub trait Point:
    Add<Self, Output = Self>
    + Sub<Self, Output = Self>
    + Mul<Self::Scalar, Output = Self>
    + Mul<NativeFloat, Output = Self>
    + Copy
    + PartialEq
    //+ PartialOrd
    + Default
    + IntoIterator
{
    type Scalar: Float
        + Default
        + PartialEq
        + From<NativeFloat>
        + Into<NativeFloat>
        + Add<NativeFloat, Output = Self::Scalar>
        + Sub<NativeFloat, Output = Self::Scalar>
        + Mul<NativeFloat, Output = Self::Scalar>
        + Div<NativeFloat, Output = Self::Scalar>
        + Sum<NativeFloat>;
    const DIM: usize;
    // Returns the component of the Point on its axis corresponding to index e.g. [0, 1, 2] -> [x, y, z]
    /// Panics if index is greater than implementors dimension
    // TODO maybe remove in favour of iterator (?)
    fn axis(&self, index: usize) -> Self::Scalar;

    // Returns the squared L2-Norm of the Point interpreted as a Vector
    // TODO this could be moved into the library because computability is ensured by its existing trait bounds
    fn squared_length(&self) -> Self::Scalar;
}
