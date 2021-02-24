use super::*;
//use num_traits::*;
use core::ops::{Add, Sub};

/// Generic floating point number, implemented for f32 and f64
pub trait SizedFloat: 
    Float  
    + Copy 
    + Default 
    + Add<NativeFloat>
    + From<NativeFloat> 
    + Into<NativeFloat>
{} //Signed + Sync + Send + 'static {}

/// Implement for any type that satisfies the bounds
impl<T> SizedFloat for T 
where 
T:  Float 
    + Copy 
    + Default 
    + Add<NativeFloat>
    + From<NativeFloat> 
    + Into<NativeFloat>
{} //Signed + Sync + Send + 'static {}

/// The Point trait is the only interface on which the library relies.
/// The associated constant DIM is necessary so that the memory layout of
/// its implementing type can be made known to the library, whenever new instances are returned.
/// TODO (TBD) specify enclosed type T (any float) of the Point into the trait instead? if so how?
pub trait Point<T>: 
            Add<Self, Output = Self>
            + Sub<Self, Output = Self> 
            + Mul<T, Output = Self>
            + Copy 
            + PartialEq 
            + Default 
            + Sized
            + IntoIterator  
where T: SizedFloat
{
    const DIM: usize;
    // Returns the component of the Point on its axis corresponding to index e.g. [0, 1, 2] -> [x, y, z]
    // TODO remove, use mutable iterator instead (?)
    fn axis(&self, index: usize) -> T;

    // Returns the squared L2-Norm of the Point interpreted as a Vector
    // TODO this could be moved into the library because computability is ensured by its existing trait bounds
    fn squared_length(&self) -> T;
}
