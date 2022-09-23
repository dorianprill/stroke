#[cfg(feature = "nalgebra")]
use num_traits::Float;
use crate::NativeFloat;
use crate::Point;
use super::*;
use core::iter::Sum;


impl<T, const N: usize> Point for nalgebra::Point<T, {N}> 
where
    T: Float
        + Copy
        + Default
        + Sized
        + Add<T, Output = T>
        + Add<NativeFloat, Output = T>
        + Sub<T, Output = T>
        + Sub<NativeFloat, Output = T>
        + Mul<T, Output = T>
        + Mul<NativeFloat, Output = T>
        + Sum<NativeFloat>
        + From<NativeFloat>
        + Into<NativeFloat>,
{
    type Scalar = NativeFloat;
    const DIM: usize = { N };

    fn axis(&self, index: usize) -> Self::Scalar {
        assert!(index <= N);
        self.0[index].into()
    }

    fn squared_length(&self) -> Self::Scalar {
        let mut sqr_dist = 0.0;
        for i in 0..N {
            sqr_dist += (self.0[i] * self.0[i]).into();
        }
        sqr_dist
    }
}
