// test/compile with:
// cargo build --features nalgebra
#[cfg(feature = "nalgebra")]
use num_traits::Float;
use crate::NativeFloat;
use crate::Point;
use super::*;
use core::{slice::IterMut, array::IntoIter, iter::Sum};

use nalgebra::Matrix;
use nalgebra::Const;
use nalgebra::ArrayStorage;
use nalgebra::SVector;

// impl<T, const N: usize> Point for nalgebra::Point<T, {N}>
impl<T, const N: usize> Point for SVector<T, {N}> 
where
    T: Float
        + Copy
        + Default
        + Sized
        + Add<T, Output = T>
        + Add<NativeFloat, Output = T>
        + Sub<T, Output = T>
        + Sub<NativeFloat, Output = T>
        //+ Mul<T, Output = T>
        + Mul<NativeFloat, Output = T>
        + Sum<NativeFloat>
        + From<NativeFloat>
        + Into<NativeFloat> + core::ops::AddAssign + core::ops::SubAssign + core::fmt::Debug, Matrix<T, Const<N>, Const<1>, ArrayStorage<T, N, 1>>: core::ops::Mul<f64>
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

    fn iter(&self) -> IntoIter<T, {N}> {
        self.iter()
    }
    fn iter_mut(&self) -> IterMut<Self::Scalar> {
        self.iter_mut()
    }
}


/// Initialize with the Default value for the underlying type
impl<T: Default + Copy, const N: usize> Default for SVector<T, N> {
    fn default() -> Self {
        SVector::from([T::default(); N])
    }
}


