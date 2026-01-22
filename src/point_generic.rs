//! Const-generic point type for users who do not want external math dependencies.

use core::iter::IntoIterator;
use core::ops::{Add, Index, IndexMut, Mul, Sub};
use core::slice;

use super::Point;
use num_traits::Float;

/// Point with dimensions of constant generic size N and of generic type T
///
/// (Implemented as Newtype Pattern on an array
/// see book or <https://www.worthe-it.co.za/blog/2020-10-31-newtype-pattern-in-rust.html>)
/// This type implements [`Point`] and `Index<usize>` so it
/// can be used with component-aware helpers like bounding boxes.
#[derive(Debug, Copy, Clone)]
pub struct PointN<T, const N: usize>([T; N]);

impl<T, const N: usize> PointN<T, N> {
    /// Construct a point from a fixed-size array of components.
    pub fn new(array: [T; N]) -> Self {
        PointN(array)
    }
}

/// Initialize with the Default value for the underlying type
impl<T: Default + Copy, const N: usize> Default for PointN<T, N> {
    fn default() -> Self {
        PointN([T::default(); N])
    }
}

impl<T, const N: usize> PartialEq for PointN<T, N>
where
    T: PartialEq,
{
    fn eq(&self, other: &Self) -> bool {
        for i in 0..N {
            if self.0[i] != other.0[i] {
                return false;
            }
        }
        true
    }
}

impl<T, const N: usize> Add for PointN<T, N>
where
    T: Add<Output = T> + Copy,
{
    type Output = Self;

    fn add(self, other: PointN<T, N>) -> PointN<T, N> {
        let mut res = self;
        for i in 0..N {
            res.0[i] = self.0[i] + other.0[i];
        }
        res
    }
}

/// This is not required by the Point trait or library but
/// convenient if you want to use the type externally
impl<T, const N: usize> Add<T> for PointN<T, N>
where
    T: Add<Output = T> + Copy,
{
    type Output = Self;

    fn add(self, _rhs: T) -> PointN<T, N> {
        let mut res = self;
        for i in 0..N {
            res.0[i] = self.0[i] + _rhs;
        }
        res
    }
}

impl<T, const N: usize> Sub for PointN<T, N>
where
    T: Sub<Output = T> + Copy,
{
    type Output = Self;

    fn sub(self, other: PointN<T, N>) -> PointN<T, N> {
        let mut res = self;
        for i in 0..N {
            res.0[i] = self.0[i] - other.0[i];
        }
        res
    }
}

/// This is not required by the Point trait or library but
/// convenient if you want to use the type externally
impl<T, const N: usize> Sub<T> for PointN<T, N>
where
    T: Sub<Output = T> + Copy,
{
    type Output = Self;

    fn sub(self, _rhs: T) -> PointN<T, N> {
        let mut res = self;
        for i in 0..N {
            res.0[i] = self.0[i] - _rhs;
        }
        res
    }
}

impl<T, const N: usize, U> Mul<U> for PointN<T, N>
where
    // The mulitplication is done by mulitpling T * U => T, this
    // trait bound for T will specify this requirement as the mul operator is
    // translated to using the first operand as self and the second as rhs.
    T: Mul<U, Output = T> + Copy, //+ SizedFloat,
    U: Clone + Copy,
{
    type Output = PointN<T, N>;

    fn mul(self, _rhs: U) -> PointN<T, N> {
        let mut res = self;
        for i in 0..res.0.len() {
            res.0[i] = res.0[i] * _rhs;
        }
        res
    }
}

impl<T, const N: usize> IntoIterator for PointN<T, N> {
    type Item = T;
    type IntoIter = core::array::IntoIter<Self::Item, N>;

    fn into_iter(self) -> Self::IntoIter {
        IntoIterator::into_iter(self.0)
    }
}

impl<'a, T, const N: usize> IntoIterator for &'a mut PointN<T, N> {
    type Item = &'a mut T;
    type IntoIter = slice::IterMut<'a, T>;

    fn into_iter(self) -> slice::IterMut<'a, T> {
        self.0.iter_mut()
    }
}

impl<T, const N: usize> Index<usize> for PointN<T, N> {
    type Output = T;

    fn index(&self, index: usize) -> &Self::Output {
        &self.0[index]
    }
}

impl<T, const N: usize> IndexMut<usize> for PointN<T, N> {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.0[index]
    }
}

impl<T, const N: usize> Point for PointN<T, N>
where
    T: Float,
{
    type Scalar = T;
    const DIM: usize = N;
}
