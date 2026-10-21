//! Array-backed points with no external math dependency.

use core::iter::IntoIterator;
use core::ops::{Add, Index, IndexMut, Mul, Sub};
use core::slice;

use super::Point;
use num_traits::Float;

/// A point backed by `[T; N]`.
/// Implements [`Point`] and coordinate indexing, including the bounds methods
/// that require `PointIndex`.
#[derive(Debug, Copy, Clone, PartialEq)]
pub struct ArrayPoint<T, const N: usize>([T; N]);

impl<T, const N: usize> ArrayPoint<T, N> {
    /// Construct a point from a fixed-size array of components.
    pub fn new(array: [T; N]) -> Self {
        ArrayPoint(array)
    }

    /// Borrow the component array.
    pub fn as_array(&self) -> &[T; N] {
        &self.0
    }

    /// Consume the point, returning its component array without copying.
    pub fn into_array(self) -> [T; N] {
        self.0
    }
}

impl<T, const N: usize> From<[T; N]> for ArrayPoint<T, N> {
    fn from(array: [T; N]) -> Self {
        Self::new(array)
    }
}

impl<T, const N: usize> From<ArrayPoint<T, N>> for [T; N] {
    fn from(point: ArrayPoint<T, N>) -> Self {
        point.into_array()
    }
}

impl<T, const N: usize> AsRef<[T; N]> for ArrayPoint<T, N> {
    fn as_ref(&self) -> &[T; N] {
        self.as_array()
    }
}

impl<'a, T, const N: usize> IntoIterator for &'a ArrayPoint<T, N> {
    type Item = &'a T;
    type IntoIter = slice::Iter<'a, T>;
    fn into_iter(self) -> Self::IntoIter {
        self.0.iter()
    }
}

/// Fill the coordinates with `T::default()`.
impl<T: Default + Copy, const N: usize> Default for ArrayPoint<T, N> {
    fn default() -> Self {
        ArrayPoint([T::default(); N])
    }
}

impl<T, const N: usize> Add for ArrayPoint<T, N>
where
    T: Add<Output = T> + Copy,
{
    type Output = Self;

    fn add(self, other: ArrayPoint<T, N>) -> ArrayPoint<T, N> {
        let mut res = self;
        for (value, rhs) in res.0.iter_mut().zip(other.0) {
            *value = *value + rhs;
        }
        res
    }
}

/// Add a scalar to each coordinate.
impl<T, const N: usize> Add<T> for ArrayPoint<T, N>
where
    T: Add<Output = T> + Copy,
{
    type Output = Self;

    fn add(self, rhs: T) -> ArrayPoint<T, N> {
        Self(self.0.map(|value| value + rhs))
    }
}

impl<T, const N: usize> Sub for ArrayPoint<T, N>
where
    T: Sub<Output = T> + Copy,
{
    type Output = Self;

    fn sub(self, other: ArrayPoint<T, N>) -> ArrayPoint<T, N> {
        let mut res = self;
        for (value, rhs) in res.0.iter_mut().zip(other.0) {
            *value = *value - rhs;
        }
        res
    }
}

/// Subtract a scalar from each coordinate.
impl<T, const N: usize> Sub<T> for ArrayPoint<T, N>
where
    T: Sub<Output = T> + Copy,
{
    type Output = Self;

    fn sub(self, rhs: T) -> ArrayPoint<T, N> {
        Self(self.0.map(|value| value - rhs))
    }
}

impl<T, const N: usize, U> Mul<U> for ArrayPoint<T, N>
where
    T: Mul<U, Output = T>,
    U: Copy,
{
    type Output = ArrayPoint<T, N>;

    fn mul(self, rhs: U) -> ArrayPoint<T, N> {
        Self(self.0.map(|value| value * rhs))
    }
}

impl<T, const N: usize> IntoIterator for ArrayPoint<T, N> {
    type Item = T;
    type IntoIter = core::array::IntoIter<Self::Item, N>;

    fn into_iter(self) -> Self::IntoIter {
        IntoIterator::into_iter(self.0)
    }
}

impl<'a, T, const N: usize> IntoIterator for &'a mut ArrayPoint<T, N> {
    type Item = &'a mut T;
    type IntoIter = slice::IterMut<'a, T>;

    fn into_iter(self) -> slice::IterMut<'a, T> {
        self.0.iter_mut()
    }
}

impl<T, const N: usize> Index<usize> for ArrayPoint<T, N> {
    type Output = T;

    fn index(&self, index: usize) -> &Self::Output {
        &self.0[index]
    }
}

impl<T, const N: usize> IndexMut<usize> for ArrayPoint<T, N> {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.0[index]
    }
}

impl<T, const N: usize> Point for ArrayPoint<T, N>
where
    T: Float,
{
    type Scalar = T;
    const DIM: usize = N;

    fn is_finite(&self) -> bool {
        self.0.iter().all(|value| value.is_finite())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn point_finiteness_checks_every_coordinate() {
        assert!(ArrayPoint::new([0.0_f64, -0.0, f64::MIN_POSITIVE]).is_finite());
        assert!(!ArrayPoint::new([0.0_f64, f64::NAN]).is_finite());
        assert!(!ArrayPoint::new([f64::NEG_INFINITY, 0.0]).is_finite());
    }
}
