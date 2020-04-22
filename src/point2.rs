use super::*;
use num_traits::Float;
use super::point::Point;

#[derive(Debug, Copy, Clone)]
pub struct Point2<T>
{
    pub(crate) x: T,
    pub(crate) y: T,
}

impl<T> Point2<T> 
where T: Add 
        + Add<T,Output=T> 
        + Sub 
        + Mul 
        + Mul<T,Output=T> 
        + Clone 
        + Float 
    {
    /// Creates a new Point2<T>, which requires that 
    /// T implements Add, Sub, Mul, and Clone
    pub fn new(x: T, y: T) -> Self {
        Point2 {
            x: x,
            y,
        }
    }

}


impl<T> PartialEq for Point2<T> 
where T: PartialOrd {
    fn eq(&self, other: &Self) -> bool {
        (self.x == other.x) && (self.y == other.y)
    }
}


impl<T> Add for Point2<T>
where
    T: Add<Output=T>,
{
    type Output = Self;

    fn add(self, other: Point2<T>) -> Point2<T> {
        Point2 {
            x: self.x + other.x,
            y: self.y + other.y
        }
    }
}



impl<T> Sub for Point2<T>
where 
    T: Sub<Output=T>
{
    type Output = Self;

    fn sub(self, other: Self) -> Self::Output {
        Point2 {
            x: self.x - other.x,
            y: self.y - other.y,
        }
    }
}

impl<T,U> Mul<U> for Point2<T>
where
    // How you have the mulitplication done is mulitpling T * U => T, this
    // trait bounds for T will specify this requirement as the mul operator is
    // translated to using the first operand as self and the second as rhs. 
    T: Mul<U,Output=T> + Copy,
    U: Clone,
{
    type Output = Point2<T>;

    fn mul(self, _rhs: U) -> Point2<T> {
        return Point2{x: self.x * _rhs.clone(), y: self.y * _rhs}
    }
}



impl Point for Point2<NativeFloat>
{
    type Scalar = NativeFloat;

    fn x(&self) -> Self::Scalar {
        self.x
    }

    fn y(&self) -> Self::Scalar {
        self.y
    }

    /// Returns the distance between self and other
    fn distance(&self, other: Self) -> Self::Scalar {
        ( ((self.x - other.x) * (self.x - other.x))
            + ((self.y - other.y) * (self.y - other.y)) ) .sqrt()
    }

    /// Interprets the Point2 as a vector and returns its norm (distance from origin)
    fn abs(&self) -> Self::Scalar {
        ((self.x * self.x) + (self.y * self.y) ).sqrt()
    }

}