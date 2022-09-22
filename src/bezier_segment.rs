use crate::bezier::Scalar;

use super::cubic_bezier::CubicBezier;
use super::line::LineSegment;
use super::quadratic_bezier::QuadraticBezier;
use super::*;
use nalgebra::SVector;

#[derive(Copy, Clone, Debug, PartialEq)]
pub enum BezierSegment<T: Scalar, const DIM: usize> {
    Linear(LineSegment<T, DIM>),
    Quadratic(QuadraticBezier<T, DIM>),
    Cubic(CubicBezier<T, DIM>),
}
impl<T: Scalar, const DIM: usize> BezierSegment<T, DIM> {
    pub fn eval<F>(&self, t: T) -> SVector<T, DIM> {
        match self {
            BezierSegment::Linear(segment) => segment.eval(t),
            BezierSegment::Quadratic(segment) => segment.eval(t),
            BezierSegment::Cubic(segment) => segment.eval(t),
        }
    }

    pub fn start(&self) -> SVector<T, DIM> {
        match self {
            BezierSegment::Linear(segment) => segment.start(),
            BezierSegment::Quadratic(segment) => segment.start(),
            BezierSegment::Cubic(segment) => segment.start(),
        }
    }

    #[inline]
    pub fn end(&self) -> SVector<T, DIM> {
        match self {
            BezierSegment::Linear(segment) => segment.end(),
            BezierSegment::Quadratic(segment) => segment.end(),
            BezierSegment::Cubic(segment) => segment.end(),
        }
    }

    #[inline]
    pub fn is_linear<F>(&self, tolerance: F) -> bool
    where
        F: Float + Into<T>,
    {
        match self {
            BezierSegment::Linear(..) => true,
            BezierSegment::Quadratic(segment) => segment.is_linear(tolerance.into()),
            BezierSegment::Cubic(segment) => segment.is_linear(tolerance.into()),
        }
    }

    #[inline]
    pub fn baseline(&self) -> LineSegment<T, DIM> {
        match self {
            BezierSegment::Linear(segment) => *segment,
            BezierSegment::Quadratic(segment) => segment.baseline(),
            BezierSegment::Cubic(segment) => segment.baseline(),
        }
    }

    /// Split this segment into two sub-segments.
    pub fn split<F>(&self, t: F) -> (BezierSegment<T, DIM>, BezierSegment<T, DIM>)
    where
        F: Float + Into<T>,
    {
        match self {
            BezierSegment::Linear(segment) => {
                let (a, b) = segment.split(t.into());
                (BezierSegment::Linear(a), BezierSegment::Linear(b))
            }
            BezierSegment::Quadratic(segment) => {
                let (a, b) = segment.split(t.into());
                (BezierSegment::Quadratic(a), BezierSegment::Quadratic(b))
            }
            BezierSegment::Cubic(segment) => {
                let (a, b) = segment.split(t.into());
                (BezierSegment::Cubic(a), BezierSegment::Cubic(b))
            }
        }
    }
}

impl<T: Scalar, const DIM: usize> From<LineSegment<T, DIM>> for BezierSegment<T, DIM> {
    fn from(s: LineSegment<T, DIM>) -> Self {
        BezierSegment::Linear(s)
    }
}

impl<T: Scalar, const DIM: usize> From<QuadraticBezier<T, DIM>> for BezierSegment<T, DIM> {
    fn from(s: QuadraticBezier<T, DIM>) -> Self {
        BezierSegment::Quadratic(s)
    }
}

impl<T: Scalar, const DIM: usize> From<CubicBezier<T, DIM>> for BezierSegment<T, DIM> {
    fn from(s: CubicBezier<T, DIM>) -> Self {
        BezierSegment::Cubic(s)
    }
}
