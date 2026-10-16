//! Lossless conversions between equivalent Bezier representations.

use crate::{Bezier, CubicBezier, LineSegment, Point, QuadraticBezier};

impl<P: Point> From<LineSegment<P>> for Bezier<P, 2> {
    fn from(curve: LineSegment<P>) -> Self {
        Self {
            control_points: curve.into_control_points(),
        }
    }
}

impl<P: Point> From<Bezier<P, 2>> for LineSegment<P> {
    fn from(curve: Bezier<P, 2>) -> Self {
        let [start, end] = curve.into_control_points();
        Self { start, end }
    }
}

impl<P: Point> From<QuadraticBezier<P>> for Bezier<P, 3> {
    fn from(curve: QuadraticBezier<P>) -> Self {
        Self {
            control_points: curve.into_control_points(),
        }
    }
}

impl<P: Point> From<Bezier<P, 3>> for QuadraticBezier<P> {
    fn from(curve: Bezier<P, 3>) -> Self {
        let [start, ctrl, end] = curve.into_control_points();
        Self { start, ctrl, end }
    }
}

impl<P: Point> From<CubicBezier<P>> for Bezier<P, 4> {
    fn from(curve: CubicBezier<P>) -> Self {
        Self {
            control_points: curve.into_control_points(),
        }
    }
}

impl<P: Point> From<Bezier<P, 4>> for CubicBezier<P> {
    fn from(curve: Bezier<P, 4>) -> Self {
        let [start, ctrl1, ctrl2, end] = curve.into_control_points();
        Self {
            start,
            ctrl1,
            ctrl2,
            end,
        }
    }
}
