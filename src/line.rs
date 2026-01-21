//! Line segment curve type.

use num_traits::{Float, NumCast};
use super::{ArrayVec, Point, PointDot, PointIndex, PointNorm};

/// LineSegment defined by a start and an endpoint, evaluable anywhere inbetween using interpolation parameter t: [0,1] in eval().
///
/// A LineSegment is equal to a linear Bezier curve, which is why there is no
/// specialized type for that case.
///
/// Methods that need component access or norms add `PointIndex`/`PointNorm`
/// bounds as required.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct LineSegment<P> {
    pub(crate) start: P,
    pub(crate) end: P,
}

impl<P> LineSegment<P>
where
    P: Point,
{
    /// Create a new line segment from `start` to `end`.
    pub fn new(start: P, end: P) -> Self {
        LineSegment { start, end }
    }

    /// Evaluate a point along the segment for `t` in `[0, 1]`.
    pub fn eval(&self, t: P::Scalar) -> P {
        self.start + (self.end - self.start) * t
    }

    /// Split the segment at `t` into two sub-segments.
    pub fn split(&self, t: P::Scalar) -> (Self, Self) {
        // compute the split point by interpolation
        let ctrl_ab = self.start + (self.end - self.start) * t;

        (
            LineSegment {
                start: self.start,
                end: ctrl_ab,
            },
            LineSegment {
                start: ctrl_ab,
                end: self.end,
            },
        )
    }

    /// Return the distance from the LineSegment to Point p by calculating the projection
    pub fn distance_to_point(&self, p: P) -> P::Scalar
    where
        P: PointNorm,
    {
        let l2 = (self.end - self.start).squared_norm();
        // if start and endpoint are approx the same, return the distance to either
        if l2 < P::Scalar::epsilon() {
            (self.start - p).squared_norm().sqrt()
        } else {
            let v1 = p - self.start;
            let v2 = self.end - self.start;
            let dot = PointDot::dot(&v1, &v2);
            // v1 and v2 will by definition always have the same number of axes and produce a value for each Item
            // dot = v1.into_iter()
            //         .zip(v2.into_iter())
            //         .map(|(x1, x2)| x1 * x2)
            //         .sum::<P::Scalar>();
            let t = (dot / l2).clamp(
                <P::Scalar as NumCast>::from(0.0).unwrap(),
                <P::Scalar as NumCast>::from(1.0).unwrap(),
            );
            let projection = self.start + (self.end - self.start) * t; // Projection falls on the segment

            (p - projection).squared_norm().sqrt()
        }
    }

    /// Return the derivative vector of the segment.
    pub fn derivative(&self) -> P {
        self.end - self.start
    }

    pub(crate) fn root(&self, a: P::Scalar, b: P::Scalar) -> ArrayVec<[P::Scalar; 1]>
    where
        [P::Scalar; 1]: tinyvec::Array<Item = P::Scalar>,
    {
        let mut r = ArrayVec::new();
        if a.abs() < P::Scalar::epsilon() {
            return r;
        }
        r.push(-b / a);
        r
    }

    /// Return the bounding box of the line as an array of (min, max) tuples for each dimension (its index)
    pub fn bounding_box(&self) -> [(P::Scalar, P::Scalar); P::DIM]
    where
        P: PointIndex,
    {
        let mut bounds = [(<P::Scalar as NumCast>::from(0.0).unwrap(), <P::Scalar as NumCast>::from(0.0).unwrap()); P::DIM];

        // find min/max for each axis
        for i in 0..P::DIM {
            if self.start[i] < self.end[i] {
                bounds[i] = (self.start[i], self.end[i]);
            } else {
                bounds[i] = (self.end[i], self.start[i]);
            }
        }

        bounds
    }
}

#[cfg(test)]
mod tests {
    use crate::{PointN, EPSILON};
    use super::*;
    /// Check whether a line segment interpolation p + t*(q-p) at t=0.5
    /// yields equal distance to the start (p)/end (q) points (up to machine accuracy).
    #[test]
    fn line_segment_interpolation() {
        let line = LineSegment {
            start: PointN::new([0f64, 1.77f64]),
            end: PointN::new([4.3f64, 3f64]),
        };

        let mid = line.eval(0.5);
        let diff =
            ((mid - line.start).squared_norm() - (mid - line.end).squared_norm()).abs();
        assert!(diff < EPSILON * 10.0);
    }

    /// Check whether classic pythagorean equality holds for sides 3, 4 with hypothenuse 5
    #[test]
    fn line_segment_distance_to_point() {
        // 3D cause why not
        let line = LineSegment {
            start: PointN::new([0f64, 1f64, 0f64]),
            end: PointN::new([3f64, 1f64, 0f64]),
        };
        // dist to start should be 4; dist to end should be 5
        let p1 = PointN::new([0f64, 5f64, 0f64]);
        assert!((line.distance_to_point(p1) - 4.0).abs() < EPSILON);
        assert!(((p1 - line.start).squared_norm().sqrt() - 4.0).abs() < EPSILON);
        assert!(((p1 - line.end).squared_norm().sqrt() - 5.0).abs() < EPSILON);
        // dist to midpoint (t=0.5) should be 1
        let p2 = PointN::new([1.5f64, 2f64, 0f64]);
        assert!(((p2 - line.eval(0.5)).squared_norm().sqrt() - 1.0).abs() < EPSILON);
    }

    #[test]
    fn line_segment_split_midpoint() {
        let line = LineSegment {
            start: PointN::new([0f64, 0f64]),
            end: PointN::new([4f64, 0f64]),
        };

        let (left, right) = line.split(0.5);
        assert_eq!(left.start, PointN::new([0.0, 0.0]));
        assert_eq!(left.end, PointN::new([2.0, 0.0]));
        assert_eq!(right.start, PointN::new([2.0, 0.0]));
        assert_eq!(right.end, PointN::new([4.0, 0.0]));
    }

    #[test]
    fn line_segment_distance_to_point_after_end() {
        let line = LineSegment {
            start: PointN::new([0f64, 0f64]),
            end: PointN::new([1f64, 0f64]),
        };
        let p = PointN::new([2f64, 0f64]);
        assert!((line.distance_to_point(p) - 1.0).abs() < EPSILON);
    }
}
