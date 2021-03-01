use super::Point;
use super::*;

/// LineSegment defined by a start and an endpoint, evaluatable
/// anywhere inbetween using interpolation parameter t: [0,1] in eval()
/// A LineSegment is equal to a linear Bezier curve, which is why there is no
/// specialized type for that case.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct LineSegment<P> {
    pub(crate) start: P,
    pub(crate) end: P,
}

impl<P> LineSegment<P>
where
    P: Point,
{
    pub fn new(start: P, end: P) -> Self {
        LineSegment { start, end }
    }

    pub fn eval(&self, t: P::Scalar) -> P {
        return self.start + (self.end - self.start) * t;
    }

    pub fn split(&self, t: P::Scalar) -> (Self, Self) {
        // compute the split point by interpolation
        let ctrl_ab = self.start + (self.start - self.end) * t;

        return (
            LineSegment {
                start: self.start,
                end: ctrl_ab,
            },
            LineSegment {
                start: ctrl_ab,
                end: self.end,
            },
        );
    }

    // DEPRECATED
    // pub fn to_line(&self) -> Line<P> {
    //     Line {
    //         origin: self.start,
    //         vector: self.end - self.start,
    //     }
    // }

    /// Return the distance from the LineSegment to Point p by calculating the projection
    pub fn distance_to_point(&self, p: P) -> P::Scalar {
        let l2 = (self.end - self.start).squared_length();
        // if start and endpoint are approx the same, return the distance to either
        if l2 < P::Scalar::from(EPSILON) {
            return (self.start - p).squared_length().sqrt();
        } else {
            let v1 = p - self.start;
            let v2 = self.end - self.start;
            let mut dot = P::Scalar::from(0.0);
            for (i, _) in v1.into_iter().enumerate() {
                dot = dot + v1.axis(i) * v2.axis(i);
            }
            let mut t = P::Scalar::from(0.0);
            if dot / l2 < P::Scalar::from(1.0) {
                t = dot / l2;
            }
            if t < P::Scalar::from(0.0) {
                t = P::Scalar::from(0.0);
            }
            let projection = self.start + (self.end - self.start) * t; // Projection falls on the segment
            return (p - projection).squared_length().sqrt();
        }
    }

    /// Sample the coordinate axis of the segment at t (expecting t between 0 and 1).
    pub fn axis(&self, t: P::Scalar, axis: usize) -> P::Scalar {
        self.start.axis(axis) + (self.end.axis(axis) - self.start.axis(axis)) * t
    }

    /// Return the derivative function.
    /// The derivative is also a bezier curve but of degree n-1 - In the case of a line the derivative vector.
    pub fn derivative(&self) -> P {
        return self.end - self.start;
    }

    pub(crate) fn root(&self, a: P::Scalar, b: P::Scalar) -> ArrayVec<[P::Scalar; 1]> {
        let mut r = ArrayVec::new();
        if a.abs() < EPSILON.into() {
            return r;
        }
        r.push(-b / a);
        return r;
    }

    /// Return the bounding box of the line as an array of (min, max) tuples for each dimension (its index)
    pub fn bounding_box(&self) -> [(P::Scalar, P::Scalar); P::DIM] {
        let mut bounds = [(P::Scalar::default(), P::Scalar::default()); P::DIM];

        // find min/max for that particular axis
        // TODO shoul be rewritten once 'Iterator' is implemented on P to get rid of .axis() method
        for (i, _) in self.start.into_iter().enumerate() {
            if self.start.axis(i) < self.end.axis(i) {
                bounds[i] = (self.start.axis(i), self.end.axis(i));
            } else {
                bounds[i] = (self.end.axis(i), self.start.axis(i));
            }
        }

        return bounds;
    }
}

#[cfg(test)]
mod tests {
    use super::point_generic::PointN;
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
        assert!((mid - line.start).squared_length() - (mid - line.end).squared_length() < EPSILON)
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
        assert!(line.distance_to_point(p1) - 4.0.abs() < EPSILON);
        assert!(((p1 - line.start).squared_length().sqrt() - 4.0).abs() < EPSILON);
        assert!(((p1 - line.end).squared_length().sqrt() - 5.0).abs() < EPSILON);
        // dist to midpoint (t=0.5) should be 1
        let p2 = PointN::new([1.5f64, 2f64, 0f64]);
        assert!(((p2 - line.eval(0.5)).squared_length().sqrt() - 1.0).abs() < EPSILON);
    }
}
