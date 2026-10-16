//! Line segment curve type.

use super::{ArrayVec, Point, PointDot, PointIndex, PointNorm};
use crate::error::{CurveConstructionError, validate_control_points};
use crate::numeric::is_negligible;
use num_traits::{Float, NumCast};

const DEFAULT_LENGTH_STEPS: usize = 64;

/// `LineSegment` defined by a start and an endpoint, evaluable anywhere inbetween using interpolation parameter t: \[0, 1\] in `eval()`.
///
/// A `LineSegment` is equal to a linear Bezier curve, which is why there is no
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
    ///
    /// # Errors
    ///
    /// Returns [`CurveConstructionError::NonFiniteControlPoint`] with the index
    /// of the first control point containing a non-finite coordinate.
    pub fn new(start: P, end: P) -> Result<Self, CurveConstructionError> {
        validate_control_points(&[start, end])?;
        Ok(LineSegment { start, end })
    }

    /// Return the start point of the segment.
    pub fn start(&self) -> P {
        let zero = <P::Scalar as NumCast>::from(0.0).unwrap();
        self.eval(zero)
    }

    /// Return the end point of the segment.
    pub fn end(&self) -> P {
        let one = <P::Scalar as NumCast>::from(1.0).unwrap();
        self.eval(one)
    }

    /// Return a segment with reversed direction.
    pub fn reverse(&self) -> Self {
        LineSegment {
            start: self.end,
            end: self.start,
        }
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

    /// Return the distance from the `LineSegment` to Point p by calculating the projection
    pub fn distance_to_point(&self, p: P) -> P::Scalar
    where
        P: PointNorm,
    {
        let l2 = (self.end - self.start).squared_norm();
        let zero = <P::Scalar as NumCast>::from(0.0).unwrap();
        // A nonzero segment remains meaningful regardless of its coordinate scale.
        if l2 == zero {
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

    /// Return the unit tangent direction at `t`.
    pub fn tangent(&self, _t: P::Scalar) -> P
    where
        P: PointNorm,
    {
        let one = <P::Scalar as NumCast>::from(1.0).unwrap();
        let zero = <P::Scalar as NumCast>::from(0.0).unwrap();
        let dir = self.derivative();
        let len = dir.squared_norm().sqrt();
        if len == zero { dir } else { dir * (one / len) }
    }

    /// Return the curvature magnitude at `t`.
    ///
    /// Lines have zero curvature, so this always returns `0`.
    pub fn curvature(&self, _t: P::Scalar) -> P::Scalar
    where
        P: PointNorm,
    {
        <P::Scalar as NumCast>::from(0.0).unwrap()
    }

    /// Return the principal normal direction at `t`.
    ///
    /// Lines have zero curvature, so this always returns `None`.
    pub fn normal(&self, _t: P::Scalar) -> Option<P>
    where
        P: PointNorm,
    {
        None
    }

    /// Approximate parameter `t` at arc length `s`.
    pub fn t_at_length_approx(&self, s: P::Scalar, _nsteps: usize) -> P::Scalar
    where
        P: PointNorm,
    {
        let zero = <P::Scalar as NumCast>::from(0.0).unwrap();
        let one = <P::Scalar as NumCast>::from(1.0).unwrap();
        let length = (self.end - self.start).squared_norm().sqrt();
        if length == zero {
            zero
        } else {
            (s / length).clamp(zero, one)
        }
    }

    /// Approximate parameter `t` at arc length `s` using a default resolution.
    pub fn t_at_length(&self, s: P::Scalar) -> P::Scalar
    where
        P: PointNorm,
    {
        self.t_at_length_approx(s, DEFAULT_LENGTH_STEPS)
    }

    /// Evaluate the point at arc length `s`.
    pub fn point_at_length_approx(&self, s: P::Scalar, nsteps: usize) -> P
    where
        P: PointNorm,
    {
        self.eval(self.t_at_length_approx(s, nsteps))
    }

    /// Evaluate the point at arc length `s` using a default resolution.
    pub fn point_at_length(&self, s: P::Scalar) -> P
    where
        P: PointNorm,
    {
        self.point_at_length_approx(s, DEFAULT_LENGTH_STEPS)
    }

    pub(crate) fn root(&self, a: P::Scalar, b: P::Scalar) -> ArrayVec<[P::Scalar; 1]>
    where
        [P::Scalar; 1]: tinyvec::Array<Item = P::Scalar>,
    {
        let mut r = ArrayVec::new();
        if is_negligible(a, a.abs().max(b.abs())) {
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
        let mut bounds = [(
            <P::Scalar as NumCast>::from(0.0).unwrap(),
            <P::Scalar as NumCast>::from(0.0).unwrap(),
        ); P::DIM];

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
    use super::*;
    use crate::CurveConstructionError;
    use crate::PointN;
    use crate::test_support::{assert_approx_equal, assert_point_approx_equal};

    const ABSOLUTE_TOLERANCE: f64 = 1e-12;
    const RELATIVE_TOLERANCE: f64 = 1e-12;

    #[test]
    fn rejects_non_finite_control_points() {
        let result = LineSegment::new(
            PointN::new([0.0_f64, 0.0]),
            PointN::new([f64::INFINITY, 0.0]),
        );

        assert_eq!(
            result,
            Err(CurveConstructionError::NonFiniteControlPoint { index: 1 })
        );
    }

    /// Check whether a line segment interpolation p + t*(q-p) at t=0.5
    /// yields equal distance to the start (p)/end (q) points (up to machine accuracy).
    #[test]
    fn line_segment_interpolation() {
        let line = LineSegment {
            start: PointN::new([0f64, 1.77f64]),
            end: PointN::new([4.3f64, 3f64]),
        };

        let mid = line.eval(0.5);
        assert_approx_equal(
            (mid - line.start).squared_norm(),
            (mid - line.end).squared_norm(),
            ABSOLUTE_TOLERANCE,
            RELATIVE_TOLERANCE,
        );
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
        assert_approx_equal(
            line.distance_to_point(p1),
            4.0,
            ABSOLUTE_TOLERANCE,
            RELATIVE_TOLERANCE,
        );
        assert_approx_equal(
            (p1 - line.start).squared_norm().sqrt(),
            4.0,
            ABSOLUTE_TOLERANCE,
            RELATIVE_TOLERANCE,
        );
        assert_approx_equal(
            (p1 - line.end).squared_norm().sqrt(),
            5.0,
            ABSOLUTE_TOLERANCE,
            RELATIVE_TOLERANCE,
        );
        // dist to midpoint (t=0.5) should be 1
        let p2 = PointN::new([1.5f64, 2f64, 0f64]);
        assert_approx_equal(
            (p2 - line.eval(0.5)).squared_norm().sqrt(),
            1.0,
            ABSOLUTE_TOLERANCE,
            RELATIVE_TOLERANCE,
        );
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
        assert_approx_equal(
            line.distance_to_point(p),
            1.0,
            ABSOLUTE_TOLERANCE,
            RELATIVE_TOLERANCE,
        );
    }

    #[test]
    fn small_scale_line_remains_non_degenerate() {
        let scale = 1e-20_f64;
        let line = LineSegment::new(PointN::new([0.0, 0.0]), PointN::new([scale, 0.0])).unwrap();

        assert_point_approx_equal(
            line.tangent(0.5),
            PointN::new([1.0, 0.0]),
            ABSOLUTE_TOLERANCE,
            RELATIVE_TOLERANCE,
        );
        assert_approx_equal(
            line.t_at_length(scale * 0.5),
            0.5,
            ABSOLUTE_TOLERANCE,
            RELATIVE_TOLERANCE,
        );
        assert_approx_equal(
            line.distance_to_point(PointN::new([scale * 0.5, scale])),
            scale,
            1e-32,
            1e-12,
        );
    }

    #[test]
    fn line_root_is_scale_invariant() {
        let line = LineSegment::new(PointN::new([0.0_f64]), PointN::new([1.0])).unwrap();
        let roots = line.root(1e-30, -5e-31);

        assert_eq!(roots.len(), 1);
        assert_approx_equal(roots[0], 0.5, ABSOLUTE_TOLERANCE, RELATIVE_TOLERANCE);
    }

    #[test]
    fn line_segment_api_parity() {
        let line = LineSegment::new(PointN::new([0.0_f64, 0.0]), PointN::new([2.0, 0.0])).unwrap();

        assert_eq!(line.start(), line.eval(0.0));
        assert_eq!(line.end(), line.eval(1.0));

        let reversed = line.reverse();
        assert_eq!(reversed.start(), line.end());
        assert_eq!(reversed.end(), line.start());

        let tangent = line.tangent(0.3);
        assert_approx_equal(tangent[0], 1.0, ABSOLUTE_TOLERANCE, RELATIVE_TOLERANCE);
        assert_approx_equal(tangent[1], 0.0, ABSOLUTE_TOLERANCE, RELATIVE_TOLERANCE);

        let curvature = line.curvature(0.3);
        assert_approx_equal(curvature, 0.0, ABSOLUTE_TOLERANCE, RELATIVE_TOLERANCE);
        assert!(line.normal(0.3).is_none());

        let t = line.t_at_length(1.0);
        assert_approx_equal(t, 0.5, ABSOLUTE_TOLERANCE, RELATIVE_TOLERANCE);

        let p = line.point_at_length(1.0);
        assert_point_approx_equal(
            p,
            PointN::new([1.0, 0.0]),
            ABSOLUTE_TOLERANCE,
            RELATIVE_TOLERANCE,
        );
    }
}
