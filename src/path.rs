use core::slice;

use crate::bezier_segment::BezierSegment;
use super::*;

/// A path composed of mixed Bezier segments (line/quadratic/cubic).
pub struct BezierPath<P, const N: usize>
where
    P: Point,
{
    segments: ArrayVec<[BezierSegment<P>; N]>,
}

impl<P, const N: usize> BezierPath<P, N>
where
    P: Point,
    [BezierSegment<P>; N]: tinyvec::Array<Item = BezierSegment<P>>,
{
    pub fn new() -> Self {
        BezierPath {
            segments: ArrayVec::new(),
        }
    }

    pub fn len(&self) -> usize {
        self.segments.len()
    }

    pub fn is_empty(&self) -> bool {
        self.segments.len() == 0
    }

    pub fn segments(&self) -> slice::Iter<'_, BezierSegment<P>> {
        self.segments.iter()
    }

    pub fn push(&mut self, segment: BezierSegment<P>) -> bool {
        if self.segments.len() < self.segments.capacity() {
            self.segments.push(segment);
            true
        } else {
            false
        }
    }

    pub fn push_line(&mut self, segment: LineSegment<P>) -> bool {
        self.push(segment.into())
    }

    pub fn push_quadratic(&mut self, segment: QuadraticBezier<P>) -> bool {
        self.push(segment.into())
    }

    pub fn push_cubic(&mut self, segment: CubicBezier<P>) -> bool {
        self.push(segment.into())
    }

    /// Evaluate a point along the path for t in [0,1]. Returns None for empty paths.
    pub fn eval(&self, t: P::Scalar) -> Option<P> {
        let (index, local_t) = self.segment_parameter(t)?;
        Some(self.segments[index].eval(local_t))
    }

    /// Return the bounding box across all segments. Returns None for empty paths.
    pub fn bounding_box(&self) -> Option<[(P::Scalar, P::Scalar); P::DIM]> {
        let mut iter = self.segments.iter();
        let first = match iter.next() {
            Some(segment) => segment.bounding_box(),
            None => return None,
        };

        let mut bounds = first;
        for segment in iter {
            let segment_bounds = segment.bounding_box();
            for dim in 0..P::DIM {
                if segment_bounds[dim].0 < bounds[dim].0 {
                    bounds[dim].0 = segment_bounds[dim].0;
                }
                if segment_bounds[dim].1 > bounds[dim].1 {
                    bounds[dim].1 = segment_bounds[dim].1;
                }
            }
        }

        Some(bounds)
    }

    fn segment_parameter(&self, t: P::Scalar) -> Option<(usize, P::Scalar)> {
        let count = self.segments.len();
        if count == 0 {
            return None;
        }

        let mut t_native: NativeFloat = t.into();
        t_native = t_native.clamp(0.0, 1.0);

        let count_native = count as NativeFloat;
        let scaled = t_native * count_native;
        if scaled >= count_native {
            return Some((count - 1, P::Scalar::from(1.0)));
        }

        let index = scaled.floor() as usize;
        let local = scaled - index as NativeFloat;
        Some((index, P::Scalar::from(local)))
    }
}

impl<P, const N: usize> Default for BezierPath<P, N>
where
    P: Point,
    [BezierSegment<P>; N]: tinyvec::Array<Item = BezierSegment<P>>,
{
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bezier_path_eval_segments() {
        let mut path: BezierPath<PointN<f64, 2>, 4> = BezierPath::new();
        path.push_line(LineSegment::new(
            PointN::new([0.0, 0.0]),
            PointN::new([1.0, 0.0]),
        ));
        path.push_line(LineSegment::new(
            PointN::new([1.0, 0.0]),
            PointN::new([1.0, 1.0]),
        ));

        let p0 = path.eval(0.0).unwrap();
        assert_eq!(p0, PointN::new([0.0, 0.0]));

        let p1 = path.eval(0.25).unwrap();
        assert!((p1 - PointN::new([0.5, 0.0])).squared_length() < EPSILON);

        let p2 = path.eval(0.75).unwrap();
        assert!((p2 - PointN::new([1.0, 0.5])).squared_length() < EPSILON);
    }

    #[test]
    fn bezier_path_bounds_union() {
        let mut path: BezierPath<PointN<f64, 2>, 4> = BezierPath::new();
        path.push_line(LineSegment::new(
            PointN::new([0.0, 0.0]),
            PointN::new([1.0, 0.0]),
        ));
        path.push_line(LineSegment::new(
            PointN::new([1.0, 0.0]),
            PointN::new([1.0, 1.0]),
        ));

        let bounds = path.bounding_box().unwrap();
        assert_eq!(bounds[0], (0.0, 1.0));
        assert_eq!(bounds[1], (0.0, 1.0));
    }

    #[test]
    fn bezier_path_clamps_out_of_range() {
        let mut path: BezierPath<PointN<f64, 2>, 2> = BezierPath::new();
        path.push_line(LineSegment::new(
            PointN::new([0.0, 0.0]),
            PointN::new([1.0, 0.0]),
        ));
        path.push_line(LineSegment::new(
            PointN::new([1.0, 0.0]),
            PointN::new([1.0, 1.0]),
        ));

        let start = path.eval(-1.0).unwrap();
        assert_eq!(start, PointN::new([0.0, 0.0]));

        let end = path.eval(2.0).unwrap();
        assert_eq!(end, PointN::new([1.0, 1.0]));
    }

    #[test]
    fn bezier_path_capacity() {
        let mut path: BezierPath<PointN<f64, 2>, 1> = BezierPath::new();
        let first = path.push_line(LineSegment::new(
            PointN::new([0.0, 0.0]),
            PointN::new([1.0, 0.0]),
        ));
        let second = path.push_line(LineSegment::new(
            PointN::new([1.0, 0.0]),
            PointN::new([2.0, 0.0]),
        ));

        assert!(first);
        assert!(!second);
        assert_eq!(path.len(), 1);
    }

}
