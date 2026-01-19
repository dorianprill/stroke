use core::slice;

use super::*;

/// A path composed of multiple B-Spline segments of the same degree/knots layout.
pub struct BSplinePath<P, const K: usize, const C: usize, const D: usize, const N: usize>
where
    P: Point,
{
    segments: ArrayVec<[BSpline<P, K, C, D>; N]>,
}

impl<P, const K: usize, const C: usize, const D: usize, const N: usize>
    BSplinePath<P, K, C, D, N>
where
    P: Point,
    [BSpline<P, K, C, D>; N]: tinyvec::Array<Item = BSpline<P, K, C, D>>,
{
    pub fn new() -> Self {
        BSplinePath {
            segments: ArrayVec::new(),
        }
    }

    pub fn len(&self) -> usize {
        self.segments.len()
    }

    pub fn is_empty(&self) -> bool {
        self.segments.len() == 0
    }

    pub fn segments(&self) -> slice::Iter<'_, BSpline<P, K, C, D>> {
        self.segments.iter()
    }

    pub fn push(&mut self, segment: BSpline<P, K, C, D>) -> bool {
        if self.segments.len() < self.segments.capacity() {
            self.segments.push(segment);
            true
        } else {
            false
        }
    }

    /// Evaluate a point along the path for t in [0,1]. Returns None for empty paths.
    pub fn eval(&self, t: P::Scalar) -> Option<P>
    where
        [(); D + 1]: Sized,
    {
        let (index, local_t) = self.segment_parameter(t)?;
        let segment = &self.segments[index];
        let (kmin, kmax) = segment.knot_domain();
        let t_mapped = kmin + (kmax - kmin) * local_t;
        segment.eval(t_mapped).ok()
    }

    /// Return a conservative bounding box based on all control points.
    /// Returns None for empty paths.
    pub fn bounding_box(&self) -> Option<[(P::Scalar, P::Scalar); P::DIM]> {
        let mut iter = self.segments.iter();
        let first = match iter.next() {
            Some(segment) => segment_control_bounds(segment),
            None => return None,
        };

        let mut bounds = first;
        for segment in iter {
            let segment_bounds = segment_control_bounds(segment);
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

impl<P, const K: usize, const C: usize, const D: usize, const N: usize> Default
    for BSplinePath<P, K, C, D, N>
where
    P: Point,
    [BSpline<P, K, C, D>; N]: tinyvec::Array<Item = BSpline<P, K, C, D>>,
{
    fn default() -> Self {
        Self::new()
    }
}

fn segment_control_bounds<P, const K: usize, const C: usize, const D: usize>(
    segment: &BSpline<P, K, C, D>,
) -> [(P::Scalar, P::Scalar); P::DIM]
where
    P: Point,
{
    let mut iter = segment.control_points();
    let first = *iter
        .next()
        .expect("BSpline must have at least one control point");
    let mut bounds = [(P::Scalar::default(), P::Scalar::default()); P::DIM];
    for dim in 0..P::DIM {
        let value = first.axis(dim);
        bounds[dim] = (value, value);
    }

    for point in iter {
        for dim in 0..P::DIM {
            let value = point.axis(dim);
            if value < bounds[dim].0 {
                bounds[dim].0 = value;
            }
            if value > bounds[dim].1 {
                bounds[dim].1 = value;
            }
        }
    }

    bounds
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bspline_path_eval_segments() {
        let knots: [f64; 4] = [0.0, 0.0, 1.0, 1.0];
        let s1: BSpline<PointN<f64, 2>, 4, 2, 1> = BSpline::new(
            knots,
            [
                PointN::new([0.0, 0.0]),
                PointN::new([1.0, 0.0]),
            ],
        )
        .unwrap();
        let s2: BSpline<PointN<f64, 2>, 4, 2, 1> = BSpline::new(
            knots,
            [
                PointN::new([1.0, 0.0]),
                PointN::new([1.0, 1.0]),
            ],
        )
        .unwrap();

        let mut path: BSplinePath<PointN<f64, 2>, 4, 2, 1, 4> = BSplinePath::new();
        path.push(s1);
        path.push(s2);

        let p1 = path.eval(0.25).unwrap();
        assert!((p1 - PointN::new([0.5, 0.0])).squared_length() < EPSILON);

        let p2 = path.eval(0.75).unwrap();
        assert!((p2 - PointN::new([1.0, 0.5])).squared_length() < EPSILON);
    }

    #[test]
    fn bspline_path_clamps_out_of_range() {
        let knots: [f64; 4] = [0.0, 0.0, 1.0, 1.0];
        let s1: BSpline<PointN<f64, 2>, 4, 2, 1> = BSpline::new(
            knots,
            [
                PointN::new([0.0, 0.0]),
                PointN::new([1.0, 0.0]),
            ],
        )
        .unwrap();
        let s2: BSpline<PointN<f64, 2>, 4, 2, 1> = BSpline::new(
            knots,
            [
                PointN::new([1.0, 0.0]),
                PointN::new([1.0, 1.0]),
            ],
        )
        .unwrap();

        let mut path: BSplinePath<PointN<f64, 2>, 4, 2, 1, 4> = BSplinePath::new();
        path.push(s1);
        path.push(s2);

        let start = path.eval(-1.0).unwrap();
        assert_eq!(start, PointN::new([0.0, 0.0]));

        let end = path.eval(2.0).unwrap();
        assert_eq!(end, PointN::new([1.0, 1.0]));
    }

    #[test]
    fn bspline_path_bounds_control_points() {
        let knots: [f64; 4] = [0.0, 0.0, 1.0, 1.0];
        let s1: BSpline<PointN<f64, 2>, 4, 2, 1> = BSpline::new(
            knots,
            [
                PointN::new([0.0, 2.0]),
                PointN::new([1.0, -1.0]),
            ],
        )
        .unwrap();
        let s2: BSpline<PointN<f64, 2>, 4, 2, 1> = BSpline::new(
            knots,
            [
                PointN::new([-1.0, 0.5]),
                PointN::new([2.0, -2.0]),
            ],
        )
        .unwrap();

        let mut path: BSplinePath<PointN<f64, 2>, 4, 2, 1, 4> = BSplinePath::new();
        path.push(s1);
        path.push(s2);

        let bounds = path.bounding_box().unwrap();
        assert_eq!(bounds[0], (-1.0, 2.0));
        assert_eq!(bounds[1], (-2.0, 2.0));
    }
}
