use alloc::vec::Vec;
use core::fmt::Debug;
use core::iter::{from_fn, IntoIterator};
use core::slice;
use itertools::Itertools;
use nalgebra::{ClosedMul, ComplexField, RealField, SVector};
use num_traits::{NumCast, Zero};
use smallvec::SmallVec;

//use crate::roots::RootFindingError;

use super::*;
use crate::spline::Spline;

//pub trait Scalar : nalgebra::Scalar + Clone + PartialOrd + SimdValue<Element = Self> + SimdPartialOrd + ClosedAdd + ClosedDiv + ClosedMul + ClosedSub + Clone + Field {}

pub trait Scalar: RealField + Copy + NumCast + Default + ClosedMul + Debug + Float {}

impl<T: RealField + Copy + NumCast + Default + ClosedMul + Debug + Float> Scalar for T {}

/// General implementation of a Bezier curve of arbitrary degree (= number of control points - 1).
/// The curve is solely defined by an array of 'control_points'. The degree is defined as degree = control_points.len() - 1.
/// Points on the curve can be evaluated with an interpolation parameter 't' in interval [0,1] using the eval() and eval_casteljau() methods.
/// Generic parameters:
/// T: Generic points 'T' as defined by there Point trait
/// const generic parameters:
/// N: Number of control points
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct Bezier<T: Scalar, const DIM: usize, const N: usize> {
    /// Control points which define the curve and hence its degree
    pub control_points: [SVector<T, DIM>; N],
}

impl<T: Scalar, const DIM: usize, const N: usize> Spline<T, DIM> for Bezier<T, DIM, N> {
    fn eval(&self, t: T) -> SVector<T, DIM> {
        self.eval(t)
    }
}

impl<T: Scalar, const DIM: usize, const N: usize> IntoIterator for Bezier<T, DIM, N> {
    type Item = SVector<T, DIM>;
    type IntoIter = core::array::IntoIter<Self::Item, N>;

    fn into_iter(self) -> Self::IntoIter {
        IntoIterator::into_iter(self.control_points)
    }
}
impl<'a, T: Scalar, const DIM: usize, const N: usize> IntoIterator for &'a mut Bezier<T, DIM, N> {
    type Item = &'a mut SVector<T, DIM>;
    type IntoIter = slice::IterMut<'a, SVector<T, DIM>>;

    fn into_iter(self) -> slice::IterMut<'a, SVector<T, DIM>> {
        self.control_points.iter_mut()
    }
}

impl<T: Scalar, const DIM: usize, const N: usize> Bezier<T, DIM, N> {
    /// Create a new Bezier curve that interpolates the `control_points`. The degree is defined as degree = control_points.len() - 1.
    /// Desired curve must have a valid number of control points and knots in relation to its degree or the constructor will return None.
    /// A B-Spline curve requires at least one more control point than the degree (`control_points.len() >
    /// degree`) and the number of knots should be equal to `control_points.len() + degree + 1`.
    pub fn new(control_points: [SVector<T, DIM>; N]) -> Self {
        Bezier { control_points }
    }

    pub fn control_points(&self) -> [SVector<T, DIM>; N] {
        self.control_points.clone()
    }

    /// Evaluate a point on the curve at point 't' which should be in the interval [0,1]
    /// This is implemented using De Casteljau's algorithm (over a temporary array with const generic sizing)
    pub fn eval(&self, t: T) -> SVector<T, DIM> {
        //let t = t.into();
        // start with a copy of the original control points array and succesively use it for evaluation
        let mut p: [SVector<T, DIM>; N] = self.control_points;
        // loop up to degree = control_points.len() -1
        for i in 1..=p.len() {
            for j in 0..p.len() - i {
                p[j] = p[j] * (-t + T::from(1.0).unwrap()) + p[j + 1] * t;
            }
        }
        p[0]
    }

    pub fn split(&self, t: T) -> (Self, Self) {
        // start with a copy of the original control points for now
        // TODO how to initialize const generic array without using unsafe?
        let mut left: [SVector<T, DIM>; N] = self.control_points;
        let mut right: [SVector<T, DIM>; N] = self.control_points;
        // these points get overriden each iteration; we save the intermediate results to 'left' and 'right'
        let mut casteljau_points: [SVector<T, DIM>; N] = self.control_points;

        for i in 1..=casteljau_points.len() {
            // save start point of level
            left[i - 1] = casteljau_points[0];
            // save end point of level
            right[right.len() - i] = casteljau_points[right.len() - i];
            // calculate next level of points (one less point each level until we reach one point, the one at t)
            for j in 0..casteljau_points.len() - i {
                casteljau_points[j] = casteljau_points[j] * (-t + T::from(1.0).unwrap())
                    + casteljau_points[j + 1] * t;
            }
        }

        (
            Bezier {
                control_points: left,
            },
            Bezier {
                control_points: right,
            },
        )
    }

    pub fn line_segments(&self, tolerance: T) -> impl Iterator<Item = SVector<T, DIM>> {
        let mut stack = SmallVec::<[Bezier<T, DIM, N>; 32]>::new();
        stack.push(*self);

        from_fn(move || {
            let mut next = None;
            while next == None {
                if let Some(curve) = stack.pop() {
                    if N == 2 || curve.is_linear(tolerance) {
                        next = Some(curve.start());
                    } else {
                        let (left, right) = curve.split(T::from(0.5).unwrap());
                        stack.push(right);
                        stack.push(left);
                    }
                } else {
                    break;
                }
            }
            next
        })
        .chain(core::iter::once(self.end()))
    }

    #[cfg(feature = "alloc")]
    pub fn calculate_offset_segments_buffered<const C: usize>(
        &self,
        tolerance: T,
        _offsets: [T; C],
        _buffers: &mut [Vec<SVector<T, DIM>>; C],
    ) {
        use itertools::*;

        for (_x, _y) in self.line_segments(tolerance).tuple_windows() {}
    }

    pub fn offset_quantized_points<const C: usize>(
        &self,
        tolerance: T,
        offsets: [T; C],
        cross_vec: SVector<T, DIM>,
    ) -> impl Iterator<Item = [SVector<T, DIM>; C]>
    where
        [(); N - 1]: Sized,
    {
        let mut initial = [[T::zero(); DIM].into(); C];
        let mut last = [[T::zero(); DIM].into(); C];
        let d_self = self.derivative();
        let initial_offset_dir = d_self.eval(T::zero()).cross(&cross_vec).normalize();
        let end_offset_dir = d_self.eval(T::one()).cross(&cross_vec).normalize();
        for (pt, offset) in initial.iter_mut().zip(offsets.into_iter()) {
            *pt = self.start() + (initial_offset_dir * offset);
        }
        for (pt, offset) in last.iter_mut().zip(offsets.into_iter()) {
            *pt = self.end() + (end_offset_dir * offset);
        }

        core::iter::once(initial)
            .chain(
                self.line_segments(tolerance)
                    .tuple_windows()
                    .map(move |(start, end)| {
                        let tangent = (end - start).normalize();
                        (start, end, tangent, tangent.cross(&cross_vec).normalize())
                    })
                    .tuple_windows()
                    .map(move |(line1, line2)| {
                        let mut points = [[T::zero(); DIM].into(); C];
                        let offset_dir = (line1.3 + line2.3).normalize();
                        for (pt, offset) in points.iter_mut().zip(offsets.into_iter()) {
                            // let start1 = line1.0 + (line1.3 * offset);
                            // let end2 = line2.1 + (line2.3 * offset);
                            // let intersect = line_intersection(
                            //     start1,
                            //     line1.2,
                            //     T::infinity(),
                            //     end2,
                            //     -line2.2,

                            //     T::infinity(),
                            // );
                            // *pt = if let Some(pt_2) = intersect {
                            //     #[cfg(test)]
                            //     dbg!(&pt_2);
                            //     pt_2
                            // } else {
                            //     #[cfg(test)]
                            //     dbg!(&line1.1);
                            //     line1.1
                            // };

                            // let mid1 = line1.1 + (line1.3 * offset);
                            // let mid2 = line2.0 + (line2.3 * offset);
                            // *pt = (mid1 + mid2) / T::from(2.0).unwrap();

                            *pt = line1.1 + (offset_dir * offset);

                            // #[cfg(test)]
                            // dbg!((*pt - line1.1).magnitude());
                        }
                        points
                    }),
            )
            .chain(core::iter::once(last))
    }

    /*
    let mut self_iter = None;
    let mut left_iter = None;
    let mut right_iter = None;
    let mut right_curve = None;
    let mut left_curve = None;
    let use_self = N == 2 || self.is_linear(tolerance);
    if use_self {
        self_iter = Some([self.start(), self.end()].into_iter());
    } else {
        let (l, r) = self.split(tolerance);
        left_curve = Some(l);
        right_curve = Some(r);
        left_iter = Some(left_curve.unwrap().line_segments(tolerance));
        right_iter = Some(right_curve.unwrap().line_segments(tolerance));
    }
    let mut left_ok = true;

    &mut from_fn(move || {
        if use_self {
            self_iter.as_mut().unwrap().next()
        } else {
            if left_ok {
                let next = left_iter.as_mut().unwrap().next();
                if next.is_some() {
                    return next;
                }
                left_ok = false;
            }
            if !left_ok {
                return right_iter.as_mut().unwrap().next()
            }
            None
        }
    })
    */

    /// Returns the derivative curve of self which has N-1 control points.
    /// The derivative of an nth degree Bézier curve is an (n-1)th degree Bézier curve,
    /// with one fewer term, and new weights w0...wn-1 derived from the
    /// original weights as n(wi+1 - wi). So for a 3rd degree curve, with four weights,
    /// the derivative has three new weights:
    ///     w0 = 3(w1-w0), w'1 = 3(w2-w1) and w'2 = 3(w3-w2).
    pub fn derivative(&self) -> Bezier<T, DIM, { N - 1 }> {
        let _x = self.control_points[0];
        let mut new_points = [[T::zero(); DIM].into(); N - 1];
        for (i, _) in self.control_points.iter().enumerate() {
            new_points[i] =
                (self.control_points[i + 1] - self.control_points[i]) * (T::from(N - 1).unwrap());
            if i == self.control_points.len() - 2 {
                break;
            }
        }
        Bezier::new(new_points)
    }

    // /// Returns the real roots of the Bezier curve along one of its coordinate
    // /// axes (i.e. the control points' axes) or a specific RootFindingError.
    // /// There are the same number of roots as the degree of the curve nroots = degree = N_points-1
    // fn real_roots(&self,
    //     axis: usize,
    //     eps: Option<T>,
    //     max_iter: Option<usize>
    // ) -> Result<ArrayVec<[T; N-1]>, RootFindingError>
    // {
    //     todo!();
    //     // Compute the axis-wise polynomial coefficients e.g. quadratic has N coefs a,b,c in at^2 + bt + c
    //     // to do this generically, we need to find the coefs of the bezier of degree n by binomial expansion
    //     // B_n(t) = sum_1_to_n ( binom(n,i) * s^(n-i) * t^i * p[i])
    //     let mut res: ArrayVec<[T; N-1]> = ArrayVec::new();
    //     let mut npascal:    [T; N] = [ T::from(0.0); N];
    //     let poly_coefs:     [T; N] = [ T::from(0.0); N];

    //     // 1. calculate the n-th row of pascals triangle on a zero-based index (all values for i in the binom(n,i) part)
    //     //    1      N = 0 (wouldn't compile due to index out of bounds)
    //     //  1  (1)   N = 1 (last 1 is always omitted)
    //     // 1  2  (1) N = 2
    //     npascal[0] = 1.0.into();
    //     for i in 1usize..N {
    //         npascal[i] = npascal[i-1] * (N - i + 1) as NativeFloat / i as NativeFloat;
    //     }
    //     // 2. calculate the coefficients to binom(n,i) and p[i] (the s^(n-i) part)

    //     // 3. find candidate points for roots of that curve (compare zero crossings)

    //     // 4. search for roots using newton-raphson algo
    //     let eps = eps.unwrap_or(1e-3.into());
    //     let max_iter = max_iter.unwrap_or(128);

    //     let mut x = T::from(0.0);

    //     let mut iter = 0;
    //     loop {
    //         let f = f(x);
    //         let d = d(x);
    //         if f < eps {
    //             return Ok(x);
    //         }
    //         // if derivative is 0
    //         if d < EPSILON {
    //             // either try to choose a better starting point
    //             if iter == 0 {
    //                 x = x + 1.0;
    //                 iter = iter + 1;
    //                 continue;
    //             // or fail
    //             } else {
    //                 return Err(RootFindingError::ZeroDerivative);
    //             }
    //         }

    //         let x1 = x - f / d;
    //         if (x - x1).abs() < eps {
    //             return Ok(x1);
    //         }

    //         x = x1;
    //         iter = iter + 1;

    //         if iter == max_iter {
    //             return Err(RootFindingError::FailedToConverge);
    //         }
    //     }
    //     res
    // }

    /// Approximates the arc length of the curve by flattening it with straight line segments.
    /// This works quite well, at ~32 segments it should already provide an error in the decimal places
    /// The accuracy gain falls off with more steps so this approximation is unfeasable if desired accuracy is greater than 1-2 decimal places
    pub fn arclen(&self, nsteps: usize) -> T {
        let stepsize = T::from(1.0).unwrap() / (T::from(nsteps).unwrap());
        let mut arclen: T = T::zero();
        for t in 1..nsteps {
            let t = T::from(t).unwrap() / (T::from(nsteps).unwrap());
            let p1 = self.eval(t);
            let p2 = self.eval(t + stepsize);

            arclen = arclen + ComplexField::sqrt((p1 - p2).magnitude_squared());
        }
        return arclen;
    }

    pub fn start(&self) -> SVector<T, DIM> {
        self.control_points[0]
    }

    pub fn end(&self) -> SVector<T, DIM> {
        self.control_points[self.control_points.len() - 1]
    }

    pub fn baseline(&self) -> Bezier<T, DIM, 2> {
        Bezier::new([self.start(), self.end()])
    }

    pub fn is_linear(&self, tolerance: T) -> bool {
        if ({ N } == 2) {
            true
        } else {
            let line = self.baseline();

            self.control_points.iter().all(|point| {
                let d = line.distance_to_point(point);
                d <= tolerance
            })
        }
    }

    /// Return the bounding box of the curve as an array of (min, max) tuples for each dimension (its index)
    pub fn bounding_box(&self) -> (SVector<T, DIM>, SVector<T, DIM>) {
        ([T::zero(); DIM].into(), [T::zero(); DIM].into())
    }
}

impl<T: Scalar, const DIM: usize> Bezier<T, DIM, 2> {
    pub fn distance_squared_to_point(&self, pt: &SVector<T, DIM>) -> T {
        let start = self.start();
        let end = self.end();
        let len_sq = (end - start).magnitude_squared();
        // #[cfg(test)]
        // dbg!(start, l2);
        // if start and endpoint are approx the same, return the distance to either

        let v1 = pt - start;
        let v2 = end - start;
        let param = if len_sq > T::from(EPSILON).unwrap() {
            v1.dot(&v2) / len_sq
        } else {
            T::from(-1.0).unwrap()
        };

        let test_pt = if param < T::zero() {
            start
        } else if param > T::one() {
            end
        } else {
            start + (v2 * param)
        };

        (pt - test_pt).magnitude_squared()
    }

    pub fn distance_to_point(&self, pt: &SVector<T, DIM>) -> T {
        Float::sqrt(self.distance_squared_to_point(pt))
    }

    pub fn intersection(&self, other: Bezier<T, DIM, 2>) -> Option<SVector<T, DIM>> {
        let dx1: SVector<T, DIM> = self.end() - self.start();
        let dx2: SVector<T, DIM> = other.end() - other.start();
        let len1 = dx1.magnitude();
        let len2 = dx2.magnitude();
        let slope1: SVector<T, DIM> = dx1 / len1;
        let slope2: SVector<T, DIM> = dx2 / len2;

        line_intersection(self.start(), slope1, len1, other.start(), slope2, len2)
    }

    fn inside(&self, pt: &SVector<T, DIM>) -> bool {
        pt.into_iter()
            .zip(self.start().into_iter().zip(self.end().into_iter()))
            .all(|(val, (start, end))| (val >= start && val <= end) || (val <= start && val >= end))
    }
}

pub fn line_intersection<T: Scalar, const DIM: usize>(
    start_a: SVector<T, DIM>,
    slope_a: SVector<T, DIM>,
    len_a: T,
    start_b: SVector<T, DIM>,
    slope_b: SVector<T, DIM>,
    len_b: T,
) -> Option<SVector<T, DIM>> {
    let t_axis: SVector<T, DIM> = (start_b - start_a).component_div(&(slope_a - slope_b));
    let mut bad = false;
    let t = t_axis.into_iter().reduce(|acc, x| {
        if x.is_nan() {
            acc
        } else {
            if acc.is_nan() {
                x
            } else {
                if Float::abs(*acc - *x) > T::from(0.05).unwrap() {
                    #[cfg(test)]
                    dbg!(acc, x, start_b, start_a, slope_a, slope_b, t_axis);
                    bad = true
                }
                acc
            }
        }
    });

    if bad || t.is_none() {
        None
    } else {
        let t = *t.unwrap();
        if t > T::from(len_a).unwrap() || t > T::from(len_b).unwrap() || t < T::zero() {
            #[cfg(test)]
            dbg!(t);
            None
        } else {
            Some(start_a + (slope_a * t))
        }
    }
}

#[cfg(test)]
mod tests {
    use nalgebra::{Vector2, Vector3};

    use super::*;

    //use crate::num_traits::{Pow};
    #[test]
    fn eval_endpoints() {
        let points = [
            Vector2::new(0f64, 1.77f64),
            [1.1f64, -1f64].into(),
            [4.3f64, 3f64].into(),
            [3.2f64, -4f64].into(),
            [7.3f64, 2.7f64].into(),
            [8.9f64, 1.7f64].into(),
        ];

        let curve = Bezier::new(points);
        let _b = curve;

        // check if start/end points match
        let start = curve.eval(0.0);
        let err_start = start - points[0];
        assert!(err_start.magnitude_squared() < EPSILON);

        let end = curve.eval(1.0);
        let err_end = end - points[points.len() - 1];
        assert!(err_end.magnitude_squared() < EPSILON);
    }

    #[test]
    fn split_equivalence() {
        // chose some arbitrary control points and construct a cubic bezier
        let bezier = Bezier {
            control_points: [
                [0f64, 1.77f64].into(),
                [2.9f64, 0f64].into(),
                [4.3f64, 3f64].into(),
                [3.2f64, -4f64].into(),
            ],
        };
        // split it at an arbitrary point
        let at = 0.5;
        let (left, right) = bezier.split(at);
        // compare left and right subcurves with parent curve
        // take the difference of the two points which must not exceed the absolute error
        let nsteps: usize = 1000;
        for t in 0..=nsteps {
            let t = t as f64 * 1f64 / (nsteps as f64);
            // check the left part of the split curve
            let mut err = bezier.eval(t / 2.0) - left.eval(t);
            assert!(err.magnitude_squared() < EPSILON);
            // check the right part of the split curve
            err = bezier.eval((t * 0.5) + 0.5) - right.eval(t);
            assert!(err.magnitude_squared() < EPSILON);
        }
    }

    #[test]
    fn test_line_segment() {
        let curve = Bezier::new([Vector2::new(0.0, 0.0), Vector2::new(1.0, 1.0)]);
        let points: Vec<_> = curve.line_segments(0.01).collect();
        assert_eq!(vec![curve.start(), curve.end()], points);
    }

    #[test]
    fn test_quadratic_high_tolerance() {
        let curve = Bezier::new([
            Vector2::new(0.0, 0.0),
            Vector2::new(0.5, 1.0),
            Vector2::new(1.0, 0.0),
        ]);
        let points: Vec<_> = curve.line_segments(1.).collect();
        assert_eq!(vec![curve.start(), curve.end()], points);
    }

    #[test]
    fn test_linear_quadratic() {
        let curve = Bezier::new([
            Vector2::new(0.0, 0.0),
            Vector2::new(0.5, 0.0),
            Vector2::new(1.0, 0.0),
        ]);
        let points: Vec<_> = curve.line_segments(0.01).collect();
        assert_eq!(vec![curve.start(), curve.end()], points);
    }

    #[test]
    fn test_intersection() {
        let curve = Bezier::new([Vector2::new(0.0, 0.0), Vector2::new(1.0, 1.0)]);
        let curve2 = Bezier::new([Vector2::new(0.0, 1.0), Vector2::new(1.0, 0.0)]);
        let intersection = curve.intersection(curve2);
        assert_eq!(intersection, Some(Vector2::new(0.5, 0.5)));
    }

    #[test]
    fn test_3d_segments() {
        let curve = Bezier::new([
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(0.5, 1.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
        ]);
        let points: Vec<_> = curve.line_segments(0.01).collect();
        let offsets = [-0.1, 0.0, 0.1];
        let offset_pts: Vec<_> = curve
            .offset_quantized_points(0.01, offsets, Vector3::new(0., 0., 1.))
            .collect();
        dbg!(offset_pts.len());
        assert_eq!(vec![curve.start(), curve.end()], points);
    }
}
