use core::slice::*;

use super::point::Point;
use super::*;
use crate::find_root::FindRoot;
use crate::roots::{root_newton_raphson, RootFindingError};

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum BSplineError {
    TooFewControlPoints,
    TooFewKnots,
    TooManyKnots,
    DecreasingKnots,
    KnotDomainViolation,
    EvaluationAlphaIsNaN,
}

#[derive(Debug, Clone, Copy, PartialEq)]
enum KnotVectorKind {
    Unclamped,
    ClampedStart,
    ClampedEnd,
    Clamped,
}

impl Default for KnotVectorKind {
    fn default() -> Self {
        KnotVectorKind::Unclamped
    }
}

impl core::fmt::Display for BSplineError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            BSplineError::TooFewControlPoints => write!(f, "Too few control points"),
            BSplineError::TooFewKnots => write!(f, "Insufficient number of knots"),
            BSplineError::TooManyKnots => write!(f, "Excess number of knots"),
            BSplineError::DecreasingKnots => write!(f, "Knots must be non-decreasing in value"),
            BSplineError::KnotDomainViolation => write!(f, "Parameter outside of knot span"),
            BSplineError::EvaluationAlphaIsNaN => {
                write!(f, "Cox-deBoor evaluation yields alpha NaN value")
            }
        }
    }
}

/// General Implementation of a BSpline with choosable degree, control points and knots.
///
/// Generic parameters:
/// P: const generic points array 'P' as defined by the Point trait
/// F: Any float value used for the knots and interpolation (usually the same as the internal generic parameter within P<F>).
/// const generic parameters:
/// C: Number of control points
/// K: Number of Knots
/// D: Degree of the piecewise function used for interpolation degree = order - 1
/// While C, K, D relate to each other in the following manner
///     K = C + D + 1
#[derive(Clone, Copy)]
pub struct BSpline<P, const K: usize, const C: usize, const D: usize>
where
    P: Point,
{
    /// Knot vector
    knots: [P::Scalar; K],
    /// Knot Vector Kind (Clamped or Unclamped)
    knot_kind: KnotVectorKind,
    /// Control points
    control_points: [P; C],
}

impl<P, const K: usize, const C: usize, const D: usize> Default for BSpline<P, K, C, D>
where
    P: Point,
{
    fn default() -> Self {
        BSpline {
            knots: core::array::from_fn(|_| P::Scalar::default()),
            knot_kind: KnotVectorKind::default(),
            control_points: core::array::from_fn(|_| P::default()),
        }
    }
}

impl<P, const K: usize, const C: usize, const D: usize> BSpline<P, { K }, { C }, { D }>
where
    P: Point,
{
    /// Create a new B-spline curve that interpolates
    /// the `control_points` using a piecewise polynomial of `degree` within intervals specified by the `knots`.
    /// The knots _must_ be sorted in non-decreasing order, the constructor enforces this which may yield undesired results.
    /// The degree is defined as `curve_order - 1`.
    /// Desired curve must have a valid number of control points and knots in relation to its degree or the constructor will return None.
    /// A B-Spline curve requires at least one more control point than the degree (`control_points.len() >
    /// degree`) and the number of knots should be equal to `control_points.len() + degree + 1`.
    /// K = C + D + 1
    pub fn new(
        knots: [P::Scalar; K],
        control_points: [P; C],
    ) -> Result<BSpline<P, { K }, { C }, { D }>, BSplineError> {
        if control_points.len() <= { D } {
            Err(BSplineError::TooFewControlPoints)
        // n_knots = n_control_points + 1 + degree
        } else if knots.len() < control_points.len() + 1 + D {
            Err(BSplineError::TooFewKnots)
        } else if knots.len() > control_points.len() + 1 + D {
            Err(BSplineError::TooManyKnots)
        // if knots array is not non-decreasing
        } else if knots.windows(2).any(|w| w[0] > w[1]) {
            Err(BSplineError::DecreasingKnots)
        } else {
            let is_clamped_front = knots[0..=D]
                .iter()
                .all(|&knot| (knot - knots[0]).abs() < P::Scalar::epsilon());
            let is_clamped_back = knots[knots.len() - D - 1..]
                .iter()
                .all(|&knot| (knot - knots[knots.len() - 1]).abs() < P::Scalar::epsilon());

            let knot_kind = match (is_clamped_front, is_clamped_back) {
                (false, false) => KnotVectorKind::Unclamped,
                (true, false) => KnotVectorKind::ClampedStart,
                (false, true) => KnotVectorKind::ClampedEnd,
                (true, true) => KnotVectorKind::Clamped,
            };

            Ok(BSpline {
                knots,
                knot_kind,
                control_points,
            })
        }
    }

    /// Compute a point on the curve at `t` using iterative de boor algorithm.
    /// The parameter **must** be in the inclusive range of values returned
    /// by `knot_domain()`. If `t` is out of bounds a KnotDomainViolation error is returned.
    pub fn eval(&self, t: P::Scalar) -> Result<P, BSplineError>
    where
        [(); D + 1]: Sized,
    {
        let (knot_start, knot_end) = self.knot_domain();

        // Ensure t is within the knot domain
        if t < knot_start || t > knot_end {
            return Err(BSplineError::KnotDomainViolation);
        }

        // Calculate the point on the B-Spline at `t` using de Boor's algorithm
        let p = self.de_boor_iterative(t)?;
        Ok(p)
    }

    /// Returns an iterator over the control points.
    pub fn control_points(&self) -> Iter<'_, P> {
        self.control_points.iter()
    }

    /// Returns an iterator over the knots.
    pub fn knots(&self) -> Iter<'_, P::Scalar> {
        self.knots.iter()
    }

    /// Returns the knot domain of the B-Spline.
    /// The knot domain is the range of values over which the B-Spline is defined.
    /// The knot domain is defined as the minimum and maximum knot values.
    /// For a clamped B-Spline the first and last knot has multiplicity D+1.
    /// For an unclamped B-Spline the first and last knot has multiplicity 1.
    pub fn knot_domain(&self) -> (P::Scalar, P::Scalar) {
        match self.knot_kind {
            KnotVectorKind::Clamped => self.knot_domain_clamped(),
            KnotVectorKind::ClampedStart => self.knot_domain_clamped_start(),
            KnotVectorKind::ClampedEnd => self.knot_domain_clamped_end(),
            KnotVectorKind::Unclamped => self.knot_domain_unclamped(),
        }
    }

    fn knot_domain_clamped_start(&self) -> (P::Scalar, P::Scalar) {
        // Start from knots[D], end at knots[knots.len() - 1]
        (self.knots[D], self.knots[self.knots.len() - 1])
    }

    fn knot_domain_clamped_end(&self) -> (P::Scalar, P::Scalar) {
        // Start from knots[0], end at knots[knots.len() - 1 - D] (C = K - D - 1)
        (self.knots[0], self.knots[self.knots.len() - 1 - D])
    }

    // Returns the knot domain for a clamped B-Spline
    // where the first and last knot has multiplicity D+1
    fn knot_domain_clamped(&self) -> (P::Scalar, P::Scalar) {
        // The valid domain is from knots[D] to knots[C] (C = K - D - 1)
        (self.knots[D], self.knots[self.knots.len() - 1 - D])
    }

    // Returns the knot domain for an unclamped B-Spline
    // where the first and last knot has multiplicity 1
    // fn knot_domain_unclamped(&self) -> (P::Scalar, P::Scalar) {
    //     (self.knots[0], self.knots[self.knots.len() - 1])
    // }
    fn knot_domain_unclamped(&self) -> (P::Scalar, P::Scalar) {
        // The valid domain is from knots[D] to knots[C] (C = K - D - 1)
        (self.knots[D], self.knots[C])
    }

    // /// Calculates the minimum distance between given 'point' and the curve.
    // /// Uses two passes with the same amount of steps in t:
    // /// 1. coarse search over the whole curve
    // /// 2. fine search around the minimum yielded by the coarse search
    // /// TODO FIXME INVESTIGATE
    // pub fn distance_to_point(&self, point: P) -> P::Scalar
    // where
    //     [(); D + 1]:,
    // {
    //     let nsteps: usize = 64;
    //     let mut tmin: P::Scalar = 0.5.into();
    //     let mut dmin: P::Scalar = (point - self.control_points[0]).squared_length();
    //     let (kstart, kend) = self.knot_domain();
    //     // 1. coarse pass
    //     for i in 0..nsteps {
    //         // calculate next step value
    //         let t: P::Scalar =
    //             kstart + (kend - kstart) * (i as NativeFloat / (nsteps as NativeFloat));
    //         // calculate distance to candidate
    //         let candidate = match self.eval(t) {
    //             Ok(val) => val,
    //             Err(_) => {
    //                 // In case of error, return zero
    //                 // this can never happen as we control the t values passed to eval()
    //                 return P::Scalar::from(0.0);
    //             }
    //         };
    //         if (candidate - point).squared_length() < dmin {
    //             tmin = t;
    //             dmin = (candidate - point).squared_length();
    //         }
    //     }
    //     // 2. fine pass
    //     for i in 0..nsteps {
    //         // calculate next step value ( a 64th of a 64th from first step)
    //         let t: P::Scalar =
    //             kstart + (kend - kstart) * (i as NativeFloat / ((nsteps * nsteps) as NativeFloat));
    //         // calculate distance to candidate centered around tmin from before
    //         let candidate = match self.eval(tmin + t - t * (nsteps as NativeFloat / 2.0)) {
    //             Ok(val) => val,
    //             Err(_) => {
    //                 // In case of error, return zero
    //                 // this can never happen as we control the t values passed to eval()
    //                 return P::Scalar::from(0.0);
    //             }
    //         };
    //         if (candidate - point).squared_length() < dmin {
    //             tmin = t;
    //             dmin = (candidate - point).squared_length();
    //         }
    //     }
    //     dmin.sqrt()
    // }

    /// Iteratively compute de Boor's B-spline algorithm, this computes the recursive
    /// de Boor algorithm tree from the bottom up. At each level we use the results
    /// from the previous one to compute this level and store the results in the
    /// array indices we no longer need to compute the current level (the left one
    /// used computing node j).
    fn de_boor_iterative(&self, t: P::Scalar) -> Result<P, BSplineError>
    where
        [(); D + 1]: Sized,
    {
        // Special case when t equals the last knot value
        if t == self.knots[self.control_points.len()] {
            return Ok(self.control_points[self.control_points.len() - 1]);
        }

        let k = self.knot_span_start_for_t(t);

        let k = match k {
            Some(r) => r,
            None => return Err(BSplineError::KnotDomainViolation),
        };

        // Create a second array for storing intermediate results
        let mut d: [P; D + 1] = [P::default(); D + 1];

        // Copy the control points needed for the first iteration
        // For a B-Spline of degree D, D+1 control points in the knot span contribute to the result
        for j in 0..=D {
            d[j] = self.control_points[j + k - D];
        }

        // // Adjust the start index to prevent negative values
        // let start_index = if k >= D { k - D } else { 0 };
        // // Calculate the number of control points to copy
        // let num_points = if k >= D { D + 1 } else { k + 1 };
        // // Copy the control points needed for the first iteration
        // for j in 0..num_points {
        //     d[j] = self.control_points[start_index + j];
        // }

        // de Boor's algorithm
        // The result overwrites the right of the two point operands used in the current iteration
        // The first iteration computes the result for the first knot span,
        // the second iteration for the next knot span, and so on.
        for r in 1..=D {
            for j in (r..=D).rev() {
                // calculate alpha, avoid division by zero
                let numerator = t - self.knots[j + k - D];
                let denominator = self.knots[j + 1 + k - r] - self.knots[j + k - D];
                let alpha: P::Scalar;
                if numerator < P::Scalar::epsilon() {
                    // we are at the start of a knot span (interpolation d[j-1] -> d[j])
                    alpha = 0.0.into();
                } else if denominator.abs() < P::Scalar::epsilon() {
                    // we are at then end of a knot span (interpolation d[j-1] -> d[j])
                    // this will lead to d[j] simply being copied to d[j]
                    alpha = 1.0.into();
                } else {
                    // inside a knot span (interpolation d[j-1] -> d[j])
                    alpha = numerator / denominator;
                }
                if alpha.is_nan() {
                    return Err(BSplineError::EvaluationAlphaIsNaN);
                }
                // The right of two points is overwritten by the result
                // to be used in the next iteration
                d[j] = d[j - 1] * (1.0 - alpha.into()) + d[j] * alpha;
            }
        }
        Ok(d[D])
    }

    fn knot_span_start_for_t(&self, t: P::Scalar) -> Option<usize> {
        match self.knot_kind {
            KnotVectorKind::Unclamped => self.knot_span_start_for_t_unclamped(t),
            KnotVectorKind::ClampedStart => self.knot_span_start_for_t_clamped_start(t),
            KnotVectorKind::Clamped => self.knot_span_start_for_t_clamped(t),
            KnotVectorKind::ClampedEnd => self.knot_span_start_for_t_clamped_end(t),
        }
    }

    /// Return the index of the start of the knot span for given `t` for an unclamped B-Spline.
    /// This is the index of the first element that satisfies:
    ///     knots[i] <= t < knots[i+1]
    /// Because the knot vector is non-decreasing, this function uses binary search.
    /// If no element greater than the value passed is found, the function returns None.
    fn knot_span_start_for_t_unclamped(&self, t: P::Scalar) -> Option<usize> {
        if t == self.knots[C] {
            return Some(C - 1);
        }

        let mut low = D;
        let mut high = C - 1;
        while low <= high {
            let mid = (low + high) / 2;
            if t >= self.knots[mid] && t < self.knots[mid + 1] {
                return Some(mid);
            }
            if t < self.knots[mid] {
                if mid == 0 {
                    return None;
                }
                high = mid - 1;
            } else {
                low = mid + 1;
            }
        }
        None
    }

    /// Return the index of the start of the knot span for given `t` for a clamped B-Spline.
    /// This is the index of the first element that satisfies:
    ///     knots[i] <= t < knots[i+1] where i >= degree and i < knots.len() - degree - 1
    /// Because the knot vector is non-decreasing, this function uses binary search.
    /// If no element greater than the value passed is found, the function returns None.
    fn knot_span_start_for_t_clamped(&self, t: P::Scalar) -> Option<usize> {
        let n = self.control_points.len() - 1; // n = number of control points - 1
        let p = D; // degree
        if t == self.knots[n + 1] {
            return Some(n);
        }
        let mut low = p;
        let mut high = n;
        while low <= high {
            let mid = (low + high) / 2;
            if t >= self.knots[mid] && t < self.knots[mid + 1] {
                return Some(mid);
            } else if t < self.knots[mid] {
                high = mid - 1;
            } else {
                low = mid + 1;
            }
        }
        None
    }

    fn knot_span_start_for_t_clamped_start(&self, t: P::Scalar) -> Option<usize> {
        let n = self.control_points.len() - 1; // n = number of control points - 1
        let p = D; // Degree of the spline
        let m = self.knots.len() - 1; // Last index of the knot vector

        // Special case when t equals the last knot value
        if t == self.knots[m] {
            return Some(n);
        }

        // Valid parameter domain is from knots[p] to knots[m]
        if t < self.knots[p] || t > self.knots[m] {
            return None;
        }

        // Initialize binary search indices
        let mut low = p;
        let mut high = m - 1; // Since we compare t with knots[mid + 1], high should be m - 1

        // Binary search to find the knot span
        while low <= high {
            let mid = (low + high) / 2;
            if t >= self.knots[mid] && t < self.knots[mid + 1] {
                return Some(mid);
            } else if t < self.knots[mid] {
                high = mid - 1;
            } else {
                low = mid + 1;
            }
        }

        None // t is not within any knot span
    }

    fn knot_span_start_for_t_clamped_end(&self, t: P::Scalar) -> Option<usize> {
        let n = self.control_points.len() - 1; // n = number of control points - 1

        // Special case when t equals knots[n + 1]
        if t == self.knots[n + 1] {
            return Some(n);
        }

        // Valid parameter domain is from knots[0] to knots[n + 1]
        if t < self.knots[0] || t > self.knots[n + 1] {
            return None;
        }

        // Initialize binary search indices
        let mut low = 0;
        let mut high = n; // Since we compare t with knots[mid + 1], high should be n

        // Binary search to find the knot span
        while low <= high {
            let mid = (low + high) / 2;
            if t >= self.knots[mid] && t < self.knots[mid + 1] {
                return Some(mid);
            } else if t < self.knots[mid] {
                high = mid - 1;
            } else {
                low = mid + 1;
            }
        }

        None // t is not within any knot span
    }

    /// Approximates the arc length of the curve by flattening it with straight line segments.
    /// This approximation is unfeasable if desired accuracy is greater than ~2 decimal places
    pub fn arclen(&self, nsteps: usize) -> P::Scalar
    where
        [(); D + 1]: Sized,
    {
        let stepsize = P::Scalar::from(1.0 / (nsteps as NativeFloat));
        let mut arclen: P::Scalar = 0.0.into();
        // evaluate the curve, t needs to be inside the knot domain!
        // we need to map [0...1] to kmin..kmax
        let (kmin, kmax) = self.knot_domain();
        for t in 0..=nsteps {
            let t = kmin
                + (P::Scalar::from(t as NativeFloat) / P::Scalar::from(nsteps as NativeFloat))
                    * (kmax - kmin);
            let p1 = self.eval(t).unwrap_or_else(|_| {
                // should never happen as we control the t values passed to eval()
                P::default()
            });
            let p2 = self.eval(t + stepsize).unwrap_or_else(|_| {
                // should never happen as we control the t values passed to eval()
                P::default()
            });
            arclen = arclen + (p1 - p2).squared_length().sqrt();
        }
        arclen
    }

    /// Returns the derivative curve of self which has N-1 control points.
    /// The derivative of an nth degree B-Spline curve is an (n-1)th degree (d) B-Spline curve,
    /// with the same knot vector, and new control points Q0...Qn-1 derived from the
    /// original control points Pi as:
    ///                 d
    /// Qi =    ----------------- (P[i+1]-P[i])
    ///         k[i+d+1] - k[i+1].
    /// with degree = curve_order - 1
    /// TODO test & verify!
    pub fn derivative(&self) -> BSpline<P, K, { C - 1 }, { D - 1 }> {
        let mut new_points: [P; C - 1] = [P::default(); C - 1];
        for (i, _) in self.control_points.iter().enumerate() {
            new_points[i] = (self.control_points[i + 1] - self.control_points[i])
                * ((D) as NativeFloat / (self.knots[i + D + 1] - self.knots[i + 1]).into());
            if i == self.control_points.len() - 2 {
                break;
            }
        }
        BSpline {
            knots: self.knots,
            knot_kind: self.knot_kind,
            control_points: new_points,
        }
    }

    fn derivative_curve(&self) -> Result<BSpline<P, { K - 2 }, { C - 1 }, { D - 1 }>, RootFindingError>
    where
        [(); D - 1]: Sized,
        [(); C - 1]: Sized,
        [(); K - 2]: Sized,
    {
        let derivative_knots: [P::Scalar; K - 2] =
            core::array::from_fn(|i| self.knots[i + 1]);
        let derivative_points: [P; C - 1] = core::array::from_fn(|i| {
            (self.control_points[i + 1] - self.control_points[i])
                * ((D) as NativeFloat / (self.knots[i + D + 1] - self.knots[i + 1]).into())
        });
        BSpline::new(derivative_knots, derivative_points)
            .map_err(|_| RootFindingError::FailedToConverge)
    }

    /// Find a root for a particular axis using Newton-Raphson on the scalar axis function.
    pub fn root_newton_axis(
        &self,
        value: P::Scalar,
        axis: usize,
        start: P::Scalar,
        eps: Option<P::Scalar>,
        max_iter: Option<usize>,
    ) -> Result<P::Scalar, RootFindingError>
    where
        [(); D + 1]: Sized,
        [(); D - 1]: Sized,
        [(); D]: Sized,
        [(); C - 1]: Sized,
        [(); (D - 1) + 1]: Sized,
        [(); K - 2]: Sized,
    {
        FindRoot::root_newton_axis(self, value, axis, start, eps, max_iter)
    }
}

impl<P, const K: usize, const C: usize, const D: usize> FindRoot<P> for BSpline<P, K, C, D>
where
    P: Point,
    [(); D + 1]: Sized,
    [(); D - 1]: Sized,
    [(); D]: Sized,
    [(); C - 1]: Sized,
    [(); (D - 1) + 1]: Sized,
    [(); K - 2]: Sized,
{
    fn parameter_domain(&self) -> (P::Scalar, P::Scalar) {
        self.knot_domain()
    }

    fn axis_value(&self, t: P::Scalar, axis: usize) -> Result<P::Scalar, RootFindingError> {
        self.eval(t)
            .map(|p| p.axis(axis))
            .map_err(|_| RootFindingError::FailedToConverge)
    }

    fn axis_derivative(&self, t: P::Scalar, axis: usize) -> Result<P::Scalar, RootFindingError> {
        let derivative = self.derivative_curve()?;
        derivative
            .eval(t)
            .map(|p| p.axis(axis))
            .map_err(|_| RootFindingError::FailedToConverge)
    }

    fn root_newton_axis(
        &self,
        value: P::Scalar,
        axis: usize,
        start: P::Scalar,
        eps: Option<P::Scalar>,
        max_iter: Option<usize>,
    ) -> Result<P::Scalar, RootFindingError> {
        let eps = eps.unwrap_or_else(|| P::Scalar::from(1e-6 as NativeFloat));
        let max_iter = max_iter.unwrap_or(64);
        let (kmin, kmax) = self.knot_domain();
        let derivative = self.derivative_curve()?;

        let clamp_t = |mut t: P::Scalar| {
            if t < kmin {
                t = kmin;
            } else if t > kmax {
                t = kmax;
            }
            t
        };

        let start = clamp_t(start);
        let fx = |x: P::Scalar| {
            let t = clamp_t(x);
            self.axis_value(t, axis).map(|v| v - value)
        };
        let dx = |x: P::Scalar| {
            let t = clamp_t(x);
            derivative
                .eval(t)
                .map(|p| p.axis(axis))
                .map_err(|_| RootFindingError::FailedToConverge)
        };

        let root = root_newton_raphson(start, fx, dx, eps, max_iter)?;
        Ok(clamp_t(root))
    }
}

#[cfg(test)]
mod tests {
    //use std;
    use super::PointN;
    use super::*;
    //use crate::num_traits::{Pow};

    #[test]
    fn degree_1_clamped_construct_and_eval_endpoints() {
        // Define the control points for a degree 1 B-Spline (essentially a line)
        let points = [
            PointN::new([0f64, 0f64]),   // Start of the line at origin
            PointN::new([10f64, 10f64]), // End of the line at (10, 10)
        ];

        // Define a uniform clamped knot vector for a degree 1 B-spline
        let knots: [f64; 4] = [0., 0., 1., 1.];

        // Create the B-Spline
        let curve: BSpline<PointN<f64, 2>, 4, 2, 1> = BSpline::new(knots, points).unwrap();

        assert!(curve.knot_kind == KnotVectorKind::Clamped);

        // evaluate along the curve to find out of bounds errors
        let (kmin, kmax) = curve.knot_domain();
        let nsteps: usize = 100;
        for t in 0..=nsteps {
            let t = kmin + (t as f64 / nsteps as f64) * (kmax - kmin);
            //dbg!(kmin, t, kmax);
            // this curve shouldn't fail to evaluate numerically unless knot index is out of bounds
            let _ = curve.eval(t).unwrap_or_else(|e| {
                panic!("Unexpected error in eval: {:?}", e);
            });
        }

        // Evaluate the B-Spline at the start (t=0) and end (t=1) of the knot vector
        let start = curve.eval(0.).unwrap();
        let end = curve.eval(1.).unwrap();

        // Check that the start and end points match the control points
        assert_eq!(start, points[0]);
        assert_eq!(end, points[1]);
    }

    #[test]
    fn degree_1_unclamped_construct_and_eval() {
        // Define the control points for a degree 1 B-Spline (essentially a line)
        let points = [
            PointN::new([0f64, 0f64]),   // Start of the line at origin
            PointN::new([10f64, 10f64]), // End of the line at (10, 10)
        ];

        // Define a uniform unclamped knot vector for a degree 1 B-spline
        let knots: [f64; 4] = [0., 1., 2., 3.];

        // Create the B-Spline
        let curve: BSpline<PointN<f64, 2>, 4, 2, 1> = BSpline::new(knots, points).unwrap();

        assert!(curve.knot_kind == KnotVectorKind::Unclamped);

        let (kmin, kmax) = curve.knot_domain();

        // evaluate along the curve to find out of bounds errors
        let nsteps: usize = 100;
        for t in 0..=nsteps {
            let t = kmin + (t as f64 / nsteps as f64) * (kmax - kmin);
            //dbg!(kmin, t, kmax);
            // this curve shouldn't fail to evaluate numerically unless knot index is out of bounds
            let _ = curve.eval(t).unwrap_or_else(|e| {
                panic!("Unexpected error in eval: {:?}", e);
            });
        }

        // Evaluate the B-Spline at the start (t=0) and end (t=1) of the knot vector
        let mid = match curve.eval(1.5) {
            Ok(p) => p,
            Err(e) => panic!("Error evaluating B-Spline: {}", e),
        };

        // Check that the start and end points have equal distance to the midpoint of the line
        assert!(
            (points[1] - mid).squared_length().sqrt() - (points[0] - mid).squared_length().sqrt()
                < EPSILON
        );
    }

    #[test]
    fn degree_2_clamped_construct_and_eval_endpoints() {
        // Define the control points for a degree 2 B-Spline
        let points = [
            PointN::new([0f64, 0f64]),   // Start of the line at origin
            PointN::new([10f64, 10f64]), // End of the line at (10, 10)
            PointN::new([20f64, 0f64]),  // End of the line at (20, 0)
        ];

        // Define a uniform clamped knot vector for a degree 2 B-spline
        let knots: [f64; 6] = [0., 0., 0., 1., 1., 1.];

        // Create the B-Spline
        let curve: BSpline<PointN<f64, 2>, 6, 3, 2> = BSpline::new(knots, points).unwrap();

        assert!(curve.knot_kind == KnotVectorKind::Clamped);

        let (kmin, kmax) = curve.knot_domain();

        // evaluate along the curve to find out of bounds errors
        let nsteps: usize = 100;
        for t in 0..=nsteps {
            let t = kmin + (t as f64 / nsteps as f64) * (kmax - kmin);
            //dbg!(kmin, t, kmax);
            // this curve shouldn't fail to evaluate numerically unless knot index is out of bounds
            let _ = curve.eval(t).unwrap_or_else(|e| {
                panic!("Unexpected error in eval: {:?}", e);
            });
        }

        // Evaluate the B-Spline at the start (t=0) and end (t=1) of the knot vector
        let start = curve.eval(0.).unwrap();
        let end = curve.eval(1.).unwrap();

        // Check that the start and end points match the control points
        assert_eq!(start, points[0]);
        assert_eq!(end, points[2]);
    }

    #[test]
    fn degree_2_unclamped_construct_and_eval() {
        // Define the control points for a degree 2 B-Spline
        let points = [
            PointN::new([0f64, 0f64]),   // Start of the line at origin
            PointN::new([10f64, 10f64]), // End of the line at (10, 10)
            PointN::new([20f64, 0f64]),  // End of the line at (20, 0)
        ];

        // Define a uniform unclamped knot vector for a degree 2 B-spline
        let knots: [f64; 6] = [0., 1., 2., 3., 4., 5.];

        // Create the B-Spline
        let curve: BSpline<PointN<f64, 2>, 6, 3, 2> = BSpline::new(knots, points).unwrap();

        assert!(curve.knot_kind == KnotVectorKind::Unclamped);

        let (kmin, kmax) = curve.knot_domain();
        // evaluate along the curve to find out of bounds errors
        let nsteps: usize = 100;
        for t in 0..=nsteps {
            let t = kmin + (t as f64 / nsteps as f64) * (kmax - kmin);
            //dbg!(kmin, t, kmax);
            // this curve shouldn't fail to evaluate numerically unless knot index is out of bounds
            let _ = curve.eval(t).unwrap_or_else(|e| {
                panic!("Unexpected error in eval: {:?}", e);
            });
        }

        // Evaluate the B-Spline at the start (t=0) mid (t=2.5) and end (t=1) of the knot vector
        let mid = match curve.eval(2.5) {
            Ok(p) => p,
            Err(e) => panic!("Error evaluating B-Spline: {}", e),
        };

        // Check that the start and end points have equal distance to the midpoint of the line
        assert!(
            (points[1] - mid).squared_length().sqrt() - (points[0] - mid).squared_length().sqrt()
                < EPSILON
        );
    }

    #[test]
    fn degree_3_clamped_construct_and_eval_endpoints() {
        // degree 3, 4 control points => 4+3+1=8 knots
        let points = [
            PointN::new([0f64, 1.77f64]),
            PointN::new([1.1f64, -1f64]),
            PointN::new([4.3f64, 3f64]),
            PointN::new([3.2f64, -4f64]),
            PointN::new([2.3f64, 1.77f64]),
        ];
        // clamped b-spline to it goes through the first and last control point
        let knots: [f64; 9] = [0., 0., 0., 0., 1., 2., 2., 2., 2.];

        // try to make a b-spline with the given parameters
        let b: Result<BSpline<PointN<f64, 2>, 9, 5, 3>, BSplineError> = BSpline::new(knots, points);

        let curve = match b {
            Err(e) => panic!("Error: {:?}", e),
            Ok(b) => b,
        };

        assert!(curve.knot_kind == KnotVectorKind::Clamped);
        // evaluate the curve, t needs to be inside the knot domain!
        // we need to map [0...1] to kmin..kmax
        let (kmin, kmax) = curve.knot_domain();
        let nsteps: usize = 100;
        for t in 0..=nsteps {
            let t = kmin + (t as f64 / nsteps as f64) * (kmax - kmin);
            //dbg!(kmin, t, kmax);
            // this curve shouldn't fail to evaluate numerically unless knot index is out of bounds
            let _ = curve.eval(t).unwrap_or_else(|e| {
                panic!("Unexpected error in eval: {:?}", e);
            });
        }

        // check if eval returns the correct error if t is out of bounds
        assert_eq!(
            curve.eval(-0.1).err(),
            Some(BSplineError::KnotDomainViolation)
        );
        assert_eq!(
            curve.eval(2.1).err(),
            Some(BSplineError::KnotDomainViolation)
        );
        // check if eval works correctly for endpoints
        assert!(curve.eval(kmin).is_ok());
        assert!(curve.eval(kmax).is_ok());
        assert_eq!(curve.eval(kmin).unwrap(), points[0]);
        assert_eq!(curve.eval(kmax).unwrap(), points[4]);
    }

    #[test]
    fn degree_3_unclamped_construct_and_eval() {
        // degree 3, 4 control points => 4+3+1=8 knots
        let points = [
            PointN::new([0f64, 1.77f64]),
            PointN::new([1.1f64, -1f64]),
            PointN::new([4.3f64, 3f64]),
            PointN::new([3.2f64, -4f64]),
            PointN::new([2.3f64, 1.77f64]),
        ];
        // unclamped b-spline (not necessarily going through the first and last control point)
        let knots: [f64; 9] = [0., 1., 2., 3., 4., 5., 6., 7., 8.];

        // try to make a b-spline with the given parameters
        let b: Result<BSpline<PointN<f64, 2>, 9, 5, 3>, BSplineError> = BSpline::new(knots, points);

        let curve = match b {
            Err(e) => panic!("Error: {:?}", e),
            Ok(b) => b,
        };

        assert!(curve.knot_kind == KnotVectorKind::Unclamped);
        // evaluate the curve, t needs to be inside the knot domain!
        // we need to map [0...1] to kmin..kmax
        let (kmin, kmax) = curve.knot_domain();
        assert!(kmin == 3.0 && kmax == 5.0);
        let nsteps: usize = 100;
        for t in 0..=nsteps {
            let t = kmin + (t as f64 / nsteps as f64) * (kmax - kmin);
            //dbg!(kmin, t, kmax);
            // this curve shouldn't fail to evaluate numerically unless knot index is out of bounds
            let _ = curve.eval(t).unwrap_or_else(|e| {
                panic!("Unexpected error in eval: {:?}", e);
            });
        }

        // check if eval returns the correct error if t is out of bounds
        assert_eq!(
            curve.eval(-0.1).err(),
            Some(BSplineError::KnotDomainViolation)
        );
        assert_eq!(
            curve.eval(8.1).err(),
            Some(BSplineError::KnotDomainViolation)
        );
        // check if eval works correctly for endpoints
        // since the curve is unclampend, the curves endpoints are not necessarily on the control points
        assert!(curve.eval(kmin).is_ok());
        assert!(curve.eval(kmax).is_ok());
    }

    #[test]
    fn degree_4_clamped_construct_and_eval() {
        // degree 4, 5 control points => 5+4+1=10 knots
        let points = [
            PointN::new([0f64, 1.77f64]),
            PointN::new([1.1f64, -1f64]),
            PointN::new([4.3f64, 3f64]),
            PointN::new([3.2f64, -4f64]),
            PointN::new([2.3f64, 1.77f64]),
        ];
        // clamped b-spline (goes through the first and last control point)
        let knots: [f64; 10] = [0., 0., 0., 0., 0., 3., 3., 3., 3., 3.];

        // try to make a b-spline with the given parameters
        let b: Result<BSpline<PointN<f64, 2>, 10, 5, 4>, BSplineError> =
            BSpline::new(knots, points);

        let curve = match b {
            Err(e) => panic!("Error: {:?}", e),
            Ok(b) => b,
        };

        assert!(curve.knot_kind == KnotVectorKind::Clamped);
        // evaluate the curve, t needs to be inside the knot domain!
        // we need to map [0...1] to kmin..kmax
        let (kmin, kmax) = curve.knot_domain();
        assert!(kmin == 0.0 && kmax == 3.0);
        let nsteps: usize = 100;
        for t in 0..=nsteps {
            let t = kmin + (t as f64 / nsteps as f64) * (kmax - kmin);
            //dbg!(kmin, t, kmax);
            // this curve shouldn't fail to evaluate numerically unless knot index is out of bounds
            let _ = curve.eval(t).unwrap_or_else(|e| {
                panic!("Unexpected error in eval: {:?}", e);
            });
        }

        // check if eval returns the correct error if t is out of bounds
        assert_eq!(
            curve.eval(-0.1).err(),
            Some(BSplineError::KnotDomainViolation)
        );
        assert_eq!(
            curve.eval(3.1).err(),
            Some(BSplineError::KnotDomainViolation)
        );
        // check if eval works correctly for endpoints
        assert!(curve.eval(kmin).is_ok());
        assert!(curve.eval(kmax).is_ok());
        assert_eq!(curve.eval(kmin).unwrap(), points[0]);
        assert_eq!(curve.eval(kmax).unwrap(), points[4]);
    }

    #[test]
    fn degree_4_unclamped_construct_and_eval() {
        // degree 4, 5 control points => 5+4+1=10 knots
        let points = [
            PointN::new([0f64, 1.77f64]),
            PointN::new([1.1f64, -1f64]),
            PointN::new([4.3f64, 3f64]),
            PointN::new([3.2f64, -4f64]),
            PointN::new([2.3f64, 1.77f64]),
            PointN::new([0.0, 0.0]),
        ];
        // unclamped b-spline (not necessarily going through the first and last control point)
        let knots: [f64; 11] = [0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10.];

        // try to make a b-spline with the given parameters
        let b: Result<BSpline<PointN<f64, 2>, 11, 6, 4>, BSplineError> =
            BSpline::new(knots, points);

        let curve = match b {
            Err(e) => panic!("Error: {:?}", e),
            Ok(b) => b,
        };

        assert!(curve.knot_kind == KnotVectorKind::Unclamped);
        // evaluate the curve, t needs to be inside the knot domain!
        // we need to map [0...1] to kmin..kmax
        let (kmin, kmax) = curve.knot_domain();
        assert!(kmin == 4.0 && kmax == 6.0);
        let nsteps: usize = 100;
        for t in 0..=nsteps {
            let t = kmin + (t as f64 / nsteps as f64) * (kmax - kmin);
            //dbg!(kmin, t, kmax);
            // this curve shouldn't fail to evaluate numerically unless knot index is out of bounds
            let _ = curve.eval(t).unwrap_or_else(|e| {
                panic!("Unexpected error in eval: {:?}", e);
            });
        }

        // check if eval returns the correct error if t is out of bounds
        assert_eq!(
            curve.eval(-0.1).err(),
            Some(BSplineError::KnotDomainViolation)
        );
        assert_eq!(
            curve.eval(9.1).err(),
            Some(BSplineError::KnotDomainViolation)
        );
    }

    #[test]
    fn root_newton_axis_linear() {
        let points = [PointN::new([0f64]), PointN::new([2f64])];
        let knots: [f64; 4] = [0.0, 0.0, 1.0, 1.0];
        let curve: BSpline<PointN<f64, 1>, 4, 2, 1> = BSpline::new(knots, points).unwrap();

        let root = curve
            .root_newton_axis(1.0, 0, 0.25, None, Some(64))
            .unwrap();
        assert!((root - 0.5).abs() < 1e-6);
    }

    #[test]
    fn root_newton_axis_zero_derivative() {
        let points = [PointN::new([1f64]), PointN::new([1f64])];
        let knots: [f64; 4] = [0.0, 0.0, 1.0, 1.0];
        let curve: BSpline<PointN<f64, 1>, 4, 2, 1> = BSpline::new(knots, points).unwrap();

        let result = curve.root_newton_axis(0.0, 0, 0.5, None, Some(8));
        assert!(matches!(result, Err(RootFindingError::ZeroDerivative)));
    }
}
