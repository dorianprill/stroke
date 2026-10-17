//! Scaled differential geometry, independent of the native parameter's units.

use crate::{Point, PointNorm};
use num_traits::{Float, NumCast, One, Zero};

/// Differences are evaluated before degree/scale factors are restored. This
/// keeps directions and curvature usable when the derivative itself overflows.
pub(crate) struct Differential<P: Point> {
    velocity: P,
    acceleration: P,
    degree: P::Scalar,
    curvature_factor: P::Scalar,
    velocity_power: i32,
    acceleration_power: i32,
}

impl<P: Point> Differential<P> {
    pub(crate) fn bezier<const N: usize>(control: [P; N], t: P::Scalar, second: bool) -> Self {
        let zero = P::Scalar::zero();
        let mut work = control;
        let mut result = Self {
            velocity: control[0] * zero,
            acceleration: control[0] * zero,
            degree: zero,
            curvature_factor: zero,
            velocity_power: 0,
            acceleration_power: 0,
        };
        if N < 2 {
            return result;
        }
        result.degree = <P::Scalar as NumCast>::from(N - 1).unwrap();
        result.velocity_power = differences(&mut work, N, 0);
        result.velocity = evaluate(work, N - 1, t);
        if second && N > 2 {
            result.acceleration_power = differences(&mut work, N - 1, result.velocity_power);
            result.acceleration = evaluate(work, N - 2, t);
            result.curvature_factor = <P::Scalar as NumCast>::from(N - 2).unwrap() / result.degree;
        }
        result
    }

    pub(crate) fn derivative(&self) -> P {
        let two = P::Scalar::one() + P::Scalar::one();
        (self.velocity * self.degree) * two.powi(self.velocity_power)
    }
}

impl<P: PointNorm> Differential<P> {
    pub(crate) fn tangent(&self) -> P {
        normalize(self.velocity)
    }

    pub(crate) fn curvature(&self) -> P::Scalar {
        let zero = P::Scalar::zero();
        let vscale = component_scale(self.velocity);
        let ascale = component_scale(self.acceleration);
        if !vscale.is_finite() || !ascale.is_finite() {
            return P::Scalar::nan();
        }
        if vscale == zero || ascale == zero || self.curvature_factor == zero {
            return zero;
        }
        let v = crate::numeric::divide(self.velocity, vscale);
        let a = crate::numeric::divide(self.acceleration, ascale);
        // The norm of the exterior product avoids subtracting nearly equal
        // Gram products. This costs O(DIM^2), with no dimension-specific types.
        let mut wedge = zero;
        for i in 0..P::DIM {
            for j in i + 1..P::DIM {
                wedge = wedge.hypot(difference_of_products(v[i], a[j], v[j], a[i]));
            }
        }
        let speed = v.norm();
        crate::numeric::scaled_product_ratio(
            &[self.curvature_factor, ascale, wedge],
            &[vscale, vscale, speed, speed, speed],
            self.acceleration_power - 2 * self.velocity_power,
        )
    }

    pub(crate) fn normal(&self) -> Option<P> {
        let zero = P::Scalar::zero();
        let vscale = component_scale(self.velocity);
        let ascale = component_scale(self.acceleration);
        if vscale == zero || ascale == zero || !vscale.is_finite() || !ascale.is_finite() {
            return None;
        }
        let direction = crate::numeric::divide(self.velocity, vscale);
        let a = crate::numeric::divide(self.acceleration, ascale);
        // The wedge detects exact collinearity without a squared-norm cutoff.
        if (0..P::DIM).all(|i| {
            (i + 1..P::DIM)
                .all(|j| difference_of_products(direction[i], a[j], direction[j], a[i]) == zero)
        }) {
            return None;
        }
        let v = normalize(direction);
        let mut perpendicular = a - v * v.dot(&a);
        // Reorthogonalize to remove the parallel error left by cancellation.
        perpendicular = perpendicular - v * v.dot(&perpendicular);
        let scale = component_scale(perpendicular);
        if scale == zero || !scale.is_finite() {
            None
        } else {
            let normal = normalize(perpendicular);
            normal.is_finite().then_some(normal)
        }
    }
}

fn component_scale<P: PointNorm>(point: P) -> P::Scalar {
    if !point.is_finite() {
        return P::Scalar::nan();
    }
    (0..P::DIM).fold(P::Scalar::zero(), |scale, i| scale.max(point[i].abs()))
}

fn normalize<P: PointNorm>(point: P) -> P {
    let scale = component_scale(point);
    if scale == P::Scalar::zero() {
        return point;
    }
    let scaled = crate::numeric::divide(point, scale);
    crate::numeric::divide(scaled, scaled.norm())
}

fn difference_of_products<T: Float>(a: T, b: T, c: T, d: T) -> T {
    let product = c * d;
    let error = (-c).mul_add(d, product);
    a.mul_add(b, -product) + error
}

fn differences<P: Point>(work: &mut [P], count: usize, mut power: i32) -> i32 {
    if work[..count]
        .windows(2)
        .any(|pair| !(pair[1] - pair[0]).is_finite())
    {
        let half = <P::Scalar as NumCast>::from(0.5).unwrap();
        for point in &mut work[..count] {
            *point = *point * half;
        }
        power += 1;
    }
    for i in 0..count - 1 {
        work[i] = work[i + 1] - work[i];
    }
    power
}

fn evaluate<P: Point, const N: usize>(mut work: [P; N], count: usize, t: P::Scalar) -> P {
    for level in 1..count {
        for i in 0..count - level {
            work[i] = crate::numeric::lerp(work[i], work[i + 1], t);
        }
    }
    work[0]
}
