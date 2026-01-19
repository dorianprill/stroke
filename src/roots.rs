/// roots.rs
/// const-generic functions to find polynomial roots using tinyvec::ArrayVec
/// available functions:
///   roots_square()
///   roots_cubic()
///   roots_newton() WIP
/// there will probably be no function for quartics, since specialized beziers only go up to cubic
use super::*;

#[allow(dead_code)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RootFindingError {
    NoRootsFound,
    FailedToConverge,
    MaxIterationsReached,
    ZeroDerivative,
}

/// Solve for the roots of the polynomial at^2 + bt + c
/// Returns an ArrayVec of roots in the order
/// needs to be called for x and y components separately
#[allow(dead_code)]
pub(crate) fn roots_square(
    a: NativeFloat,
    b: NativeFloat,
    c: NativeFloat,
) -> ArrayVec<[NativeFloat; 2]> {
    let mut result = ArrayVec::new();

    // check if can be handled below quadratic order
    if a.abs() < EPSILON {
        if b.abs() < EPSILON {
            // no solutions
            return result;
        }
        // is linear equation
        result.push(-c / b);
        return result;
    }
    // is quadratic equation
    let delta = b * b - a * c * 4.0;
    if delta > 0.0 {
        let sqrt_delta = delta.sqrt();
        result.push((-b - sqrt_delta) / (a * 2.0));
        result.push((-b + sqrt_delta) / (a * 2.0));
    } else if delta.abs() < EPSILON {
        result.push(-b / (a * 2.0));
    }
    result
}

/// Compute the real roots of the cubic bezier function with
/// parameters of the form a*t^3 + b*t^2 + c*t + d for each dimension
/// using cardano's algorithm (code adapted from github.com/nical/lyon)
/// returns an ArrayVec of the present roots (max 3)
#[allow(clippy::many_single_char_names)] // this is math, get over it
#[allow(dead_code)]
pub(crate) fn roots_cubic(
    a: NativeFloat,
    b: NativeFloat,
    c: NativeFloat,
    d: NativeFloat,
) -> ArrayVec<[NativeFloat; 3]> {
    let mut result = ArrayVec::new();
    let pi: NativeFloat = core::f32::consts::PI.into();

    // check if can be handled below cubic order
    if a.abs() < EPSILON {
        if b.abs() < EPSILON {
            if c.abs() < EPSILON {
                // no solutions
                return result;
            }
            // is linear equation
            result.push(-d / c);
            return result;
        }
        // is quadratic equation
        let delta = c * c - b * d * 4.0;
        if delta > 0.0 {
            let sqrt_delta = delta.sqrt();
            result.push((-c - sqrt_delta) / (b * 2.0));
            result.push((-c + sqrt_delta) / (b * 2.0));
        } else if delta.abs() < EPSILON {
            result.push(-c / (b * 2.0));
        }
        return result;
    }

    // is cubic equation -> use cardano's algorithm
    let frac_1_3 = NativeFloat::from(1.0 / 3.0);

    let bn = b / a;
    let cn = c / a;
    let dn = d / a;

    let delta0: NativeFloat = (cn * 3.0 - bn * bn) / 9.0;
    let delta1: NativeFloat = (bn * cn * 9.0 - dn * 27.0 - bn * bn * bn * 2.0) / 54.0;
    let delta_01: NativeFloat = delta0 * delta0 * delta0 + delta1 * delta1;

    if delta_01 >= NativeFloat::from(0.0) {
        let delta_p_sqrt: NativeFloat = delta1 + delta_01.sqrt();
        let delta_m_sqrt: NativeFloat = delta1 - delta_01.sqrt();

        let s = delta_p_sqrt.signum() * delta_p_sqrt.abs().powf(frac_1_3);
        let t = delta_m_sqrt.signum() * delta_m_sqrt.abs().powf(frac_1_3);

        result.push(-bn * frac_1_3 + (s + t));

        // Don't add the repeated root when s + t == 0.
        if (s - t).abs() < EPSILON && (s + t).abs() >= EPSILON {
            result.push(-bn * frac_1_3 - (s + t) / 2.0);
        }
    } else {
        let theta = (delta1 / (-delta0 * delta0 * delta0).sqrt()).acos();
        let two_sqrt_delta0 = (-delta0).sqrt() * 2.0;
        result.push(two_sqrt_delta0 * Float::cos(theta * frac_1_3) - bn * frac_1_3);
        result.push(two_sqrt_delta0 * Float::cos((theta + 2.0 * pi) * frac_1_3) - bn * frac_1_3);
        result.push(two_sqrt_delta0 * Float::cos((theta + 4.0 * pi) * frac_1_3) - bn * frac_1_3);
    }

    result
}

// /// Find a single (any) root of the function f(x) = 0 close to a given a start value, f and its derivative f'.
// /// This function cannot predict which root is going to be found.
// /// Uses the Newton-Raphson method because it is suited for splines as they are cont. differentiable.
// Parameters:
//   start:    Starting point for the search
//   f:        The function for which to find the root
//   d:        d=f', the derivative of f
//   eps:      Desired accuracy of solution
//   max_iter: Number of iterations until a solution shall be found
pub(crate) fn root_newton_raphson<F, Func, Deriv>(
    start: F,
    f: Func,
    d: Deriv,
    eps: F,
    max_iter: usize,
) -> Result<F, RootFindingError>
where
    F: Float,
    Func: Fn(F) -> Result<F, RootFindingError>,
    Deriv: Fn(F) -> Result<F, RootFindingError>,
{
    let mut x = start;
    for _ in 0..max_iter {
        let fx = f(x)?;
        if fx.abs() <= eps {
            return Ok(x);
        }
        let dx = d(x)?;
        if dx.abs() <= eps {
            return Err(RootFindingError::ZeroDerivative);
        }
        let x1 = x - fx / dx;
        if (x1 - x).abs() <= eps {
            return Ok(x1);
        }
        if x1.is_nan() {
            return Err(RootFindingError::FailedToConverge);
        }
        x = x1;
    }
    Err(RootFindingError::MaxIterationsReached)
}
