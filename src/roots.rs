//! Polynomial root helpers using `tinyvec::ArrayVec`.
//!
//! Available functions:
//! - `root_newton_raphson()`
use num_traits::Float;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RootFindingError {
    NoRootsFound,
    FailedToConverge,
    MaxIterationsReached,
    ZeroDerivative,
}

/// Find a single (any) root of the function f(x) = 0 close to a given a start value, f and its derivative f'.
/// This function cannot predict which root is going to be found.
/// Uses the Newton-Raphson method because it is suited for splines as they are cont. differentiable.
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
