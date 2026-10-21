//! Newton's method with parameter bounds and coordinate-error checks.
use num_traits::Float;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RootFindingError {
    /// The axis index is outside the point dimension.
    InvalidAxis,
    /// A target, initial guess, or parameter domain is invalid or non-finite.
    InvalidParameter,
    /// The absolute coordinate residual tolerance is not finite and positive.
    InvalidAccuracy,
    /// At least one Newton update must be allowed.
    InvalidIterationLimit,
    /// Function, derivative, or Newton-step arithmetic is non-finite.
    NonFiniteComputation,
    /// The solver returned no isolated roots.
    NoRootsFound,
    /// Parameter rounding or clamping prevented further progress.
    FailedToConverge,
    /// The search used all permitted Newton updates.
    MaxIterationsReached,
    /// A Newton update encountered a zero derivative.
    ZeroDerivative,
}

impl core::fmt::Display for RootFindingError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.write_str(match self {
            Self::InvalidAxis => "root-search axis is out of bounds",
            Self::InvalidParameter => "invalid root-search parameter",
            Self::InvalidAccuracy => "root-search accuracy must be finite and positive",
            Self::InvalidIterationLimit => "root-search iteration limit must be positive",
            Self::NonFiniteComputation => "root-search arithmetic produced a non-finite value",
            Self::NoRootsFound => "no roots found",
            Self::FailedToConverge => "root search failed to converge",
            Self::MaxIterationsReached => "root search reached its iteration limit",
            Self::ZeroDerivative => "root search encountered a zero derivative",
        })
    }
}

impl core::error::Error for RootFindingError {}

/// Search within the parameter domain using Newton's method.
/// Success requires a finite function value within `eps` of zero at the returned
/// parameter. A small parameter step alone does not satisfy the tolerance.
pub(crate) fn root_newton_raphson<F, Func, Deriv>(
    start: F,
    domain: (F, F),
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
    if !eps.is_finite() || eps <= F::zero() {
        return Err(RootFindingError::InvalidAccuracy);
    }
    if max_iter == 0 {
        return Err(RootFindingError::InvalidIterationLimit);
    }
    let (min, max) = domain;
    if !start.is_finite() || !min.is_finite() || !max.is_finite() || min > max {
        return Err(RootFindingError::InvalidParameter);
    }
    let mut x = start.clamp(min, max);
    let mut fx = f(x)?;
    if !fx.is_finite() {
        return Err(RootFindingError::NonFiniteComputation);
    }
    if fx.abs() <= eps {
        return Ok(x);
    }
    for _ in 0..max_iter {
        let dx = d(x)?;
        if !dx.is_finite() {
            return Err(RootFindingError::NonFiniteComputation);
        }
        // A small nonzero derivative can give a valid Newton step. Its units
        // differ from the coordinate-error tolerance, so compare only with zero.
        if dx == F::zero() {
            return Err(RootFindingError::ZeroDerivative);
        }
        let x1 = x - fx / dx;
        if !x1.is_finite() {
            return Err(RootFindingError::NonFiniteComputation);
        }
        let next = x1.clamp(min, max);
        if next == x {
            return Err(RootFindingError::FailedToConverge);
        }
        x = next;
        fx = f(x)?;
        if !fx.is_finite() {
            return Err(RootFindingError::NonFiniteComputation);
        }
        // Check even the final permitted update before exhausting the budget.
        if fx.abs() <= eps {
            return Ok(x);
        }
    }
    Err(RootFindingError::MaxIterationsReached)
}
