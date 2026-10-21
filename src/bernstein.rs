//! Bernstein weights and parameter derivatives, shared by the Bezier types.

use crate::error::BezierBasisError;
use crate::numeric::Scaled;
use num_traits::{Float, NumCast};

pub(crate) fn basis<T: Float, const N: usize, const R: usize>(
    t: T,
) -> Result<[[T; N]; R + 1], BezierBasisError>
where
    [(); R + 1]: Sized,
{
    if !t.is_finite() {
        return Err(BezierBasisError::NonFiniteParameter);
    }

    let zero = Scaled::new(T::zero());
    let one = Scaled::new(T::one());
    let parameter = Scaled::new(t);
    let complement = one.difference(parameter);
    let mut table = [[zero; N]; R + 1];
    // A constructed Bezier has at least one control point.
    table[0][0] = one;

    for degree in 1..N {
        let mut next = [[zero; N]; R + 1];
        for i in 0..=degree {
            let left = if i == 0 { zero } else { table[0][i - 1] };
            let right = if i == degree { zero } else { table[0][i] };
            next[0][i] = left.product(parameter).sum(right.product(complement));
        }
        if R > 0 {
            let factor = Scaled::new(
                <T as NumCast>::from(degree).ok_or(BezierBasisError::NonFiniteComputation)?,
            );
            // d^r B(i,n) = n * (d^(r-1) B(i-1,n-1) - d^(r-1) B(i,n-1)).
            // Keep separate exponents until the final degree, since a small
            // lower-degree weight can still give a representable derivative.
            for order in 1..=R.min(degree) {
                for i in 0..=degree {
                    let left = if i == 0 {
                        zero
                    } else {
                        table[order - 1][i - 1]
                    };
                    let right = if i == degree {
                        zero
                    } else {
                        table[order - 1][i]
                    };
                    next[order][i] = left.difference(right).product(factor);
                }
            }
        }
        table = next;
    }

    let result = table.map(|row| row.map(Scaled::value));
    if result.iter().flatten().all(|value| value.is_finite()) {
        Ok(result)
    } else {
        Err(BezierBasisError::NonFiniteComputation)
    }
}
