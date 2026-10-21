#![allow(incomplete_features)]
#![feature(generic_const_exprs)]

use core::{array, fmt::Debug};
use num_traits::{Float, NumCast};
use stroke::{
    ArrayPoint, BSpline, Bezier, BezierBasisError, CubicBezier, LineSegment, QuadraticBezier,
};

fn scalar<T: Float>(value: f64) -> T {
    <T as NumCast>::from(value).unwrap()
}

fn curve<T: Float, const N: usize>() -> Bezier<ArrayPoint<T, 2>, N> {
    Bezier::new(array::from_fn(|i| {
        ArrayPoint::from([scalar(i as f64), scalar((i * i) as f64)])
    }))
    .unwrap()
}

#[track_caller]
fn near<T: Float + Debug>(actual: T, expected: T) {
    let tolerance = scalar::<T>(64.0) * T::epsilon() * expected.abs().max(T::one());
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:?} != {expected:?}"
    );
}

fn cubic_formula<T: Float + Debug>() {
    let curve = curve::<T, 4>();
    for t in [-2.0, -0.25, 0.0, 0.125, 0.5, 0.875, 1.0, 1.25, 3.0] {
        let t = scalar::<T>(t);
        let s = T::one() - t;
        let three = scalar::<T>(3.0);
        let six = scalar::<T>(6.0);
        let twelve = scalar::<T>(12.0);
        let eighteen = scalar::<T>(18.0);
        let expected = [
            [s * s * s, three * s * s * t, three * s * t * t, t * t * t],
            [
                -three * s * s,
                three * s * s - six * s * t,
                six * s * t - three * t * t,
                three * t * t,
            ],
            [six * s, six * t - twelve * s, six * s - twelve * t, six * t],
            [-six, eighteen, -eighteen, six],
            [T::zero(); 4],
            [T::zero(); 4],
        ];
        let actual = curve.basis_functions_with_derivatives::<5>(t).unwrap();
        assert_eq!(actual[0], curve.basis_functions(t).unwrap());
        assert_eq!(
            curve.basis_functions_with_derivatives::<0>(t).unwrap(),
            [actual[0]]
        );
        for (actual, expected) in actual
            .into_iter()
            .flatten()
            .zip(expected.into_iter().flatten())
        {
            near(actual, expected);
        }
        // The highest derivative is constant, even while extrapolating.
        assert_eq!(actual[3], [-six, eighteen, -eighteen, six]);
        assert_eq!(actual[4], [T::zero(); 4]);
        assert_eq!(actual[5], [T::zero(); 4]);
    }
}

#[test]
fn analytic_cubic_inside_and_outside_domain_f32_and_f64() {
    cubic_formula::<f32>();
    cubic_formula::<f64>();
}

fn endpoints<T: Float + Debug, const N: usize>() {
    let curve = curve::<T, N>();
    for (t, endpoint) in [(T::zero(), 0), (-T::zero(), 0), (T::one(), N - 1)] {
        let actual = curve.basis_functions(t).unwrap();
        assert_eq!(
            actual,
            array::from_fn(|i| if i == endpoint { T::one() } else { T::zero() })
        );
    }
}

#[test]
fn endpoint_weights_are_exact_for_constants_and_higher_degrees() {
    endpoints::<f32, 1>();
    endpoints::<f32, 4>();
    endpoints::<f32, 151>();
    endpoints::<f64, 1>();
    endpoints::<f64, 4>();
    endpoints::<f64, 65>();
}

fn polynomial_reproduction<T: Float + Debug, const R: usize>()
where
    [(); R + 1]: Sized,
{
    const N: usize = 10;
    let curve = curve::<T, N>();
    for t in [-0.5, 0.0, 0.125, 0.5, 0.875, 1.0, 1.5] {
        let t = scalar::<T>(t);
        let basis = curve.basis_functions_with_derivatives::<R>(t).unwrap();
        for power in 0..=R {
            // Bernstein coefficients of t^power: (i)_power / (degree)_power.
            let coefficients: [T; N] = array::from_fn(|i| {
                if i < power {
                    return T::zero();
                }
                (0..power).fold(T::one(), |coefficient, k| {
                    coefficient * scalar::<T>((i - k) as f64) / scalar::<T>((N - 1 - k) as f64)
                })
            });
            for (order, row) in basis.iter().enumerate() {
                let expected = if order > power {
                    T::zero()
                } else {
                    (0..order).fold(t.powi((power - order) as i32), |v, k| {
                        v * scalar::<T>((power - k) as f64)
                    })
                };
                let (actual, magnitude) = row.iter().zip(coefficients).fold(
                    (T::zero(), T::zero()),
                    |(sum, magnitude), (&weight, coefficient)| {
                        let term = weight * coefficient;
                        (sum + term, magnitude + term.abs())
                    },
                );
                // Relative error is ill-conditioned at zeros and under cancellation.
                // Bound absolute error using the magnitudes of the summed terms.
                let tolerance = scalar::<T>((64 * N) as f64)
                    * T::epsilon()
                    * magnitude.max(expected.abs()).max(T::one());
                assert!(
                    (actual - expected).abs() <= tolerance,
                    "t={t:?}, power={power}, order={order}, actual={actual:?}, expected={expected:?}"
                );
            }
        }
        if t >= T::zero() && t <= T::one() {
            assert!(basis[0].iter().all(|&weight| weight >= T::zero()));
            near(
                basis[0].iter().copied().fold(T::zero(), |a, b| a + b),
                T::one(),
            );
        }
    }
}

#[test]
fn basis_reproduces_polynomials_and_their_derivatives() {
    polynomial_reproduction::<f32, 5>();
    polynomial_reproduction::<f64, 5>();
}

#[test]
fn weights_reconstruct_values_and_derivatives() {
    let controls = [[-2.0, 3.0], [4.0, -1.0], [2.0, 5.0], [-3.0, -2.0]].map(ArrayPoint::from);
    let curve = Bezier::new(controls).unwrap();
    let first = curve.derivative();
    let second = first.derivative();
    let third = second.derivative();
    for t in [-0.5, 0.0, 0.3, 0.5, 1.0, 1.5] {
        let table = curve.basis_functions_with_derivatives::<4>(t).unwrap();
        let expected = [
            curve.eval(t),
            first.eval(t),
            second.eval(t),
            third.eval(t),
            ArrayPoint::from([0.0; 2]),
        ];
        for (weights, expected) in table.into_iter().zip(expected) {
            let actual = weights
                .into_iter()
                .zip(controls)
                .fold(ArrayPoint::from([0.0; 2]), |p, (w, control)| {
                    p + control * w
                });
            for axis in 0..2 {
                near(actual[axis], expected[axis]);
            }
        }
    }
}

fn spline_equivalence<T: Float + Debug>() {
    let curve = curve::<T, 4>();
    let spline = BSpline::<_, 8, 4, 3>::new_clamped_uniform(curve.into_control_points()).unwrap();
    for t in [0.0, 0.1, 0.5, 0.9, 1.0] {
        let t = scalar::<T>(t);
        let span = spline.knot_span(t).unwrap();
        let bernstein = curve.basis_functions_with_derivatives::<3>(t).unwrap();
        let bspline = spline
            .basis_functions_with_derivatives::<3>(span, t)
            .unwrap();
        for (a, b) in bernstein
            .into_iter()
            .flatten()
            .zip(bspline.into_iter().flatten())
        {
            near(a, b);
        }
    }
}

#[test]
fn single_span_bspline_has_the_same_basis_and_derivatives() {
    spline_equivalence::<f32>();
    spline_equivalence::<f64>();
}

#[test]
fn specializations_share_values_derivatives_and_errors() {
    let generic_line = curve::<f64, 2>();
    let generic_quadratic = curve::<f64, 3>();
    let generic_cubic = curve::<f64, 4>();
    let line = LineSegment::from(generic_line);
    let quadratic = QuadraticBezier::from(generic_quadratic);
    let cubic = CubicBezier::from(generic_cubic);
    for t in [-2.0, 0.0, 0.125, 1.0, 3.0] {
        assert_eq!(line.basis_functions(t), generic_line.basis_functions(t));
        assert_eq!(
            quadratic.basis_functions(t),
            generic_quadratic.basis_functions(t)
        );
        assert_eq!(cubic.basis_functions(t), generic_cubic.basis_functions(t));
        assert_eq!(
            line.basis_functions_with_derivatives::<5>(t),
            generic_line.basis_functions_with_derivatives::<5>(t)
        );
        assert_eq!(
            quadratic.basis_functions_with_derivatives::<5>(t),
            generic_quadratic.basis_functions_with_derivatives::<5>(t)
        );
        assert_eq!(
            cubic.basis_functions_with_derivatives::<5>(t),
            generic_cubic.basis_functions_with_derivatives::<5>(t)
        );
    }
    for t in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert_eq!(
            line.basis_functions(t),
            Err(BezierBasisError::NonFiniteParameter)
        );
        assert_eq!(
            quadratic.basis_functions(t),
            Err(BezierBasisError::NonFiniteParameter)
        );
        assert_eq!(
            cubic.basis_functions(t),
            Err(BezierBasisError::NonFiniteParameter)
        );
        assert_eq!(
            line.basis_functions_with_derivatives::<5>(t),
            Err(BezierBasisError::NonFiniteParameter)
        );
        assert_eq!(
            quadratic.basis_functions_with_derivatives::<5>(t),
            Err(BezierBasisError::NonFiniteParameter)
        );
        assert_eq!(
            cubic.basis_functions_with_derivatives::<5>(t),
            Err(BezierBasisError::NonFiniteParameter)
        );
    }
}

#[test]
fn constant_basis_and_derivatives_hold_at_any_finite_parameter() {
    let curve = curve::<f64, 1>();
    for t in [-f64::MAX, -2.0, 0.0, 1.0, f64::MAX] {
        assert_eq!(curve.basis_functions(t).unwrap(), [1.0]);
        assert_eq!(
            curve.basis_functions_with_derivatives::<4>(t).unwrap(),
            [[1.0], [0.0], [0.0], [0.0], [0.0]]
        );
    }
    for t in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert_eq!(
            curve.basis_functions(t),
            Err(BezierBasisError::NonFiniteParameter)
        );
        assert_eq!(
            curve.basis_functions_with_derivatives::<4>(t),
            Err(BezierBasisError::NonFiniteParameter)
        );
    }
}

#[test]
fn finite_extrapolation_reports_unrepresentable_weights() {
    let line = curve::<f64, 2>();
    for t in [-f64::MAX, f64::MAX] {
        assert!(
            line.basis_functions(t)
                .unwrap()
                .into_iter()
                .all(f64::is_finite)
        );
        assert_eq!(
            line.basis_functions_with_derivatives::<2>(t).unwrap()[1],
            [-1.0, 1.0]
        );
    }
    let quadratic = curve::<f64, 3>();
    for t in [-f64::MAX, -1e200, 1e200, f64::MAX] {
        assert_eq!(
            quadratic.basis_functions(t),
            Err(BezierBasisError::NonFiniteComputation)
        );
    }
    assert!(
        curve::<f64, 4>()
            .basis_functions(1e100)
            .unwrap()
            .into_iter()
            .all(f64::is_finite)
    );
}

#[test]
fn unrepresentable_high_order_derivatives_return_errors() {
    let curve = curve::<f32, 33>();
    assert!(curve.basis_functions(0.0).is_ok());
    assert!(curve.basis_functions_with_derivatives::<1>(0.0).is_ok());
    assert_eq!(
        curve.basis_functions_with_derivatives::<32>(0.0),
        Err(BezierBasisError::NonFiniteComputation)
    );
}

#[test]
fn scaled_recurrence_preserves_representable_subnormal_derivatives() {
    let curve = curve::<f32, 151>();
    let table = curve.basis_functions_with_derivatives::<1>(0.5).unwrap();
    let smallest = f32::from_bits(1);
    assert_eq!(table[0][0], 0.0); // 2^-150 rounds below the smallest subnormal.
    assert_eq!(table[0][1], 75.0 * smallest);
    assert_eq!(table[1][0], -150.0 * smallest);
    assert_eq!(table[1][150], 150.0 * smallest);
    near(table[0].iter().sum(), 1.0_f32);
    assert!(table[0][75] > 0.0); // The unscaled central binomial coefficient overflows f32.
}

#[test]
fn subnormal_parameters_keep_representable_weights() {
    let tiny = f32::from_bits(1);
    let curve = curve::<f32, 3>();
    assert_eq!(curve.basis_functions(tiny).unwrap(), [1.0, 2.0 * tiny, 0.0]);
    assert_eq!(
        curve.basis_functions_with_derivatives::<1>(tiny).unwrap()[1][2],
        2.0 * tiny
    );
}

#[test]
fn basis_errors_support_standard_error_handling() {
    for error in [
        BezierBasisError::NonFiniteParameter,
        BezierBasisError::NonFiniteComputation,
    ] {
        let standard: &dyn core::error::Error = &error;
        assert!(!standard.to_string().is_empty());
        assert!(standard.source().is_none());
    }
}

#[test]
#[cfg(feature = "nalgebra")]
fn nalgebra_points_use_the_same_basis() {
    use nalgebra::SVector;
    let generic = curve::<f64, 4>();
    let curve = Bezier::new(
        generic
            .into_control_points()
            .map(|p| SVector::<f64, 2>::from_row_slice(p.as_array())),
    )
    .unwrap();
    assert_eq!(curve.basis_functions(-0.5), generic.basis_functions(-0.5));
    assert_eq!(
        curve.basis_functions_with_derivatives::<5>(1.5),
        generic.basis_functions_with_derivatives::<5>(1.5)
    );
}
