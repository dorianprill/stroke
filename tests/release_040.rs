//! Public-API regressions for the final 0.4.0 numerical/API pass.
#![allow(incomplete_features)]
#![feature(generic_const_exprs)]

use stroke::{
    ArrayPoint, BSpline, BSplineError, BSplinePath, Bezier, CubicBezier, LineSegment,
    QuadraticBezier,
};

#[test]
fn array_arithmetic_preserves_generic_components() {
    let a = ArrayPoint::from([1_i32, 2, 3]);
    let b = ArrayPoint::from([4, 5, 6]);
    assert_eq!((a + b).into_array(), [5, 7, 9]);
    assert_eq!((b - a).into_array(), [3, 3, 3]);
    assert_eq!((a * 3).into_array(), [3, 6, 9]);
    assert_eq!((a + 1).into_array(), [2, 3, 4]);
    assert_eq!((a - 1).into_array(), [0, 1, 2]);
    assert_eq!((ArrayPoint::from([] as [i32; 0]) * 3).into_array(), []);
    #[derive(Debug, PartialEq)]
    struct Owned(i32);
    impl core::ops::Mul<i32> for Owned {
        type Output = Self;
        fn mul(self, rhs: i32) -> Self {
            Self(self.0 * rhs)
        }
    }
    assert_eq!((ArrayPoint::from([Owned(2)]) * 3).into_array(), [Owned(6)]);
}

macro_rules! scalar_tests {
    ($module:ident, $float:ty, $huge:expr, $tiny:expr, $scales:expr, $knot_scales:expr) => {
        mod $module {
            use super::*;
            type F = $float;
            type P = ArrayPoint<F, 2>;

            fn near(actual: F, expected: F) {
                assert!(
                    actual.is_finite(),
                    "non-finite result {actual}, expected {expected}"
                );
                if expected == 0.0 {
                    assert_eq!(actual, 0.0);
                } else {
                    assert!(
                        (actual - expected).abs() / expected.abs() <= 128.0 * F::EPSILON,
                        "{actual} differs from {expected}"
                    );
                }
            }

            fn point_near(actual: P, expected: [F; 2]) {
                for axis in 0..2 {
                    near(actual[axis], expected[axis]);
                }
            }

            #[test]
            fn constant_and_linear_differentials_are_total() {
                let zero = P::from([0.0, 0.0]);
                let constant = Bezier::new([P::from([$huge, $huge])]).unwrap();
                assert_eq!(constant.derivative_at(0.5), zero);
                assert_eq!(constant.tangent(0.5), zero);
                assert_eq!(constant.curvature(0.5), 0.0);
                assert_eq!(constant.normal(0.5), None);
                let line = Bezier::new([zero, P::from([1.0, 0.0])]).unwrap();
                assert_eq!(line.derivative().eval(0.5), P::from([1.0, 0.0]));
                assert_eq!(line.curvature(0.5), 0.0);
                assert_eq!(line.normal(0.5), None);
                let spline =
                    BSpline::<P, 2, 1, 0>::new_clamped_uniform([constant.start()]).unwrap();
                assert_eq!(spline.tangent(0.5).unwrap(), zero);
                assert_eq!(spline.curvature(0.5).unwrap(), 0.0);
                assert_eq!(spline.normal(0.5).unwrap(), None);
                let mut path = BSplinePath::<P, 2, 1, 0, 1>::new();
                assert!(path.push(spline));
                assert_eq!(path.tangent(0.5), Some(zero));
                assert_eq!(path.curvature(0.5), Some(0.0));
                assert_eq!(path.normal(0.5), None);
            }

            #[test]
            fn specialized_interpolation_and_splits_preserve_finite_endpoints() {
                let p = P::from([$huge, $huge]);
                let constant = QuadraticBezier::new(p, p, p).unwrap();
                let cubic = CubicBezier::new(p, p, p, p).unwrap();
                for t in [0.0, 0.25, 0.5, 0.75, 1.0] {
                    assert_eq!(constant.eval(t), p);
                    assert_eq!(constant.eval_casteljau(t), p);
                    assert_eq!(cubic.eval(t), p);
                    assert_eq!(cubic.eval_casteljau(t), p);
                }
                let q = QuadraticBezier::new(
                    P::from([-$huge, 0.0]),
                    P::from([$huge, 0.0]),
                    P::from([0.0, 0.0]),
                )
                .unwrap();
                let c = CubicBezier::new(
                    P::from([-$huge, 0.0]),
                    P::from([$huge, 0.0]),
                    P::from([0.0, 0.0]),
                    P::from([0.0, 0.0]),
                )
                .unwrap();
                for t in [0.0, 0.25, 0.5, 0.75, 1.0] {
                    let qpoint = q.eval(t);
                    let cpoint = c.eval(t);
                    assert!(qpoint[0].is_finite() && cpoint[0].is_finite());
                    assert_eq!(q.eval_casteljau(t), qpoint);
                    assert_eq!(c.eval_casteljau(t), cpoint);
                    assert_eq!(q.split(t).0.end(), qpoint);
                    assert_eq!(q.split(t).1.start(), qpoint);
                    assert_eq!(c.split(t).0.end(), cpoint);
                    assert_eq!(c.split(t).1.start(), cpoint);
                }
                near(q.eval(0.5)[0], $huge * 0.25);
                near(c.eval(0.5)[0], $huge * 0.25);
                let cubic = CubicBezier::new(
                    P::from([1.0, 2.0]),
                    P::from([1.0, 2.0]),
                    P::from([1.0, 2.0]),
                    P::from([1.0, 2.0]),
                )
                .unwrap();
                assert_eq!(cubic.derivative_at(0.5), P::from([0.0, 0.0]));
            }

            #[test]
            fn tangents_survive_subnormals_and_overflowing_lengths() {
                let line = LineSegment::new(P::from([0.0, 0.0]), P::from([$tiny, 0.0])).unwrap();
                point_near(line.tangent(0.5), [1.0, 0.0]);
                let tiny = F::from_bits(1);
                let line = LineSegment::new(P::from([0.0, 0.0]), P::from([tiny, 0.0])).unwrap();
                point_near(line.tangent(0.5), [1.0, 0.0]);
                let line = LineSegment::new(P::from([-F::MAX, -F::MAX]), P::from([F::MAX, F::MAX]))
                    .unwrap();
                let unit = 1.0 / (2.0 as F).sqrt();
                point_near(line.tangent(0.5), [unit, unit]);
                point_near(line.reverse().tangent(0.5), [-unit, -unit]);
                let returning = Bezier::new([
                    P::from([-F::MAX, 0.0]),
                    P::from([F::MAX, 0.0]),
                    P::from([-F::MAX, 0.0]),
                ])
                .unwrap();
                assert_eq!(returning.derivative_at(0.5), P::from([0.0, 0.0]));
                assert_eq!(returning.tangent(0.5), P::from([0.0, 0.0]));
            }

            #[test]
            fn geometry_obeys_uniform_scale_and_reversal() {
                let unit = 1.0 / (2.0 as F).sqrt();
                for scale in $scales {
                    let controls = [
                        P::from([0.0, 0.0]),
                        P::from([scale, 0.0]),
                        P::from([scale, scale]),
                    ];
                    let q = QuadraticBezier::new(controls[0], controls[1], controls[2]).unwrap();
                    let generic = Bezier::new(controls).unwrap();
                    let cubic = CubicBezier::new(
                        controls[0],
                        P::from([scale * (2.0 / 3.0), 0.0]),
                        P::from([scale, scale / 3.0]),
                        controls[2],
                    )
                    .unwrap();
                    let spline =
                        BSpline::<P, 6, 3, 2>::new([-1.0, -1.0, -1.0, 3.0, 3.0, 3.0], controls)
                            .unwrap();
                    let expected_curvature = (2.0 as F).sqrt() / scale;
                    near(q.curvature(0.5), expected_curvature);
                    near(generic.curvature(0.5), expected_curvature);
                    near(cubic.curvature(0.5), expected_curvature);
                    near(spline.curvature(1.0).unwrap(), expected_curvature);
                    near(q.reverse().curvature(0.5), expected_curvature);
                    point_near(q.tangent(0.5), [unit, unit]);
                    point_near(cubic.tangent(0.5), [unit, unit]);
                    point_near(spline.tangent(1.0).unwrap(), [unit, unit]);
                    point_near(q.normal(0.5).unwrap(), [-unit, unit]);
                    point_near(cubic.normal(0.5).unwrap(), [-unit, unit]);
                    point_near(spline.normal(1.0).unwrap().unwrap(), [-unit, unit]);
                    point_near(q.reverse().normal(0.5).unwrap(), [-unit, unit]);
                }
                let stationary = Bezier::new([
                    P::from([0.0, 0.0]),
                    P::from([1.0, 1.0]),
                    P::from([0.0, 0.0]),
                ])
                .unwrap();
                assert_eq!(stationary.tangent(0.5), P::from([0.0, 0.0]));
                assert_eq!(stationary.curvature(0.5), 0.0);
                assert_eq!(stationary.normal(0.5), None);

                let scale: F = $huge;
                let q = QuadraticBezier::new(
                    P::from([0.0, 0.0]),
                    P::from([scale, 0.0]),
                    P::from([scale, scale]),
                )
                .unwrap();
                near(q.curvature(0.5), (2.0 as F).sqrt() / scale);
                point_near(q.normal(0.5).unwrap(), [-unit, unit]);
            }

            #[test]
            fn degeneracy_uses_distance_units_and_exact_zero() {
                for scale in $scales {
                    let q = QuadraticBezier::new(
                        P::from([0.0, 0.0]),
                        P::from([0.0, 0.0]),
                        P::from([scale, 0.0]),
                    )
                    .unwrap();
                    let c = CubicBezier::new(q.start(), q.start(), q.start(), q.end()).unwrap();
                    for tolerance in [0.0, scale * 0.5, -scale * 0.5] {
                        assert!(!q.is_a_point(tolerance));
                        assert!(!c.is_a_point(tolerance));
                        assert!(q.is_linear(tolerance));
                        assert!(c.is_linear(tolerance));
                    }
                    assert!(q.is_a_point(scale));
                    assert!(c.is_a_point(-scale));
                    assert!(!q.is_linear(scale));
                    for tolerance in [F::NAN, F::INFINITY, F::NEG_INFINITY] {
                        assert!(!q.is_a_point(tolerance));
                        assert!(!c.is_a_point(tolerance));
                        assert!(!q.is_linear(tolerance));
                        assert!(!c.is_linear(tolerance));
                    }
                }
                let tiny = F::from_bits(1);
                let q = QuadraticBezier::new(
                    P::from([0.0, 0.0]),
                    P::from([0.0, 0.0]),
                    P::from([tiny, 0.0]),
                )
                .unwrap();
                assert!(!q.is_a_point(0.0));
                assert!(q.is_linear(0.0));
                // Approximate coincidence uses the first point as its reference.
                let c = CubicBezier::new(
                    P::from([0.0, 0.0]),
                    P::from([0.0, 0.0]),
                    P::from([-0.75, 0.0]),
                    P::from([0.75, 0.0]),
                )
                .unwrap();
                assert!(c.is_a_point(1.0));
                assert!(!c.is_a_point(0.5));
            }

            #[test]
            fn cubic_basis_matches_analytic_derivatives_in_native_units() {
                let points = [
                    P::from([0.0, 0.0]),
                    P::from([1.0, 0.0]),
                    P::from([1.0, 1.0]),
                    P::from([0.0, 1.0]),
                ];
                for width in $knot_scales {
                    let spline = BSpline::<P, 8, 4, 3>::new(
                        [0.0, 0.0, 0.0, 0.0, width, width, width, width],
                        points,
                    )
                    .unwrap();
                    for u in [0.0, 0.125, 0.5, 0.875, 1.0] {
                        let v = 1.0 - u;
                        let expected = [
                            [v * v * v, 3.0 * u * v * v, 3.0 * u * u * v, u * u * u],
                            [
                                -3.0 * v * v,
                                3.0 * v * v - 6.0 * u * v,
                                6.0 * u * v - 3.0 * u * u,
                                3.0 * u * u,
                            ],
                            [6.0 * v, -12.0 + 18.0 * u, 6.0 - 18.0 * u, 6.0 * u],
                            [-6.0, 18.0, -18.0, 6.0],
                        ];
                        let t = width * u;
                        let table = spline.basis_functions_with_derivatives::<3>(3, t).unwrap();
                        assert_eq!(spline.basis_functions(3, t).unwrap(), table[0]);
                        for order in 0..=3 {
                            for j in 0..4 {
                                let mut native = expected[order][j];
                                for _ in 0..order {
                                    native /= width;
                                }
                                near(table[order][j], native);
                            }
                        }
                        let reconstructed = points
                            .into_iter()
                            .zip(table[0])
                            .fold(P::from([0.0, 0.0]), |sum, (p, n)| sum + p * n);
                        let actual = spline.eval(t).unwrap();
                        point_near(reconstructed, actual.into_array());
                    }
                }
                let wide = BSpline::<P, 4, 2, 1>::new(
                    [-$huge, -$huge, $huge, $huge],
                    [points[0], points[1]],
                )
                .unwrap();
                assert_eq!(wide.basis_functions(1, 0.0).unwrap(), [0.5, 0.5]);
                let table = wide.basis_functions_with_derivatives::<1>(1, 0.0).unwrap();
                near(table[1][0], -0.5 / ($huge as F));
                near(table[1][1], 0.5 / ($huge as F));
            }
        }
    };
}

scalar_tests!(
    double,
    f64,
    1e308,
    1e-310,
    [1e-300, 1e-100, 1.0, 1e100, 1e300],
    [1e-100, 1.0, 1e100]
);
scalar_tests!(
    single,
    f32,
    3e38,
    1e-40,
    [1e-30, 1e-10, 1.0, 1e10, 1e30],
    [1e-10, 1.0, 1e10]
);

#[test]
fn basis_rejects_invalid_spans_and_preserves_one_sided_ownership() {
    type P = ArrayPoint<f64, 1>;
    let points = [0.0, 1.0, 10.0, 9.0].map(|x| P::from([x]));
    let spline = BSpline::<P, 6, 4, 1>::new([0.0, 0.0, 0.5, 0.5, 1.0, 1.0], points).unwrap();
    for span in [0, 2, 4, 5, usize::MAX - 1, usize::MAX] {
        assert_eq!(
            spline.basis_functions(span, 0.25),
            Err(BSplineError::InvalidSpan)
        );
        assert_eq!(
            spline.basis_functions_with_derivatives::<1>(span, 0.25),
            Err(BSplineError::InvalidSpan)
        );
    }
    assert_eq!(
        spline.basis_functions(1, 0.75),
        Err(BSplineError::InvalidSpan)
    );
    assert_eq!(
        spline.basis_functions(1, 0.5),
        Err(BSplineError::InvalidSpan)
    );
    assert_eq!(spline.basis_functions(3, 0.5).unwrap(), [1.0, 0.0]);
    assert_eq!(spline.basis_functions(3, 1.0).unwrap(), [0.0, 1.0]);
    assert_eq!(spline.tangent(0.5).unwrap(), P::from([-1.0]));
    assert_eq!(spline.tangent(1.0).unwrap(), P::from([-1.0]));
    for t in [f64::NAN, f64::INFINITY, -1.0, 2.0] {
        assert_eq!(
            spline.basis_functions(1, t),
            Err(BSplineError::KnotDomainViolation)
        );
    }
    let spline = BSpline::<P, 2, 1, 0>::new([2.0, 3.0], [P::from([7.0])]).unwrap();
    assert_eq!(spline.basis_functions(0, 2.0).unwrap(), [1.0]);
    assert_eq!(
        spline
            .basis_functions_with_derivatives::<0>(0, 3.0)
            .unwrap(),
        [[1.0]]
    );
}

#[test]
fn basis_handles_overflowing_width_and_reports_unrepresentable_derivatives() {
    type P = ArrayPoint<f64, 1>;
    let points = [P::from([0.0]), P::from([1.0])];
    let spline = BSpline::<P, 4, 2, 1>::new([-1e308, -1e308, 1e308, 1e308], points).unwrap();
    assert_eq!(spline.basis_functions(1, 0.0).unwrap(), [0.5, 0.5]);
    let derivatives = spline
        .basis_functions_with_derivatives::<1>(1, 0.0)
        .unwrap();
    assert_eq!(derivatives[1], [-5e-309, 5e-309]);
    let spline = BSpline::<P, 4, 2, 1>::new([0.0, 0.0, 1e-310, 1e-310], points).unwrap();
    assert_eq!(spline.basis_functions(1, 0.0).unwrap(), [1.0, 0.0]);
    assert_eq!(
        spline.basis_functions_with_derivatives::<1>(1, 0.0),
        Err(BSplineError::NonFiniteComputation)
    );
}

#[test]
fn geometry_is_independent_of_knot_units() {
    type P = ArrayPoint<f64, 2>;
    let controls = [
        P::from([0.0, 0.0]),
        P::from([1.0, 0.0]),
        P::from([1.0, 1.0]),
    ];
    let q = Bezier::new(controls).unwrap();
    // Powers of two retain an exactly representable midpoint even below the
    // normal range; arbitrary decimal subnormals can have an odd ulp count.
    let tiny = f64::MIN_POSITIVE * 0.125;
    for (start, end, mid) in [(-1e308, 1e308, 0.0), (0.0, tiny, tiny * 0.5)] {
        let spline =
            BSpline::<P, 6, 3, 2>::new([start, start, start, end, end, end], controls).unwrap();
        assert_eq!(spline.tangent(mid).unwrap(), q.tangent(0.5));
        assert_eq!(spline.curvature(mid).unwrap(), q.curvature(0.5));
        assert_eq!(spline.normal(mid).unwrap(), q.normal(0.5));
    }
}

#[test]
fn almost_straight_and_collinear_curves_have_consistent_geometry() {
    type P = ArrayPoint<f64, 3>;
    let offset = 1e-12;
    let curve = Bezier::new([
        P::from([0.0, 0.0, 0.0]),
        P::from([0.5, 0.0, 0.0]),
        P::from([1.0, offset, 0.0]),
    ])
    .unwrap();
    let expected = 2.0 * offset / (1.0 + offset * offset).powf(1.5);
    assert!((curve.curvature(0.5) / expected - 1.0).abs() < 1e-14);
    let normal = curve.normal(0.5).unwrap();
    assert!((normal[0] + offset).abs() < 1e-26);
    assert!((normal[1] - 1.0).abs() < 1e-14);
    assert_eq!(normal[2], 0.0);
    let straight = Bezier::new([
        P::from([0.0, 0.0, 0.0]),
        P::from([1.0, 3.0, 0.0]),
        P::from([3.0, 9.0, 0.0]),
    ])
    .unwrap();
    assert_eq!(straight.curvature(0.5), 0.0);
    assert_eq!(straight.normal(0.5), None);

    let controls = [
        P::from([0.0, 0.0, 0.0]),
        P::from([1.0, 2.0, 0.0]),
        P::from([1.0, 2.0, 3.0]),
    ];
    let expected = (24.0_f64.hypot(12.0)) / 14.0_f64.powf(1.5);
    for scale in [1e-200, 1.0, 1e200] {
        let curve = Bezier::new(controls.map(|p| p * scale)).unwrap();
        assert!((curve.curvature(0.5) * scale / expected - 1.0).abs() < 1e-14);
        let normal = curve.normal(0.5).unwrap();
        let v = [1.0, 2.0, 3.0];
        let perpendicular = [-2.0 - 8.0 / 14.0, -4.0 - 16.0 / 14.0, 6.0 - 24.0 / 14.0];
        let norm = perpendicular.iter().fold(0.0_f64, |norm, &x| norm.hypot(x));
        let mut dot = 0.0;
        for i in 0..3 {
            assert!((normal[i] - perpendicular[i] / norm).abs() < 1e-14);
            dot += v[i] * normal[i];
        }
        assert!(dot.abs() < 1e-14);
    }
}

#[test]
fn basis_keeps_small_spans_next_to_distant_knots() {
    type P = ArrayPoint<f64, 1>;
    let points = [0.0, 1.0, 2.0, 3.0, 4.0].map(|x| P::from([x]));
    let tiny = f64::MIN_POSITIVE * 0.125;
    let curve =
        BSpline::<P, 9, 5, 3>::new([0.0, 0.0, 0.0, 0.0, tiny, 1.0, 1.0, 1.0, 1.0], points).unwrap();
    for t in [0.0, tiny * 0.5, tiny, 0.5, 1.0] {
        let span = curve.knot_span(t).unwrap();
        let basis = curve.basis_functions(span, t).unwrap();
        assert!((basis.iter().sum::<f64>() - 1.0).abs() < 1e-14);
        assert!(basis.iter().all(|&value| value >= 0.0));
        let value = basis
            .into_iter()
            .enumerate()
            .map(|(j, n)| points[span - 3 + j][0] * n)
            .sum::<f64>();
        assert!((value - curve.eval(t).unwrap()[0]).abs() < 1e-14);
    }
    // Uniform unclamped cubic basis at the start and midpoint of its domain.
    let curve = BSpline::<P, 8, 4, 3>::new(
        [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0],
        [points[0], points[1], points[2], points[3]],
    )
    .unwrap();
    for (t, coefficients) in [
        (3.0, [1.0, 4.0, 1.0, 0.0].map(|x| x / 6.0)),
        (3.5, [1.0, 23.0, 23.0, 1.0].map(|x| x / 48.0)),
    ] {
        let table = curve.basis_functions_with_derivatives::<3>(3, t).unwrap();
        for (actual, expected) in table[0].into_iter().zip(coefficients) {
            assert!((actual - expected).abs() < 1e-14);
        }
        for row in &table[1..] {
            assert!(row.iter().sum::<f64>().abs() < 1e-14);
        }
    }
}

#[cfg(feature = "nalgebra")]
#[test]
fn nalgebra_uses_the_same_scaled_geometry() {
    use nalgebra::SVector;
    type P = SVector<f64, 2>;
    let scale = 1e-200;
    let curve =
        QuadraticBezier::new(P::new(0.0, 0.0), P::new(scale, 0.0), P::new(scale, scale)).unwrap();
    assert!((curve.tangent(0.5).norm() - 1.0).abs() < 1e-14);
    assert!((curve.normal(0.5).unwrap().norm() - 1.0).abs() < 1e-14);
    assert!((curve.curvature(0.5) * scale - 2.0_f64.sqrt()).abs() < 1e-14);
}
