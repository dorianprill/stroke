//! Input contracts and edge cases for local root and distance queries.

use stroke::{
    BSpline, Bezier, CubicBezier, FindRoot, LineSegment, PointN, PointNorm, RootFindingError,
};

#[test]
fn newton_validates_inputs_and_checks_the_last_update() {
    type P = PointN<f64, 1>;
    let curve = Bezier::new([P::from([0.0]), P::from([1.0])]).unwrap();
    let spline = BSpline::<P, 4, 2, 1>::new_clamped_uniform(curve.into_control_points()).unwrap();
    for accuracy in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert_eq!(
            curve.root_newton_axis(0.5, 0, 0.0, Some(accuracy), None),
            Err(RootFindingError::InvalidAccuracy)
        );
        assert_eq!(
            spline.root_newton_axis(0.5, 0, 0.0, Some(accuracy), None),
            Err(RootFindingError::InvalidAccuracy)
        );
    }
    for invalid in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert_eq!(
            curve.root_newton_axis(invalid, 0, 0.0, None, None),
            Err(RootFindingError::InvalidParameter)
        );
        assert_eq!(
            curve.root_newton_axis(0.5, 0, invalid, None, None),
            Err(RootFindingError::InvalidParameter)
        );
        assert_eq!(
            spline.root_newton_axis(invalid, 0, 0.0, None, None),
            Err(RootFindingError::InvalidParameter)
        );
        assert_eq!(
            spline.root_newton_axis(0.5, 0, invalid, None, None),
            Err(RootFindingError::InvalidParameter)
        );
    }
    assert_eq!(
        curve.root_newton_axis(0.5, 1, 0.0, None, None),
        Err(RootFindingError::InvalidAxis)
    );
    assert_eq!(
        spline.root_newton_axis(0.5, 1, 0.0, None, None),
        Err(RootFindingError::InvalidAxis)
    );
    assert_eq!(curve.axis_value(0.5, 1), Err(RootFindingError::InvalidAxis));
    assert_eq!(
        curve.axis_derivative(0.5, 1),
        Err(RootFindingError::InvalidAxis)
    );
    assert_eq!(
        spline.axis_value(0.5, 1),
        Err(RootFindingError::InvalidAxis)
    );
    assert_eq!(
        spline.axis_derivative(0.5, 1),
        Err(RootFindingError::InvalidAxis)
    );
    assert_eq!(
        curve.root_newton_axis(0.5, 0, 0.0, None, Some(0)),
        Err(RootFindingError::InvalidIterationLimit)
    );
    assert_eq!(
        spline.root_newton_axis(0.5, 0, 0.0, None, Some(0)),
        Err(RootFindingError::InvalidIterationLimit)
    );
    for guess in [-1.0, 0.0, 2.0] {
        assert_eq!(
            curve.root_newton_axis(0.5, 0, guess, None, Some(1)),
            Ok(0.5)
        );
        assert_eq!(
            spline.root_newton_axis(0.5, 0, guess, None, Some(1)),
            Ok(0.5)
        );
    }
    assert_eq!(
        curve.root_newton_axis(2.0, 0, 0.5, None, None),
        Err(RootFindingError::FailedToConverge)
    );
    assert_eq!(
        spline.root_newton_axis(2.0, 0, 0.5, None, None),
        Err(RootFindingError::FailedToConverge)
    );
}

#[test]
fn newton_keeps_coordinate_and_derivative_units_separate() {
    type P = PointN<f64, 1>;
    let curve = Bezier::new([P::from([0.0]), P::from([1e-12])]).unwrap();
    let spline = BSpline::<P, 4, 2, 1>::new(
        [2.0, 2.0, 1000002.0, 1000002.0],
        curve.into_control_points(),
    )
    .unwrap();
    assert_eq!(
        curve.root_newton_axis(5e-13, 0, 0.0, Some(1e-14), Some(1)),
        Ok(0.5)
    );
    let t = spline
        .root_newton_axis(5e-13, 0, 2.0, Some(1e-14), Some(1))
        .unwrap();
    assert!((spline.eval(t).unwrap()[0] - 5e-13).abs() <= 1e-14);
    let constant = Bezier::new([P::from([7.0])]).unwrap();
    assert_eq!(
        constant.root_newton_axis(7.0, 0, 0.25, None, None),
        Ok(0.25)
    );
    assert_eq!(
        constant.root_newton_axis(8.0, 0, 0.25, None, None),
        Err(RootFindingError::ZeroDerivative)
    );
    let constant = BSpline::<P, 2, 1, 0>::new_clamped_uniform([P::from([7.0])]).unwrap();
    assert_eq!(
        constant.root_newton_axis(7.0, 0, 0.25, None, None),
        Ok(0.25)
    );
    assert_eq!(
        constant.root_newton_axis(8.0, 0, 0.25, None, None),
        Err(RootFindingError::ZeroDerivative)
    );
    let overflowing_derivative = Bezier::new([P::from([-f64::MAX]), P::from([f64::MAX])]).unwrap();
    assert_eq!(
        overflowing_derivative.root_newton_axis(0.0, 0, 0.0, Some(1e-12), None),
        Err(RootFindingError::NonFiniteComputation)
    );
}

#[test]
fn newton_stagnation_and_discontinuities_cannot_fake_a_root() {
    type P = PointN<f64, 1>;
    let curve = Bezier::new([P::from([1.0]), P::from([1e20])]).unwrap();
    assert_eq!(
        curve.root_newton_axis(0.0, 0, 0.0, Some(1e-12), None),
        Err(RootFindingError::FailedToConverge)
    );
    let jump = BSpline::<P, 6, 4, 1>::new(
        [0.0, 0.0, 0.5, 0.5, 1.0, 1.0],
        [
            P::from([0.0]),
            P::from([1.0]),
            P::from([3.0]),
            P::from([4.0]),
        ],
    )
    .unwrap();
    assert!(
        jump.root_newton_axis(2.0, 0, 0.25, Some(1e-12), Some(12))
            .is_err()
    );
}

#[test]
fn refinement_never_worsens_the_best_coarse_candidate() {
    type P = PointN<f64, 2>;
    let curve = Bezier::new([
        P::from([0.0, 0.0]),
        P::from([4.0, 5.0]),
        P::from([-3.0, 4.0]),
        P::from([1.0, -1.0]),
    ])
    .unwrap();
    let cubic = CubicBezier::from(curve);
    let spline = BSpline::<P, 8, 4, 3>::new_clamped_uniform(curve.into_control_points()).unwrap();
    for query in [
        P::from([-1.0, 2.0]),
        P::from([5.0, 0.0]),
        P::from([0.0, 0.0]),
    ] {
        for segments in [0, 1, 2, 3, 8, 64] {
            let count = segments.max(1);
            let best = (0..=count)
                .map(|i| (curve.eval(i as f64 / count as f64) - query).norm())
                .fold(f64::INFINITY, f64::min);
            assert!(curve.distance_to_point_approx(query, segments) <= best);
            assert!(cubic.distance_to_point_approx(query, segments) <= best);
            assert!(
                spline.distance_to_point_approx(query, segments).unwrap()
                    <= best + f64::EPSILON * 16.0
            );
        }
    }
}

macro_rules! projection_edges {
    ($name:ident, $scalar:ty, $accuracy:expr) => {
        #[test]
        fn $name() {
            type P = PointN<$scalar, 2>;
            let max = <$scalar>::MAX;
            let line = LineSegment::new(P::from([-max, -max]), P::from([max, max])).unwrap();
            assert_eq!(line.distance_to_point(P::from([1.0, 1.0])), 0.0);
            assert!(
                (line.distance_to_point(P::from([1.0, 2.0])) - (0.5 as $scalar).sqrt()).abs()
                    <= $accuracy
            );
            assert_eq!(line.distance_to_point(P::from([max * 0.5, max * 0.5])), 0.0);
            let eps = <$scalar>::EPSILON;
            let short =
                LineSegment::new(P::from([1.0, 1.0]), P::from([1.0 + eps, 1.0 + eps])).unwrap();
            let distance = short.distance_to_point(P::from([1.0, 1.0 + eps]));
            assert!((distance / eps - (0.5 as $scalar).sqrt()).abs() <= $accuracy);
            let tiny = <$scalar>::from_bits(1);
            let tiny_line =
                LineSegment::new(P::from([0.0, 0.0]), P::from([tiny * 4.0, 0.0])).unwrap();
            assert_eq!(
                tiny_line.distance_to_point(P::from([tiny, tiny * 2.0])),
                tiny * 2.0
            );
            let point = LineSegment::new(P::from([0.0, 0.0]), P::from([0.0, 0.0])).unwrap();
            assert_eq!(point.distance_to_point(P::from([3.0, 4.0])), 5.0);
            for invalid in [<$scalar>::NAN, <$scalar>::INFINITY, <$scalar>::NEG_INFINITY] {
                assert!(line.distance_to_point(P::from([invalid, 0.0])).is_nan());
                let curve = Bezier::new([P::from([0.0, 0.0])]).unwrap();
                assert!(curve.distance_to_point(P::from([invalid, 0.0])).is_nan());
            }
        }
    };
}

projection_edges!(projection_edges_f32, f32, 3e-6_f32);
projection_edges!(projection_edges_f64, f64, 2e-14_f64);

#[test]
fn projection_matches_reference_geometry_at_multiple_scales() {
    type P = PointN<f64, 3>;
    let a = P::from([-3.25, 1.5, 4.0]);
    let b = P::from([7.75, 3.5, -2.0]);
    let direction = b - a;
    for i in -10..=10 {
        for j in -10..=10 {
            let query = P::from([i as f64 * 0.25, j as f64 * 0.5, (i - j) as f64 * 0.3]);
            let dot = (0..3)
                .map(|axis| (query[axis] - a[axis]) * direction[axis])
                .sum::<f64>();
            let t = (dot / direction.squared_norm()).clamp(0.0, 1.0);
            let expected = (query - (a + direction * t)).norm();
            for scale in [1e-200, 1.0, 1e200] {
                let line = LineSegment::new(a * scale, b * scale).unwrap();
                for segment in [line, line.reverse()] {
                    let actual = segment.distance_to_point(query * scale) / scale;
                    assert!(
                        (actual - expected).abs() <= 128.0 * f64::EPSILON * expected.max(1.0),
                        "{actual} != {expected}"
                    );
                }
            }
        }
    }
}

#[cfg(feature = "nalgebra")]
#[test]
fn query_fixes_support_nalgebra_and_mixed_coordinate_scales() {
    use nalgebra::SVector;
    type V = SVector<f64, 3>;
    let line = LineSegment::new(V::new(0.0, 0.0, 0.0), V::new(1e200, 1e200, 0.0)).unwrap();
    assert_eq!(line.distance_to_point(V::new(5e199, 5e199, 1e-200)), 1e-200);
    let curve = Bezier::new([
        V::new(0.0, 0.0, 0.0),
        V::new(0.0, 0.0, 0.0),
        V::new(1e12, 0.0, 0.0),
    ])
    .unwrap();
    assert!(
        curve
            .root_newton_axis(-1.0, 0, 1e-6, Some(1e-6), None)
            .is_err()
    );
    assert_eq!(curve.distance_to_point(curve.eval(0.5)), 0.0);
}
