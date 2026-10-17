//! Regressions for residual-checked roots and scale-safe distance queries.

use stroke::{ArrayPoint, BSpline, Bezier, CubicBezier, LineSegment, QuadraticBezier};

macro_rules! numeric_queries {
    ($root:ident, $candidates:ident, $line:ident, $scalar:ty, $large:expr, $small:expr, $eps:expr) => {
        #[test]
        fn $root() {
            type P = ArrayPoint<$scalar, 1>;
            let points = [P::from([0.0]), P::from([0.0]), P::from([1e12])];
            let curve = Bezier::new(points).unwrap();
            let spline = BSpline::<P, 6, 3, 2>::new_clamped_uniform(points).unwrap();
            assert!(
                curve
                    .root_newton_axis(-1.0, 0, 1e-6, Some(1e-6), None)
                    .is_err()
            );
            assert!(
                spline
                    .root_newton_axis(-1.0, 0, 1e-6, Some(1e-6), None)
                    .is_err()
            );
            // A small parameter step on a steep curve is not evidence of a small residual.
            for target in [0.0, 1.0, 4.0] {
                let t = curve
                    .root_newton_axis(target, 0, 1e-6, Some($eps), Some(128))
                    .unwrap();
                assert!((curve.eval(t)[0] - target).abs() <= $eps);
                let t = spline
                    .root_newton_axis(target, 0, 1e-6, Some($eps), Some(128))
                    .unwrap();
                assert!((spline.eval(t).unwrap()[0] - target).abs() <= $eps);
            }
        }

        #[test]
        fn $candidates() {
            type P = ArrayPoint<$scalar, 2>;
            let quadratic = QuadraticBezier::new(
                P::from([0.0, 0.0]),
                P::from([1.0, 0.0]),
                P::from([2.0, 0.0]),
            )
            .unwrap();
            let cubic = CubicBezier::new(
                P::from([0.0, 0.0]),
                P::from([1.0, 0.0]),
                P::from([2.0, 0.0]),
                P::from([3.0, 0.0]),
            )
            .unwrap();
            let generic = Bezier::from(quadratic);
            let spline =
                BSpline::<P, 6, 3, 2>::new_clamped_uniform(quadratic.into_control_points())
                    .unwrap();
            for t in [0.0, 0.25, 0.5, 0.75, 1.0] {
                assert_eq!(generic.distance_to_point(generic.eval(t)), 0.0);
                assert_eq!(quadratic.distance_to_point(quadratic.eval(t)), 0.0);
                assert_eq!(cubic.distance_to_point(cubic.eval(t)), 0.0);
                assert_eq!(
                    spline.distance_to_point(spline.eval(t).unwrap()).unwrap(),
                    0.0
                );
            }
            for scale in [$large, $small] {
                let point = P::from([scale, scale]);
                let curve = Bezier::new([P::from([0.0, 0.0])]).unwrap();
                let actual = curve.distance_to_point(point);
                assert!(actual.is_finite() && actual > 0.0);
                assert!((actual / scale - (2.0 as $scalar).sqrt()).abs() <= $eps);
            }
        }

        #[test]
        fn $line() {
            type P = ArrayPoint<$scalar, 2>;
            for scale in [$large, $small, 1.0] {
                let segment = LineSegment::new(P::from([0.0, 0.0]), P::from([scale, 0.0])).unwrap();
                for line in [segment, segment.reverse()] {
                    assert_eq!(line.distance_to_point(P::from([0.0, 0.0])), 0.0);
                    assert_eq!(line.distance_to_point(P::from([scale, 0.0])), 0.0);
                    for (point, expected) in [
                        (P::from([scale * 0.5, scale]), scale),
                        (P::from([-scale, 0.0]), scale),
                        (P::from([scale * 2.0, 0.0]), scale),
                    ] {
                        let actual = line.distance_to_point(point);
                        assert!(actual.is_finite() && actual > 0.0);
                        assert!(
                            (actual / expected - 1.0).abs() <= $eps,
                            "{actual} != {expected}"
                        );
                    }
                }
            }
            let extreme = <$scalar>::MAX;
            let line = LineSegment::new(P::from([-extreme, 0.0]), P::from([extreme, 0.0])).unwrap();
            assert_eq!(line.distance_to_point(P::from([0.0, 1.0])), 1.0);
            assert_eq!(line.reverse().distance_to_point(P::from([0.0, 1.0])), 1.0);
        }
    };
}

numeric_queries!(
    roots_f32,
    candidates_f32,
    lines_f32,
    f32,
    1e30_f32,
    1e-30_f32,
    2e-6_f32
);
numeric_queries!(
    roots_f64,
    candidates_f64,
    lines_f64,
    f64,
    1e200_f64,
    1e-200_f64,
    1e-12_f64
);
