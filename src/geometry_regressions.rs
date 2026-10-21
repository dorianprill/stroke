//! Regression coverage for representable path boundaries and coordinate ranges.

use crate::test_support::assert_approx_equal;
use crate::{ArcLengthError, ArrayPoint, BSpline, BSplinePath, Bezier, BezierPath, LineSegment};

macro_rules! path_boundary_test {
    ($name:ident, $scalar:ty, $accuracy:expr) => {
        #[test]
        fn $name() {
            type P = ArrayPoint<$scalar, 1>;
            for count in 2..=32 {
                let mut path = BezierPath::<P, 32>::new();
                let mut bspline_path = BSplinePath::<P, 4, 2, 1, 32>::new();
                for i in 0..count {
                    let start = P::new([i as $scalar * 100.0]);
                    let end = P::new([i as $scalar * 100.0 + 1.0]);
                    assert!(path.push_line(LineSegment::new(start, end).unwrap()));
                    assert!(
                        bspline_path
                            .push(BSpline::new([0.0, 0.0, 1.0, 1.0], [start, end]).unwrap(),)
                    );
                }
                for index in 1..count {
                    let boundary = index as $scalar / count as $scalar;
                    let expected = P::new([index as $scalar * 100.0]);
                    assert_eq!(path.eval(boundary), Some(expected));
                    assert_eq!(bspline_path.eval(boundary), Some(expected));
                    let s = index as $scalar;
                    assert_eq!(path.point_at_length(s, $accuracy).unwrap(), Some(expected));
                    assert_eq!(
                        bspline_path.point_at_length(s, $accuracy).unwrap(),
                        Some(expected)
                    );
                    for t in [
                        path.t_at_length(s, $accuracy).unwrap().unwrap(),
                        bspline_path.t_at_length(s, $accuracy).unwrap().unwrap(),
                    ] {
                        assert_eq!(t, boundary);
                        assert_eq!(path.eval(t), Some(expected));
                        assert_eq!(bspline_path.eval(t), Some(expected));
                    }
                    // The nearest representable parameters on each side must
                    // stay on that side of the jump.
                    let left_min = (index - 1) as $scalar * 100.0;
                    for point in [
                        path.eval(boundary.next_down()),
                        bspline_path.eval(boundary.next_down()),
                    ] {
                        let x = point.unwrap()[0];
                        assert!(x >= left_min && x <= left_min + 1.0);
                    }
                    for point in [
                        path.eval(boundary.next_up()),
                        bspline_path.eval(boundary.next_up()),
                    ] {
                        let x = point.unwrap()[0];
                        assert!(x >= expected[0] && x <= expected[0] + 1.0);
                    }
                    // Direct length lookup must select the same segment even
                    // when the path parameter loses segment-local precision.
                    for distance in [s.next_down(), s.next_up()] {
                        let segment_start = if distance < s { left_min } else { expected[0] };
                        for point in [
                            path.point_at_length(distance, $accuracy).unwrap(),
                            bspline_path.point_at_length(distance, $accuracy).unwrap(),
                        ] {
                            let x = point.unwrap()[0];
                            assert!(x >= segment_start && x <= segment_start + 1.0);
                        }
                        let t = path.t_at_length(distance, $accuracy).unwrap().unwrap();
                        let x = path.eval(t).unwrap()[0];
                        assert!(x >= segment_start && x <= segment_start + 1.0);
                        let t = bspline_path
                            .t_at_length(distance, $accuracy)
                            .unwrap()
                            .unwrap();
                        let x = bspline_path.eval(t).unwrap()[0];
                        assert!(x >= segment_start && x <= segment_start + 1.0);
                    }
                }
            }
        }
    };
}

path_boundary_test!(path_boundaries_preserve_ownership_f64, f64, 1e-12);
path_boundary_test!(path_boundaries_preserve_ownership_f32, f32, 1e-5);

macro_rules! translated_box_test {
    ($name:ident, $scalar:ty, $origins:expr) => {
        #[test]
        fn $name() {
            for origin in $origins {
                for scale in [0.25 as $scalar, 1.0, 4.0, 16.0] {
                    // Quadratic B(t) = origin + 2*scale*t*(1-t).
                    // Its exact maximum at t=1/2 remains representable here.
                    let points = [
                        ArrayPoint::new([origin]),
                        ArrayPoint::new([origin + scale]),
                        ArrayPoint::new([origin]),
                    ];
                    let curve = Bezier::new(points).unwrap();
                    let spline = BSpline::<ArrayPoint<$scalar, 1>, 6, 3, 2>::new(
                        [0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
                        points,
                    )
                    .unwrap();
                    let expected = [(origin, origin + scale * 0.5)];
                    assert_eq!(curve.bounding_box().into_axes(), expected);
                    assert_eq!(spline.bounding_box().into_axes(), expected);
                    let mut path = BSplinePath::<ArrayPoint<$scalar, 1>, 6, 3, 2, 1>::new();
                    assert!(path.push(spline));
                    assert_eq!(path.bounding_box().unwrap().into_axes(), expected);
                    assert_eq!(curve.control_box()[0].1, origin + scale);
                }
                // Translation invariance also covers two non-dyadic extrema.
                let base = Bezier::new([
                    ArrayPoint::new([0.0 as $scalar]),
                    ArrayPoint::new([3.0]),
                    ArrayPoint::new([-2.0]),
                    ArrayPoint::new([0.0]),
                ])
                .unwrap()
                .bounding_box();
                let translated = Bezier::new([
                    ArrayPoint::new([origin]),
                    ArrayPoint::new([origin + 3.0]),
                    ArrayPoint::new([origin - 2.0]),
                    ArrayPoint::new([origin]),
                ])
                .unwrap()
                .bounding_box();
                assert_eq!(translated[0], (origin + base[0].0, origin + base[0].1));
            }
        }
    };
}

translated_box_test!(
    translated_boxes_remain_tight_f64,
    f64,
    [0.0_f64, 1e15, -1e15]
);
translated_box_test!(translated_boxes_remain_tight_f32, f32, [0.0_f32, 1e6, -1e6]);

macro_rules! partial_line_test {
    ($name:ident, $scalar:ty, $large:expr, $accuracy:expr, $relative:expr) => {
        #[test]
        fn $name() {
            type P = ArrayPoint<$scalar, 2>;
            let large: $scalar = $large;
            let line = LineSegment::new(P::new([-large, -large]), P::new([large, large])).unwrap();
            assert_eq!(
                line.arclen($accuracy),
                Err(ArcLengthError::NonFiniteComputation)
            );
            for (start, end) in [(0.0, 0.25), (0.25, 0.5), (0.5, 0.75), (0.75, 1.0)] {
                let expected = large * ((2.0 as $scalar).sqrt() * 0.5);
                for curve in [line, line.reverse()] {
                    assert_approx_equal(
                        curve.arclen_between(start, end, $accuracy).unwrap() / large,
                        expected / large,
                        $relative,
                        $relative,
                    );
                }
            }
            for t in [0.0, 0.5, 1.0] {
                assert_eq!(line.arclen_between(t, t, $accuracy), Ok(0.0));
            }
            let one_axis =
                LineSegment::new(ArrayPoint::new([-large]), ArrayPoint::new([large])).unwrap();
            assert_approx_equal(
                one_axis.arclen_between(0.25, 0.75, $accuracy).unwrap() / large,
                1.0,
                $relative,
                $relative,
            );
            for t in [0.0, 0.25, 0.5, 0.75, 1.0] {
                assert_approx_equal(
                    one_axis.eval(t)[0] / large,
                    2.0 * t - 1.0,
                    $relative,
                    $relative,
                );
            }
            // A subnormal parameter interval can describe a normal displacement.
            let fraction = <$scalar>::from_bits(1);
            assert_approx_equal(
                one_axis.arclen_between(0.0, fraction, $accuracy).unwrap(),
                large * (2.0 * fraction),
                0.0,
                $relative,
            );
            // Here the full displacement is finite but its norm overflows.
            let diagonal = LineSegment::new(P::new([0.0, 0.0]), P::new([large, large])).unwrap();
            let half = diagonal.arclen_between(0.0, 0.5, $accuracy).unwrap();
            assert_approx_equal(
                half / large,
                (2.0 as $scalar).sqrt() * 0.5,
                $relative,
                $relative,
            );
            assert_eq!(
                line.arclen_between(0.8, 0.2, $accuracy),
                Err(ArcLengthError::InvalidParameter)
            );
        }
    };
}

partial_line_test!(
    partial_line_lengths_avoid_overflow_f64,
    f64,
    1.4e308,
    1e290,
    1e-14
);
partial_line_test!(
    partial_line_lengths_avoid_overflow_f32,
    f32,
    3e38,
    1e30,
    1e-6
);

#[test]
fn recentered_boxes_preserve_endpoints_of_mixed_coordinate_scales() {
    for points in [
        [
            ArrayPoint::new([1e15_f64]),
            ArrayPoint::new([0.0]),
            ArrayPoint::new([-0.01]),
        ],
        [
            ArrayPoint::new([-1e15]),
            ArrayPoint::new([0.0]),
            ArrayPoint::new([0.01]),
        ],
    ] {
        let curve = Bezier::new(points).unwrap();
        let bounds = curve.bounding_box();
        assert!(bounds.contains(&points[0]));
        assert!(bounds.contains(&points[2]));
        assert_eq!(bounds, curve.reverse().bounding_box());
    }
}

#[test]
fn length_lookup_skips_zero_length_segments_without_losing_boundaries() {
    type P = ArrayPoint<f64, 1>;
    let mut path = BezierPath::<P, 6>::new();
    let mut bspline_path = BSplinePath::<P, 4, 2, 1, 6>::new();
    for (i, length) in [0.0, 1.0, 0.0, 0.0, 2.0, 0.0].into_iter().enumerate() {
        let start = P::new([i as f64 * 100.0]);
        let end = P::new([i as f64 * 100.0 + length]);
        assert!(path.push_line(LineSegment::new(start, end).unwrap()));
        assert!(bspline_path.push(BSpline::new([0.0, 0.0, 1.0, 1.0], [start, end]).unwrap()));
    }
    for (s, expected) in [(0.0, 0.0), (1.0, 400.0), (3.0, 500.0)] {
        let point = Some(P::new([expected]));
        assert_eq!(path.point_at_length(s, 1e-12).unwrap(), point);
        assert_eq!(bspline_path.point_at_length(s, 1e-12).unwrap(), point);
        assert_eq!(
            path.eval(path.t_at_length(s, 1e-12).unwrap().unwrap()),
            point
        );
        assert_eq!(
            bspline_path.eval(bspline_path.t_at_length(s, 1e-12).unwrap().unwrap()),
            point
        );
    }
    for i in 0..6 {
        let start = i as f64 / 6.0;
        let end = (i + 1) as f64 / 6.0;
        let expected = [0.0, 1.0, 0.0, 0.0, 2.0, 0.0][i];
        assert_eq!(path.arclen_between(start, end, 1e-12), Ok(expected));
        assert_eq!(bspline_path.arclen_between(start, end, 1e-12), Ok(expected));
    }
}
