//! Cross-representation regression tests for length, extraction and bounds.

use super::*;
use crate::bspline::BSplineError;
use crate::test_support::{assert_approx_equal, assert_point_approx_equal};

type P = PointN<f64, 2>;

fn ph_cubic() -> CubicBezier<P> {
    // C(t) = (t - t^3/3, t^2), whose speed is exactly 1 + t^2.
    CubicBezier::new(
        P::new([0.0, 0.0]),
        P::new([1.0 / 3.0, 0.0]),
        P::new([2.0 / 3.0, 1.0 / 3.0]),
        P::new([2.0 / 3.0, 1.0]),
    )
    .unwrap()
}

#[test]
fn stable_norm_preserves_representable_extreme_magnitudes() {
    for scale in [1e-200_f64, 1e200] {
        let p = P::new([3.0 * scale, 4.0 * scale]);
        assert_approx_equal(p.norm() / scale, 5.0, 1e-14, 1e-14);
    }
    for scale in [1e-30_f32, 1e30] {
        let p = PointN::new([3.0 * scale, 4.0 * scale]);
        assert_approx_equal(p.norm() / scale, 5.0, 1e-5, 1e-5);
    }
    assert!(P::new([0.0, f64::NAN]).norm().is_nan());
    assert_eq!(P::new([0.0, f64::INFINITY]).norm(), f64::INFINITY);
    assert_eq!(P::new([0.0, 0.0]).norm(), 0.0);
}

#[test]
fn line_length_matches_every_representation() {
    let points = [P::new([0.0, 0.0]), P::new([1.5, 2.0]), P::new([3.0, 4.0])];
    let line = LineSegment::new(points[0], points[2]).unwrap();
    let quad = QuadraticBezier::new(points[0], points[1], points[2]).unwrap();
    let generic = Bezier::new(points).unwrap();
    let spline = BSpline::<P, 6, 3, 2>::new([0.0, 0.0, 0.0, 1.0, 1.0, 1.0], points).unwrap();
    for length in [
        line.arclen(1e-12).unwrap(),
        quad.arclen(1e-12).unwrap(),
        generic.arclen(1e-12).unwrap(),
        spline.arclen(1e-12).unwrap(),
        bezier_segment::BezierSegment::from(quad)
            .arclen(1e-12)
            .unwrap(),
    ] {
        assert_approx_equal(length, 5.0, 1e-12, 1e-12);
    }
}

#[test]
fn quadrature_matches_analytic_parabola_length() {
    let points = [P::new([0.0, 0.0]), P::new([0.5, 0.0]), P::new([1.0, 1.0])];
    let quad = QuadraticBezier::new(points[0], points[1], points[2]).unwrap();
    let generic = Bezier::new(points).unwrap();
    let expected = 5.0_f64.sqrt() * 0.5 + 2.0_f64.asinh() * 0.25;
    assert_approx_equal(quad.arclen(1e-12).unwrap(), expected, 1e-12, 1e-12);
    assert_approx_equal(generic.arclen(1e-12).unwrap(), expected, 1e-12, 1e-12);
    assert!(quad.arclen_polyline(4) < expected);
}

#[test]
fn analytic_cubic_length_inverse_and_partial_integral() {
    let curve = ph_cubic();
    assert_approx_equal(curve.arclen(1e-12).unwrap(), 4.0 / 3.0, 1e-12, 1e-12);
    for t in [0.0_f64, 0.07, 0.31, 0.89, 1.0] {
        let length = t + t.powi(3) / 3.0;
        assert_approx_equal(
            curve.arclen_between(0.0, t, 1e-12).unwrap(),
            length,
            1e-12,
            1e-12,
        );
        assert_approx_equal(curve.t_at_length(length, 1e-12).unwrap(), t, 2e-12, 2e-12);
        assert_point_approx_equal(
            curve.point_at_length(length, 1e-12).unwrap(),
            curve.eval_casteljau(t),
            2e-12,
            2e-12,
        );
    }
    assert_eq!(curve.t_at_length(-1.0, 1e-12), Ok(0.0));
    assert_eq!(curve.t_at_length(10.0, 1e-12), Ok(1.0));
}

#[test]
fn backtracking_and_non_dyadic_cusps_have_full_length() {
    // B(t) = 2t - 3t^2 has a turning point at t = 1/3.
    let curve = Bezier::new([
        PointN::new([0.0_f64]),
        PointN::new([1.0]),
        PointN::new([-1.0]),
    ])
    .unwrap();
    assert_approx_equal(curve.arclen(1e-12).unwrap(), 5.0 / 3.0, 1e-12, 1e-12);
    assert_approx_equal(
        curve.t_at_length(1.0 / 3.0, 1e-12).unwrap(),
        1.0 / 3.0,
        1e-6,
        0.0,
    );

    let cubic = Bezier::new([
        PointN::new([0.0_f64]),
        PointN::new([2.0]),
        PointN::new([-1.0]),
        PointN::new([0.0]),
    ])
    .unwrap();
    let a = (5.0 - 7.0_f64.sqrt()) / 9.0;
    let b = (5.0 + 7.0_f64.sqrt()) / 9.0;
    let expected = 2.0 * (cubic.eval(a)[0] - cubic.eval(b)[0]);
    assert_approx_equal(cubic.arclen(1e-12).unwrap(), expected, 1e-12, 1e-12);
}

#[test]
fn length_is_invariant_under_reverse_split_and_scale() {
    let curve = ph_cubic();
    let (left, right) = curve.split(0.37);
    let length = curve.arclen(1e-12).unwrap();
    assert_approx_equal(curve.reverse().arclen(1e-12).unwrap(), length, 1e-12, 1e-12);
    assert_approx_equal(
        left.arclen(5e-13).unwrap() + right.arclen(5e-13).unwrap(),
        length,
        1e-12,
        1e-12,
    );
    for scale in [1e-200, 1e200] {
        let scaled = Bezier::new([
            P::new([0.0, 0.0]),
            P::new([scale / 3.0, 0.0]),
            P::new([2.0 * scale / 3.0, scale / 3.0]),
            P::new([2.0 * scale / 3.0, scale]),
        ])
        .unwrap();
        assert_approx_equal(
            scaled.arclen(scale * 1e-12).unwrap() / scale,
            length,
            1e-12,
            1e-12,
        );
    }
}

#[test]
fn f32_quadrature_and_inverse_work_at_distinct_scales() {
    for scale in [1e-25_f32, 1.0, 1e25] {
        let curve = Bezier::new([
            PointN::new([0.0, 0.0]),
            PointN::new([scale / 3.0, 0.0]),
            PointN::new([2.0 * scale / 3.0, scale / 3.0]),
            PointN::new([2.0 * scale / 3.0, scale]),
        ])
        .unwrap();
        let length = curve.arclen(scale * 1e-5).unwrap();
        assert_approx_equal(length / scale, 4.0 / 3.0, 1e-5, 1e-5);
        let t = curve.t_at_length(length * 0.5, scale * 1e-5).unwrap();
        assert_approx_equal(
            curve.arclen_between(0.0, t, scale * 1e-5).unwrap() / scale,
            2.0 / 3.0,
            2e-5,
            2e-5,
        );
    }
}

#[test]
fn constant_curves_and_degree_zero_splines_work() {
    let curve = Bezier::new([P::new([3.0, 4.0])]).unwrap();
    assert_eq!(curve.arclen(1e-12), Ok(0.0));
    assert_eq!(curve.t_at_length(1.0, 1e-12), Ok(0.0));
    assert_eq!(curve.bounding_box().into_axes(), [(3.0, 3.0), (4.0, 4.0)]);
    assert_eq!(curve.tangent(0.5), P::new([0.0, 0.0]));
    let spline =
        BSpline::<P, 3, 2, 0>::new([2.0, 3.0, 4.0], [P::new([0.0, 0.0]), P::new([10.0, 20.0])])
            .unwrap();
    assert_eq!(spline.arclen(1e-12), Ok(0.0));
    assert_eq!(spline.arclen_polyline(8), 0.0);
    assert_eq!(spline.bezier_spans().count(), 2);
    assert_eq!(
        spline.bounding_box().into_axes(),
        [(0.0, 10.0), (0.0, 20.0)]
    );
}

#[test]
fn invalid_accuracy_parameters_and_overflow_are_reported() {
    let curve = ph_cubic();
    for accuracy in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert_eq!(curve.arclen(accuracy), Err(ArcLengthError::InvalidAccuracy));
    }
    for s in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert_eq!(
            curve.t_at_length(s, 1e-12),
            Err(ArcLengthError::InvalidParameter)
        );
    }
    for (start, end) in [(-0.1, 0.5), (0.8, 0.3), (0.0, 1.1), (f64::NAN, 1.0)] {
        assert_eq!(
            curve.arclen_between(start, end, 1e-12),
            Err(ArcLengthError::InvalidParameter)
        );
    }
    let huge = LineSegment::new(P::new([-f64::MAX, 0.0]), P::new([f64::MAX, 0.0])).unwrap();
    assert_eq!(
        huge.arclen(1e-12),
        Err(ArcLengthError::NonFiniteComputation)
    );
}

#[test]
fn integration_reports_exhausted_convergence_budget() {
    let result = crate::arc_length::integrate(0.0_f64, 1.0, 1e-16, |t| {
        Ok(1.0 / (t - 0.123_456_789).abs().sqrt())
    });
    assert_eq!(result, Err(ArcLengthError::FailedToConverge));
}

#[test]
fn tight_boxes_exclude_off_curve_control_points() {
    let points = [P::new([0.0, 0.0]), P::new([1.0, 2.0]), P::new([2.0, 0.0])];
    let quad = QuadraticBezier::new(points[0], points[1], points[2]).unwrap();
    let generic = Bezier::new(points).unwrap();
    let spline = BSpline::<P, 6, 3, 2>::new([0.0, 0.0, 0.0, 1.0, 1.0, 1.0], points).unwrap();
    for bounds in [
        quad.bounding_box(),
        generic.bounding_box(),
        spline.bounding_box(),
    ] {
        assert_approx_equal(bounds[1].1, 1.0, 1e-13, 0.0);
        assert!(!bounds.contains(&points[1]));
        assert!(bounds.contains(&P::new([1.0, 1.0])));
    }
    for bounds in [
        quad.control_box(),
        generic.control_box(),
        spline.control_box(),
    ] {
        assert_eq!(bounds.into_axes(), [(0.0, 2.0), (0.0, 2.0)]);
    }
}

#[test]
fn tight_boxes_and_roots_remain_scale_independent() {
    for scale in [1e-150_f64, 1.0, 1e150] {
        let curve = Bezier::new([
            PointN::new([0.0]),
            PointN::new([4.0 * scale]),
            PointN::new([-4.0 * scale]),
            PointN::new([0.0]),
        ])
        .unwrap();
        let bounds = curve.bounding_box();
        let expected = 2.0 / 3.0_f64.sqrt();
        assert_approx_equal(bounds[0].0 / scale, -expected, 1e-12, 1e-12);
        assert_approx_equal(bounds[0].1 / scale, expected, 1e-12, 1e-12);
        let roots = curve.derivative_roots(0);
        assert_eq!(roots.len(), 2);
        assert_approx_equal(roots[0], 0.5 - 3.0_f64.sqrt() / 6.0, 1e-12, 0.0);
        for step in 0..=1000 {
            let p = curve.eval(step as f64 / 1000.0);
            assert!(p[0] >= bounds[0].0 - scale * 1e-14);
            assert!(p[0] <= bounds[0].1 + scale * 1e-14);
        }
    }
}

#[test]
fn derivative_roots_include_tangency_and_flat_axes() {
    // Derivative is 3(t - 1/2)^2, with an even-multiplicity zero.
    let curve = Bezier::new([
        P::new([0.0, 3.0]),
        P::new([0.25, 3.0]),
        P::new([0.0, 3.0]),
        P::new([0.25, 3.0]),
    ])
    .unwrap();
    let roots = curve.derivative_roots(0);
    assert_eq!(roots.len(), 1);
    assert_approx_equal(roots[0], 0.5, 1e-12, 0.0);
    assert!(curve.derivative_roots(1).is_empty());
    assert_eq!(curve.bounding_box().into_axes(), [(0.0, 0.25), (3.0, 3.0)]);
}

fn check_extraction<const K: usize, const C: usize, const D: usize>(spline: &BSpline<P, K, C, D>)
where
    [(); D + 1]: Sized,
{
    let bounds = spline.bounding_box();
    for (range, curve) in spline.bezier_spans() {
        for step in 0..=64 {
            let local = step as f64 / 64.0;
            let t = range.start * (1.0 - local) + range.end * local;
            assert_point_approx_equal(spline.eval(t).unwrap(), curve.eval(local), 2e-12, 2e-12);
            let point = curve.eval(local);
            for axis in 0..P::DIM {
                assert!(point[axis] >= bounds[axis].0 - 1e-12);
                assert!(point[axis] <= bounds[axis].1 + 1e-12);
            }
        }
    }
}

#[test]
fn extraction_matches_clamped_unclamped_and_partly_clamped_splines() {
    let points = [
        P::new([0.0, 0.0]),
        P::new([1.0, 3.0]),
        P::new([2.0, -2.0]),
        P::new([3.0, 1.0]),
    ];
    for knots in [
        [0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
        [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0],
        [0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 4.0],
        [0.0, 1.0, 2.0, 3.0, 4.0, 4.0, 4.0, 4.0],
    ] {
        let spline = BSpline::<P, 8, 4, 3>::new(knots, points).unwrap();
        assert_eq!(spline.knot_domain(), (knots[3], knots[4]));
        check_extraction(&spline);
    }
    let spline = BSpline::<P, 7, 4, 2>::new([0.0, 0.0, 0.0, 0.3, 1.0, 1.0, 1.0], points).unwrap();
    check_extraction(&spline);
    assert_eq!(spline.bezier_spans().count(), 2);
}

#[test]
fn repeated_knots_and_nonuniform_multispans_preserve_geometry() {
    let points = [
        P::new([0.0, 0.0]),
        P::new([0.5, 2.0]),
        P::new([1.0, -3.0]),
        P::new([2.0, 4.0]),
        P::new([3.0, -2.0]),
        P::new([4.0, 1.0]),
    ];
    let spline =
        BSpline::<P, 10, 6, 3>::new([0.0, 0.0, 0.0, 0.0, 0.2, 0.2, 1.0, 1.0, 1.0, 1.0], points)
            .unwrap();
    check_extraction(&spline);
    let sum: f64 = spline
        .bezier_spans()
        .map(|(_, curve)| curve.arclen(5e-13).unwrap())
        .sum();
    assert_approx_equal(spline.arclen(1e-12).unwrap(), sum, 1e-12, 1e-12);
    assert_approx_equal(spline.reverse().arclen(1e-12).unwrap(), sum, 1e-12, 1e-12);
    let target = sum * 0.4;
    let t = spline.t_at_length(target, 1e-11).unwrap();
    let mut prefix = 0.0;
    for (range, curve) in spline.bezier_spans() {
        if t >= range.end {
            prefix += curve.arclen(1e-12).unwrap();
        } else if t > range.start {
            prefix += curve
                .arclen_between(0.0, (t - range.start) / (range.end - range.start), 1e-12)
                .unwrap();
        }
    }
    assert_approx_equal(prefix, target, 2e-11, 1e-12);
}

#[test]
fn knot_reparameterization_does_not_change_length_or_bounds() {
    let points = [P::new([0.0, 0.0]), P::new([0.5, 0.0]), P::new([1.0, 1.0])];
    let reference = Bezier::new(points).unwrap();
    for (a, b) in [(0.0, 1e-200), (10.0, 10.000_000_001), (-f64::MAX, f64::MAX)] {
        let spline = BSpline::<P, 6, 3, 2>::new([a, a, a, b, b, b], points).unwrap();
        assert_approx_equal(
            spline.arclen(1e-12).unwrap(),
            reference.arclen(1e-12).unwrap(),
            1e-12,
            1e-12,
        );
        assert_point_approx_equal(
            spline.eval(a * 0.5 + b * 0.5).unwrap(),
            reference.eval(0.5),
            1e-12,
            1e-12,
        );
        assert_eq!(spline.bounding_box(), reference.bounding_box());
    }
}

#[test]
fn discontinuous_spans_include_both_limits_without_jump_lengths() {
    let points = [
        P::new([0.0, 0.0]),
        P::new([1.0, 2.0]),
        P::new([2.0, 0.0]),
        P::new([20.0, 0.0]),
        P::new([21.0, 2.0]),
        P::new([22.0, 0.0]),
    ];
    let spline =
        BSpline::<P, 9, 6, 2>::new([0.0, 0.0, 0.0, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0], points).unwrap();
    let mut spans = spline.bezier_spans();
    let (_, left) = spans.next().unwrap();
    let (_, right) = spans.next().unwrap();
    assert_eq!(left.end(), P::new([2.0, 0.0]));
    assert_eq!(right.start(), P::new([20.0, 0.0]));
    assert_eq!(spline.eval(0.5).unwrap(), right.start());
    let length = left.arclen(1e-12).unwrap();
    assert_approx_equal(spline.arclen(2e-12).unwrap(), 2.0 * length, 2e-12, 1e-12);
    assert_approx_equal(
        spline.arclen_polyline(16),
        left.arclen_polyline(16) + right.arclen_polyline(16),
        1e-12,
        1e-12,
    );
    let bounds = spline.bounding_box();
    assert_approx_equal(bounds[1].1, 1.0, 1e-12, 0.0);
    assert_eq!(bounds[0], (0.0, 22.0));
}

#[test]
fn disconnected_and_empty_paths_have_consistent_length_semantics() {
    let mut path = BezierPath::<P, 3>::new();
    assert_eq!(path.arclen(1e-12), Ok(0.0));
    assert_eq!(path.t_at_length(1.0, 1e-12), Ok(None));
    assert_eq!(path.point_at_length(1.0, 1e-12), Ok(None));
    assert!(path.bounding_box().is_none());
    path.push_line(LineSegment::new(P::new([0.0, 0.0]), P::new([1.0, 0.0])).unwrap());
    path.push_line(LineSegment::new(P::new([100.0, 0.0]), P::new([103.0, 0.0])).unwrap());
    assert_eq!(path.arclen(1e-12), Ok(4.0));
    assert_eq!(path.arclen_polyline(8), 4.0);
    assert_eq!(
        path.point_at_length(1.0, 1e-12).unwrap(),
        Some(P::new([100.0, 0.0]))
    );
    assert_point_approx_equal(
        path.point_at_length(2.5, 1e-12).unwrap().unwrap(),
        P::new([101.5, 0.0]),
        1e-12,
        1e-12,
    );
    assert_eq!(path.t_at_length(20.0, 1e-12), Ok(Some(1.0)));
}

#[test]
fn bspline_path_uses_tight_segment_boxes_and_additive_lengths() {
    let points = [P::new([0.0, 0.0]), P::new([1.0, 2.0]), P::new([2.0, 0.0])];
    let segment = BSpline::<P, 6, 3, 2>::new([2.0, 2.0, 2.0, 5.0, 5.0, 5.0], points).unwrap();
    let mut path = BSplinePath::<P, 6, 3, 2, 2>::new();
    assert_eq!(path.arclen(1e-12), Ok(0.0));
    assert_eq!(path.point_at_length(1.0, 1e-12), Ok(None));
    path.push(segment);
    path.push(segment.reverse());
    let length = segment.arclen(1e-12).unwrap();
    assert_approx_equal(path.arclen(2e-12).unwrap(), length * 2.0, 2e-12, 1e-12);
    assert_approx_equal(path.bounding_box().unwrap()[1].1, 1.0, 1e-12, 0.0);
    assert_eq!(path.control_box().unwrap()[1].1, 2.0);
    assert_approx_equal(
        path.t_at_length(length * 1.5, 1e-11).unwrap().unwrap(),
        0.75,
        1e-11,
        0.0,
    );
}

#[test]
fn boxes_support_union_axes_and_validated_construction() {
    let a = BoundingBox::new([0.0_f64, 1.0], [2.0, 3.0]).unwrap();
    let b = BoundingBox::new([-2.0, 2.0], [1.0, 4.0]).unwrap();
    assert_eq!(a.union(&b).into_axes(), [(-2.0, 2.0), (1.0, 4.0)]);
    assert_eq!(a.axis(1), (1.0, 3.0));
    assert!(BoundingBox::new([2.0], [1.0]).is_none());
    assert!(BoundingBox::new([f64::NAN], [1.0]).is_none());
}

#[test]
fn specialized_boxes_handle_large_finite_coordinate_ranges() {
    for scale in [1e-200_f64, 1e200, f64::MAX * 0.25] {
        let cubic = CubicBezier::new(
            P::new([0.0, 0.0]),
            P::new([scale, 2.0 * scale]),
            P::new([2.0 * scale, 2.0 * scale]),
            P::new([3.0 * scale, 0.0]),
        )
        .unwrap();
        assert_approx_equal(cubic.bounding_box()[1].1 / scale, 1.5, 1e-12, 1e-12);
        let quad = QuadraticBezier::new(
            P::new([0.0, 0.0]),
            P::new([scale, 2.0 * scale]),
            P::new([2.0 * scale, 0.0]),
        )
        .unwrap();
        assert_approx_equal(quad.bounding_box()[1].1 / scale, 1.0, 1e-12, 1e-12);
    }
}

#[test]
fn higher_degree_tight_bounds_and_lengths_match_independent_samples() {
    let points = [
        P::new([0.0, 0.0]),
        P::new([0.5, 4.0]),
        P::new([2.0, -3.0]),
        P::new([3.0, 4.0]),
        P::new([4.5, -2.0]),
        P::new([5.0, 1.0]),
    ];
    let curve = Bezier::new(points).unwrap();
    let spline = BSpline::<P, 12, 6, 5>::new(
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
        points,
    )
    .unwrap();
    assert_eq!(curve.bounding_box(), spline.bounding_box());
    assert_approx_equal(
        curve.arclen(1e-12).unwrap(),
        spline.arclen(1e-12).unwrap(),
        1e-12,
        1e-12,
    );
    let (left, right) = curve.split(0.13);
    assert_approx_equal(
        curve.arclen(1e-12).unwrap(),
        left.arclen(5e-13).unwrap() + right.arclen(5e-13).unwrap(),
        1e-12,
        1e-12,
    );
    let bounds = curve.bounding_box();
    let mut sampled_min = [f64::INFINITY; 2];
    let mut sampled_max = [f64::NEG_INFINITY; 2];
    for step in 0..=10000 {
        let point = curve.eval(step as f64 / 10000.0);
        for axis in 0..2 {
            sampled_min[axis] = sampled_min[axis].min(point[axis]);
            sampled_max[axis] = sampled_max[axis].max(point[axis]);
            assert!(point[axis] >= bounds[axis].0 - 1e-12);
            assert!(point[axis] <= bounds[axis].1 + 1e-12);
        }
    }
    for axis in 0..2 {
        assert!((bounds[axis].0 - sampled_min[axis]).abs() < 1e-6);
        assert!((bounds[axis].1 - sampled_max[axis]).abs() < 1e-6);
    }
}

#[test]
fn partial_length_apis_share_path_and_knot_parameter_semantics() {
    let line = LineSegment::new(P::new([0.0, 0.0]), P::new([3.0, 4.0])).unwrap();
    assert_approx_equal(
        line.arclen_between(0.2, 0.8, 1e-12).unwrap(),
        3.0,
        1e-12,
        1e-12,
    );
    let mut path = BezierPath::<P, 2>::new();
    path.push_line(line);
    path.push_line(line.reverse());
    assert_approx_equal(
        path.arclen_between(0.2, 0.7, 1e-12).unwrap(),
        5.0,
        1e-12,
        1e-12,
    );

    let points = [P::new([0.0, 0.0]), P::new([1.5, 2.0]), P::new([3.0, 4.0])];
    let spline = BSpline::<P, 6, 3, 2>::new([2.0, 2.0, 2.0, 5.0, 5.0, 5.0], points).unwrap();
    assert_approx_equal(
        spline.arclen_between(2.6, 4.4, 1e-12).unwrap(),
        3.0,
        1e-12,
        1e-12,
    );
    let mut path = BSplinePath::<P, 6, 3, 2, 2>::new();
    path.push(spline);
    path.push(spline.reverse());
    assert_approx_equal(
        path.arclen_between(0.2, 0.7, 1e-12).unwrap(),
        5.0,
        1e-12,
        1e-12,
    );
}

#[test]
fn bspline_construction_rejects_empty_domain_and_excess_multiplicity() {
    assert!(matches!(
        BSpline::<P, 4, 2, 1>::new([0.0; 4], [P::new([0.0, 0.0]); 2]),
        Err(BSplineError::EmptyKnotDomain)
    ));
    assert!(matches!(
        BSpline::<P, 6, 4, 1>::new([0.0, 0.0, 0.5, 0.5, 0.5, 1.0], [P::new([0.0, 0.0]); 4]),
        Err(BSplineError::ExcessKnotMultiplicity)
    ));
    let curve = BSpline::<P, 6, 3, 2>::default();
    assert!(curve.eval(0.5).is_ok());
    assert_eq!(curve.arclen(1e-12), Ok(0.0));
}

#[cfg(feature = "nalgebra")]
#[test]
fn nalgebra_points_share_length_and_bounds_apis() {
    use nalgebra::SVector;
    let curve = Bezier::new([
        SVector::<f64, 2>::new(0.0, 0.0),
        SVector::<f64, 2>::new(1.0, 2.0),
        SVector::<f64, 2>::new(2.0, 0.0),
    ])
    .unwrap();
    assert_approx_equal(curve.bounding_box()[1].1, 1.0, 1e-12, 0.0);
    assert_approx_equal(
        curve
            .t_at_length(curve.arclen(1e-12).unwrap() * 0.5, 1e-12)
            .unwrap(),
        0.5,
        1e-12,
        0.0,
    );
    assert!(
        curve
            .control_box()
            .contains(&SVector::<f64, 2>::new(1.0, 2.0))
    );
}
