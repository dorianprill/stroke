#![allow(incomplete_features)]
#![feature(generic_const_exprs)]

#[path = "../examples/support/oriented_reference.rs"]
mod reference;

use core::array;
use reference::{dot, identity, random, rotated};
use stroke::{ArrayPoint, BSpline, OrientedBoundingBoxError, Point};

fn check_faces<const DIM: usize, const K: usize, const C: usize, const D: usize>(
    curve: &BSpline<ArrayPoint<f64, DIM>, K, C, D>,
    budget: usize,
) where
    [(); <ArrayPoint<f64, DIM> as Point>::DIM]: Sized,
    [(); D + 1]: Sized,
    [(); C + 64]: Sized,
{
    let fit = curve
        .bounding_box_oriented_approx_with_budget(budget)
        .unwrap();
    assert!(fit.evaluations() <= budget);
    assert!(
        fit.log_volume()
            <= curve
                .bounding_box_oriented_approx_with_budget(1)
                .unwrap()
                .log_volume()
    );
    let center = fit.center();
    let half = fit.half_extents();
    for (range, piece) in curve.bezier_spans() {
        let controls = piece.into_control_points().map(ArrayPoint::into_array);
        for (axis, &extent) in fit.axes().iter().zip(&half) {
            let projected = controls.map(|p| (0..DIM).map(|d| axis[d] * (p[d] - center[d])).sum());
            let (lo, hi) = reference::range(projected);
            // Direct face checks use polynomial roots, without contains()'s allowance.
            assert!(
                lo >= -extent && hi <= extent,
                "{DIM}D degree {D}: {lo}..{hi}, half {extent}"
            );
        }
        for i in 0..=32 {
            assert!(fit.contains(&piece.eval(i as f64 / 32.0)));
            if i < 32 {
                let t = range.start * (1.0 - i as f64 / 32.0) + range.end * (i as f64 / 32.0);
                let actual = curve.eval(t).unwrap();
                assert!(fit.contains(&actual));
                let extracted = piece.eval(i as f64 / 32.0);
                for d in 0..DIM {
                    assert!((actual[d] - extracted[d]).abs() < 1e-10);
                }
            }
        }
    }
    for (i, axis) in fit.axes().iter().enumerate() {
        for (j, other) in fit.axes().iter().enumerate() {
            assert!(
                (dot(axis, other) - if i == j { 1.0 } else { 0.0 }).abs()
                    < 64.0 * DIM as f64 * f64::EPSILON
            );
        }
    }
}

#[test]
fn independent_faces_for_clamped_unclamped_and_repeated_knots_in_2d_3d_4d() {
    let mut state = 0x000b_511e;
    for _ in 0..8 {
        let controls = array::from_fn::<_, 7, _>(|_| {
            ArrayPoint::from(array::from_fn::<_, 2, _>(|_| random(&mut state)))
        });
        check_faces(
            &BSpline::<_, 11, 7, 3>::new_clamped_uniform(controls).unwrap(),
            256,
        );
        let controls = array::from_fn::<_, 7, _>(|_| {
            ArrayPoint::from(array::from_fn::<_, 3, _>(|_| random(&mut state)))
        });
        check_faces(
            &BSpline::<_, 11, 7, 3>::new(array::from_fn(|i| i as f64), controls).unwrap(),
            256,
        );
        let controls = array::from_fn::<_, 9, _>(|_| {
            ArrayPoint::from(array::from_fn::<_, 4, _>(|_| random(&mut state)))
        });
        check_faces(
            &BSpline::<_, 14, 9, 4>::new(
                [
                    0.0, 0.0, 0.0, 0.0, 0.0, 0.2, 0.5, 0.5, 0.8, 1.0, 1.0, 1.0, 1.0, 1.0,
                ],
                controls,
            )
            .unwrap(),
            256,
        );
    }
}

#[test]
fn discontinuities_preserve_both_limits_and_share_one_frame() {
    let curve = BSpline::<_, 12, 9, 2>::new(
        [0.0, 0.0, 0.0, 0.2, 0.2, 0.2, 0.8, 0.8, 0.8, 1.0, 1.0, 1.0],
        [
            [-30.0_f64, 10.0, 2.0],
            [-4.0, 3.0, -2.0],
            [20.0, -5.0, 1.0],
            [100.0, 50.0, 3.0],
            [103.0, 47.0, -4.0],
            [105.0, 60.0, 5.0],
            [-60.0, -40.0, -5.0],
            [-55.0, -30.0, 2.0],
            [-40.0, -25.0, 0.0],
        ]
        .map(ArrayPoint::from),
    )
    .unwrap();
    check_faces(&curve, 256);
    let fit = curve.bounding_box_oriented_approx_with_budget(256).unwrap();
    assert!(fit.diagnostics().unwrap().representative_points >= 6);
    for (_, piece) in curve.bezier_spans() {
        assert!(fit.contains(&piece.eval(0.0)));
        assert!(fit.contains(&piece.eval(1.0)));
    }
}

fn known_box<const DIM: usize, const C: usize>()
where
    [(); <ArrayPoint<f64, DIM> as Point>::DIM]: Sized,
    [(); C + 1]: Sized,
    [(); C + 64]: Sized,
{
    assert_eq!(C, 1 << DIM);
    let widths = [4.0, 0.7, 0.2, 0.09];
    let mut axes = identity::<DIM>();
    for i in 0..DIM {
        for j in i + 1..DIM {
            axes = rotated(axes, i, j, 0.19 * (i + j + 1) as f64);
        }
    }
    let controls = array::from_fn::<_, C, _>(|mask| {
        let local = array::from_fn(|i| {
            if mask & (1 << i) == 0 {
                -0.5 * widths[i]
            } else {
                0.5 * widths[i]
            }
        });
        ArrayPoint::from(axes.map(|axis| dot(&axis, &local)))
    });
    // Degree-zero spans are isolated points. Their convex hull is exactly this
    // rotated cuboid, whose volume lower-bounds every enclosing box.
    let curve = BSpline::<_, { C + 1 }, C, 0>::new(array::from_fn(|i| i as f64), controls).unwrap();
    check_faces::<DIM, { C + 1 }, C, 0>(&curve, 1024);
    let fit = curve
        .bounding_box_oriented_approx_with_budget(1024)
        .unwrap();
    let exact: f64 = widths[..DIM].iter().product();
    let ratio = fit.log_volume().exp() / exact;
    assert!(
        (1.0..1.01).contains(&ratio),
        "{DIM}D known-box ratio {ratio}"
    );
}

#[test]
fn known_minimum_volume_boxes_in_2d_3d_and_4d() {
    known_box::<2, 4>();
    known_box::<3, 8>();
    known_box::<4, 16>();
}

#[test]
fn single_span_quality_matches_independent_polynomial_reference() {
    let frame = rotated(identity(), 0, 1, 0.63);
    let controls =
        [[0.0, 0.0], [1.0, 3.0], [3.0, 3.0], [4.0, 0.0]].map(|p| frame.map(|a| dot(&a, &p)));
    let curve = BSpline::<_, 8, 4, 3>::new_clamped_uniform(controls.map(ArrayPoint::from)).unwrap();
    let fit = curve.bounding_box_oriented_approx_with_budget(256).unwrap();
    let numerical =
        reference::reference(&controls).min(reference::reference_at(&controls, *fit.axes()));
    assert!((fit.log_volume() - numerical).exp() < 1.0001);
}

#[test]
fn more_than_64_spans_are_not_truncated() {
    let mut controls = [ArrayPoint::from([0.0_f64, 0.0]); 80];
    for (i, p) in controls.iter_mut().enumerate() {
        *p = ArrayPoint::from([i as f64, (i as f64 * 0.3).sin()]);
    }
    controls[79] = ArrayPoint::from([-100.0, 200.0]);
    let curve = BSpline::<_, 81, 80, 0>::new(array::from_fn(|i| i as f64), controls).unwrap();
    check_faces(&curve, 128);
    let info = curve
        .bounding_box_oriented_approx_with_budget(128)
        .unwrap()
        .diagnostics()
        .unwrap();
    assert_eq!(info.representative_points, 80);
    assert!(!info.representative_limit_reached);
    assert_eq!(info.representative_error, 0.0);

    let line = BSpline::<_, 82, 80, 1>::new_clamped_uniform(controls).unwrap();
    check_faces(&line, 128);
    assert!(
        line.bounding_box_oriented_approx_with_budget(128)
            .unwrap()
            .contains(&controls[79])
    );
}

#[test]
fn higher_degree_in_six_dimensions() {
    let mut state = 0x0006_db51;
    let controls = array::from_fn::<_, 10, _>(|_| {
        ArrayPoint::from(array::from_fn::<_, 6, _>(|_| random(&mut state)))
    });
    let curve = BSpline::<_, 19, 10, 8>::new_clamped_uniform(controls).unwrap();
    let fit = curve.bounding_box_oriented_approx_with_budget(128).unwrap();
    for (_, span) in curve.bezier_spans() {
        for i in 0..=256 {
            assert!(fit.contains(&span.eval(i as f64 / 256.0)));
        }
    }
}

#[test]
fn deterministic_prefix_budgets_include_preprocessing_and_jumps() {
    let curve = BSpline::<_, 10, 7, 2>::new(
        [0.0, 0.0, 0.0, 0.3, 0.3, 0.3, 0.7, 1.0, 1.0, 1.0],
        [
            [0.0_f64; 3],
            [1.0, 3.0, 2.0],
            [4.0, -2.0, 1.0],
            [7.0, 1.0, 4.0],
            [3.0, -1.0, 3.0],
            [9.0, 4.0, -2.0],
            [10.0, 2.0, 1.0],
        ]
        .map(ArrayPoint::from),
    )
    .unwrap();
    let mut best = f64::INFINITY;
    let mut preprocessing = None;
    for budget in [1, 2, 3, 4, 5, 8, 17, 63, 64, 65, 256, 1024] {
        let fit = curve
            .bounding_box_oriented_approx_with_budget(budget)
            .unwrap();
        assert_eq!(
            fit,
            curve
                .bounding_box_oriented_approx_with_budget(budget)
                .unwrap()
        );
        assert!(fit.evaluations() <= budget);
        assert!(fit.log_volume() <= best);
        best = fit.log_volume();
        if budget > 1 {
            let info = fit.diagnostics().unwrap();
            assert!(info.model_log_volume_gap.is_finite() && info.model_log_volume_gap >= 0.0);
            let fixed = (
                info.representative_points,
                info.point_set_evaluations,
                info.representative_error,
            );
            if let Some(previous) = preprocessing {
                assert_eq!(previous, fixed);
            }
            preprocessing = Some(fixed);
        } else {
            assert!(fit.diagnostics().is_none());
        }
    }
}

#[test]
fn constants_one_dimension_and_request_errors() {
    let constant =
        BSpline::<_, 7, 4, 2>::new_clamped_uniform([ArrayPoint::from([2.0_f64, -3.0]); 4]).unwrap();
    let fit = constant
        .bounding_box_oriented_approx_with_budget(64)
        .unwrap();
    assert_eq!(fit.evaluations(), 0);
    assert_eq!(fit.half_extents(), [0.0; 2]);
    assert!(fit.diagnostics().is_none());
    let line = BSpline::<_, 5, 3, 1>::new_clamped_uniform(
        [[0.0_f64], [3.0], [-2.0]].map(ArrayPoint::from),
    )
    .unwrap();
    let fit = line.bounding_box_oriented_approx_with_budget(64).unwrap();
    assert_eq!(fit.evaluations(), 1);
    assert!(fit.diagnostics().is_none());
    assert_eq!(
        line.bounding_box_oriented_approx_with_budget(0),
        Err(OrientedBoundingBoxError::ZeroBudget)
    );
    let zero =
        BSpline::<_, 4, 2, 1>::new_clamped_uniform([ArrayPoint::<f64, 0>::from([]); 2]).unwrap();
    assert_eq!(
        zero.bounding_box_oriented_approx_with_budget(64),
        Err(OrientedBoundingBoxError::ZeroDimension)
    );
    let sub = f64::from_bits(1);
    let extreme = BSpline::<_, 4, 2, 1>::new_clamped_uniform(
        [[1e300_f64, 0.0], [1e300, sub]].map(ArrayPoint::from),
    )
    .unwrap();
    assert_eq!(
        extreme.bounding_box_oriented_approx_with_budget(64),
        Err(OrientedBoundingBoxError::NonFiniteComputation)
    );
}

#[test]
fn coordinate_and_knot_scales_f32_and_f64() {
    let base = [
        [0.0_f64; 3],
        [1.0, 3.0, 2.0],
        [4.0, -2.0, 1.0],
        [5.0, 1.0, 4.0],
        [3.0, -1.0, 2.0],
    ];
    for scale in [1e-200, 1.0, 1e200] {
        for knots in [
            [-1e300, -1e300, -1e300, -1e299, 1e299, 1e300, 1e300, 1e300],
            [0.0, 0.0, 0.0, 1e-300, 2e-300, 3e-300, 3e-300, 3e-300],
        ] {
            let curve = BSpline::<_, 8, 5, 2>::new(
                knots,
                base.map(|p| ArrayPoint::from(p.map(|x| x * scale))),
            )
            .unwrap();
            let fit = curve.bounding_box_oriented_approx_with_budget(128).unwrap();
            for (_, span) in curve.bezier_spans() {
                for i in 0..=64 {
                    assert!(fit.contains(&span.eval(i as f64 / 64.0)));
                }
            }
        }
    }
    for (scale, offset) in [(1e-20_f32, 0.0_f32), (1.0, 1e6), (1e20, 0.0)] {
        let curve = BSpline::<_, 8, 5, 2>::new_clamped_uniform(
            base.map(|p| ArrayPoint::from(p.map(|x| x as f32 * scale + offset))),
        )
        .unwrap();
        let fit = curve.bounding_box_oriented_approx_with_budget(128).unwrap();
        for (_, span) in curve.bezier_spans() {
            for i in 0..=64 {
                assert!(fit.contains(&span.eval(i as f32 / 64.0)));
            }
        }
    }
}

#[cfg(feature = "nalgebra")]
#[test]
fn nalgebra_uses_the_shared_search() {
    let curve = BSpline::<_, 7, 4, 2>::new_clamped_uniform([
        nalgebra::SVector::<f64, 3>::new(0.0, 0.0, 0.0),
        nalgebra::SVector::<f64, 3>::new(1.0, 3.0, 2.0),
        nalgebra::SVector::<f64, 3>::new(4.0, -2.0, 1.0),
        nalgebra::SVector::<f64, 3>::new(5.0, 1.0, 4.0),
    ])
    .unwrap();
    let fit = curve.bounding_box_oriented_approx_with_budget(256).unwrap();
    for i in 0..=128 {
        assert!(fit.contains(&curve.eval(i as f64 / 128.0).unwrap()));
    }
}
