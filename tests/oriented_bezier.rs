#![allow(incomplete_features)]
#![feature(generic_const_exprs)]

#[path = "../examples/support/oriented_reference.rs"]
mod reference;

use core::array;
use reference::{dot, identity, random, rotated};
use stroke::{
    ArrayPoint, Bezier, CubicBezier, OrientedBoundingBoxError, Point, PointIndex, QuadraticBezier,
};

fn check_specialized<P: PointIndex>(controls: [P; 4])
where
    P::Scalar: core::fmt::Debug,
    [(); P::DIM]: Sized,
{
    let generic_quadratic = Bezier::new([controls[0], controls[1], controls[2]]).unwrap();
    let generic_cubic = Bezier::new(controls).unwrap();
    let quadratic = QuadraticBezier::from(generic_quadratic);
    let cubic = CubicBezier::from(generic_cubic);
    assert_eq!(
        quadratic.bounding_box_oriented_approx(),
        generic_quadratic.bounding_box_oriented_approx()
    );
    assert_eq!(
        cubic.bounding_box_oriented_approx(),
        generic_cubic.bounding_box_oriented_approx()
    );
    for budget in [0, 1, 8, 64, 256] {
        assert_eq!(
            quadratic.bounding_box_oriented_approx_with_budget(budget),
            generic_quadratic.bounding_box_oriented_approx_with_budget(budget)
        );
        assert_eq!(
            cubic.bounding_box_oriented_approx_with_budget(budget),
            generic_cubic.bounding_box_oriented_approx_with_budget(budget)
        );
    }
}

#[test]
fn specialized_beziers_match_generic_boxes_and_diagnostics() {
    check_specialized([[0.0_f32, 0.0], [1.0, 3.0], [4.0, -2.0], [5.0, 1.0]].map(ArrayPoint::from));
    check_specialized(
        [
            [0.0_f64; 3],
            [1.0, 3.0, 2.0],
            [4.0, -2.0, 1.0],
            [5.0, 1.0, 4.0],
        ]
        .map(ArrayPoint::from),
    );
    check_specialized(
        [
            [0.0_f64; 4],
            [1.0, 3.0, 2.0, 4.0],
            [4.0, -2.0, 1.0, -1.0],
            [5.0, 1.0, 4.0, 2.0],
        ]
        .map(ArrayPoint::from),
    );
    check_specialized([[-2.0_f64], [3.0], [1.0], [0.0]].map(ArrayPoint::from));
    check_specialized([ArrayPoint::from([2.0_f64, -3.0, 4.0]); 4]);
}

#[test]
fn specialized_beziers_forward_generic_errors() {
    check_specialized([ArrayPoint::<f64, 0>::from([]); 4]);
    let sub = f64::from_bits(1);
    let controls =
        [[1e300_f64, 0.0], [1e300, sub], [1e300, 0.0], [1e300, sub]].map(ArrayPoint::from);
    assert_eq!(
        Bezier::new(controls)
            .unwrap()
            .bounding_box_oriented_approx(),
        Err(OrientedBoundingBoxError::NonFiniteComputation)
    );
    check_specialized(controls);
}

fn check_faces<const D: usize, const N: usize>(controls: [[f64; D]; N])
where
    [(); <ArrayPoint<f64, D> as Point>::DIM]: Sized,
{
    let curve = Bezier::new(controls.map(ArrayPoint::from)).unwrap();
    let fit = curve.bounding_box_oriented_approx_with_budget(256).unwrap();
    assert!(fit.evaluations() <= 256);
    let center = fit.center();
    let half = fit.half_extents();
    for (i, &extent) in half.iter().enumerate() {
        let axis = fit.axes()[i];
        let projected = controls.map(|p| (0..D).map(|j| axis[j] * (p[j] - center[j])).sum());
        let (lo, hi) = reference::range(projected);
        // Check the fitted faces directly, without contains()'s roundoff allowance.
        assert!(lo >= -extent, "{D}D minimum {lo} < {}", -extent);
        assert!(hi <= extent, "{D}D maximum {hi} > {extent}");
        for j in 0..D {
            let product: f64 = (0..D).map(|k| axis[k] * fit.axes()[j][k]).sum();
            let expected = if i == j { 1.0 } else { 0.0 };
            assert!((product - expected).abs() < 32.0 * D as f64 * f64::EPSILON);
        }
    }
    let baseline = curve.bounding_box_oriented_approx_with_budget(1).unwrap();
    assert!(fit.log_volume() <= baseline.log_volume());
    for i in 0..=128 {
        assert!(fit.contains(&curve.eval(i as f64 / 128.0)));
    }
    for mask in 0..1_usize << D {
        assert!(fit.contains(&ArrayPoint::from(
            fit.corner(array::from_fn(|i| mask & (1 << i) != 0))
        )));
    }
}

#[test]
fn polynomial_oracle_has_known_interior_and_repeated_extrema() {
    for (actual, expected) in [
        (reference::range([0.0, 2.0, 0.0]), (0.0, 1.0)),
        (reference::range([0.0, 1.0, 1.0, 0.0]), (0.0, 0.75)),
        (
            reference::range([0.0, 0.0, 8.0 / 3.0, 0.0, 0.0]),
            (0.0, 1.0),
        ),
        (reference::range([1.0, -1.0, 1.0, -1.0, 1.0]), (0.0, 1.0)),
    ] {
        assert!((actual.0 - expected.0).abs() < 1e-12);
        assert!((actual.1 - expected.1).abs() < 1e-12);
    }
}

#[test]
fn independent_cubic_faces_in_2d_3d_and_4d() {
    let mut state = 0x51_de;
    for _ in 0..16 {
        check_faces(array::from_fn::<_, 4, _>(|_| {
            array::from_fn::<_, 2, _>(|_| random(&mut state))
        }));
        check_faces(array::from_fn::<_, 4, _>(|_| {
            array::from_fn::<_, 3, _>(|_| random(&mut state))
        }));
        check_faces(array::from_fn::<_, 4, _>(|_| {
            array::from_fn::<_, 4, _>(|_| random(&mut state))
        }));
    }
}

#[test]
fn full_dimensional_4d_quartics_against_power_basis_root_isolation() {
    let mut state = 0x4d_cafe;
    for _ in 0..16 {
        check_faces(array::from_fn::<_, 5, _>(|_| {
            array::from_fn::<_, 4, _>(|_| random(&mut state))
        }));
    }
}

#[test]
fn fixed_fixture_quality_against_independent_2d_reference() {
    let frame = rotated(identity(), 0, 1, 0.63);
    let controls =
        [[0.0, 0.0], [1.0, 3.0], [3.0, 3.0], [4.0, 0.0]].map(|p| frame.map(|a| dot(&a, &p)));
    let curve = Bezier::new(controls.map(ArrayPoint::from)).unwrap();
    let fit = curve.bounding_box_oriented_approx_with_budget(256).unwrap();
    let numerical = reference::reference(&controls);
    let numerical = numerical.min(reference::reference_at(&controls, *fit.axes()));
    assert!((fit.log_volume() - numerical).exp() < 1.0001);
}

fn check_samples<T: num_traits::Float, const D: usize, const N: usize>(controls: [[T; D]; N])
where
    [(); <ArrayPoint<T, D> as Point>::DIM]: Sized,
{
    let curve = Bezier::new(controls.map(ArrayPoint::from)).unwrap();
    let fit = curve.bounding_box_oriented_approx_with_budget(128).unwrap();
    for i in 0..=128 {
        let t = <T as num_traits::NumCast>::from(i as f64 / 128.0).unwrap();
        assert!(fit.contains(&curve.eval(t)));
    }
}

#[test]
fn budgets_determinism_constants_lines_and_closed_curves() {
    let curve = Bezier::new(
        [
            [0.0_f64, 0.0, 0.0],
            [1.0, 3.0, 2.0],
            [4.0, -2.0, 1.0],
            [5.0, 1.0, 4.0],
        ]
        .map(ArrayPoint::from),
    )
    .unwrap();
    for budget in [1, 2, 3, 4, 8, 64, 256] {
        let a = curve
            .bounding_box_oriented_approx_with_budget(budget)
            .unwrap();
        let b = curve
            .bounding_box_oriented_approx_with_budget(budget)
            .unwrap();
        assert_eq!(a, b);
        assert!(a.evaluations() <= budget);
        for i in 0..=128 {
            assert!(a.contains(&curve.eval(i as f64 / 128.0)));
        }
    }
    let constant = Bezier::new([ArrayPoint::from([2.0_f64, -3.0, 4.0, 5.0])]).unwrap();
    let fit = constant
        .bounding_box_oriented_approx_with_budget(1)
        .unwrap();
    assert_eq!(fit.evaluations(), 0);
    assert_eq!(fit.half_extents(), [0.0; 4]);
    assert_eq!(fit.center(), [2.0, -3.0, 4.0, 5.0]);
    assert_eq!(fit.log_volume(), f64::NEG_INFINITY);
    check_samples([[-2.0_f64], [3.0]]);
    check_samples([[0.0_f64; 4], [1.0, 2.0, 3.0, 4.0]]);
    check_faces([[0.0_f64; 3], [3.0, 4.0, 1.0], [-3.0, 4.0, -1.0], [0.0; 3]]);
    check_faces([[0.0_f64; 2], [4.0, 3.0], [-3.0, -2.0], [1.0, 1.0]]);
}

fn check_prefix<const D: usize, const N: usize>(controls: [[f64; D]; N])
where
    [(); <ArrayPoint<f64, D> as Point>::DIM]: Sized,
{
    let curve = Bezier::new(controls.map(ArrayPoint::from)).unwrap();
    let mut best = f64::INFINITY;
    let mut preprocessing = None;
    for budget in [1, 2, 3, 4, 5, 8, 17, 63, 64, 65, 256, 1024] {
        let fit = curve
            .bounding_box_oriented_approx_with_budget(budget)
            .unwrap();
        assert!(fit.evaluations() <= budget);
        assert!(
            fit.log_volume() <= best,
            "{D}D budget {budget} worsened the incumbent"
        );
        best = fit.log_volume();
        if budget > 1 {
            let info = fit.diagnostics().unwrap();
            assert!((2..=65).contains(&info.representative_points));
            assert!(info.model_log_volume_gap.is_finite() && info.model_log_volume_gap >= 0.0);
            assert!(info.point_set_evaluations <= 5000);
            let fixed = (
                info.representative_points,
                info.point_set_evaluations,
                info.representative_error,
            );
            if let Some(previous) = preprocessing {
                assert_eq!(previous, fixed);
            }
            preprocessing = Some(fixed);
        }
    }
}

#[test]
fn budgets_extend_a_fixed_prefix_in_every_dimension() {
    check_prefix([[0.0_f64, 0.0], [3.0, 4.0], [2.0, -3.0], [6.0, 1.0]]);
    check_prefix([
        [0.0_f64; 3],
        [1.0, 3.0, 2.0],
        [4.0, -2.0, 1.0],
        [5.0, 1.0, 4.0],
    ]);
    check_prefix([
        [0.0_f64; 4],
        [1.0, 3.0, 2.0, 4.0],
        [4.0, -2.0, 1.0, -1.0],
        [5.0, 1.0, 4.0, 2.0],
        [3.0, -1.0, 3.0, 1.0],
    ]);
}

#[test]
fn tiny_translated_and_large_geometry_f32_and_f64() {
    for scale in [1e-20_f32, 1.0, 1e20] {
        check_samples(
            [
                [0.0_f32; 3],
                [1.0, 3.0, 2.0],
                [4.0, -2.0, 1.0],
                [5.0, 1.0, 4.0],
            ]
            .map(|p| p.map(|x| x * scale)),
        );
    }
    for scale in [1e-200_f64, 1.0, 1e200] {
        check_samples(
            [
                [0.0_f64; 3],
                [1.0, 3.0, 2.0],
                [4.0, -2.0, 1.0],
                [5.0, 1.0, 4.0],
            ]
            .map(|p| p.map(|x| x * scale)),
        );
    }
    check_samples(
        [
            [0.0_f64; 3],
            [1.0, 3.0, 2.0],
            [4.0, -2.0, 1.0],
            [5.0, 1.0, 4.0],
        ]
        .map(|p| p.map(|x| x + 1e12)),
    );
    check_samples(
        [
            [0.0_f32; 3],
            [1.0, 3.0, 2.0],
            [4.0, -2.0, 1.0],
            [5.0, 1.0, 4.0],
        ]
        .map(|p| p.map(|x| x + 1e6)),
    );
    let sub = f64::from_bits(1);
    check_samples([[0.0_f64, 0.0], [sub, 2.0 * sub], [2.0 * sub, sub]]);
}

#[test]
fn higher_degree_and_dimension_without_special_cases() {
    let mut state = 0x0d69;
    check_samples(array::from_fn::<_, 9, _>(|_| {
        array::from_fn::<_, 6, _>(|_| random(&mut state))
    }));
}

#[test]
fn a_line_has_known_zero_minimum_volume_in_4d() {
    let curve = Bezier::new([
        ArrayPoint::from([0.0_f64; 4]),
        ArrayPoint::from([1.0, 2.0, 3.0, 4.0]),
    ])
    .unwrap();
    let fit = curve.bounding_box_oriented_approx_with_budget(256).unwrap();
    let mut half = fit.half_extents();
    half.sort_unstable_by(f64::total_cmp);
    // The mathematical optimum has three zero widths and length sqrt(30).
    // Numerical padding keeps the reported widths positive.
    assert!((half[3] - 30.0_f64.sqrt() / 2.0).abs() < 1e-9);
    assert!(half[..3].iter().all(|&x| x < 1e-9));
}

#[test]
fn invalid_requests_and_unrepresentable_padding() {
    let curve = Bezier::new([ArrayPoint::<f64, 0>::from([])]).unwrap();
    assert_eq!(
        curve.bounding_box_oriented_approx_with_budget(1),
        Err(OrientedBoundingBoxError::ZeroDimension)
    );
    let curve = Bezier::new([ArrayPoint::from([0.0_f64])]).unwrap();
    assert_eq!(
        curve.bounding_box_oriented_approx_with_budget(0),
        Err(OrientedBoundingBoxError::ZeroBudget)
    );
    let sub = f64::from_bits(1);
    let curve = Bezier::new([
        ArrayPoint::from([1e300_f64, 0.0]),
        ArrayPoint::from([1e300, sub]),
    ])
    .unwrap();
    assert_eq!(
        curve.bounding_box_oriented_approx_with_budget(256),
        Err(OrientedBoundingBoxError::NonFiniteComputation)
    );
}

#[cfg(feature = "nalgebra")]
#[test]
fn nalgebra_points_use_the_same_generic_search() {
    let curve = Bezier::new([
        nalgebra::SVector::<f64, 3>::new(0.0, 0.0, 0.0),
        nalgebra::SVector::<f64, 3>::new(1.0, 3.0, 2.0),
        nalgebra::SVector::<f64, 3>::new(4.0, -2.0, 1.0),
        nalgebra::SVector::<f64, 3>::new(5.0, 1.0, 4.0),
    ])
    .unwrap();
    let fit = curve.bounding_box_oriented_approx_with_budget(256).unwrap();
    for i in 0..=128 {
        assert!(fit.contains(&curve.eval(i as f64 / 128.0)));
    }
    check_specialized(curve.into_control_points());
}
