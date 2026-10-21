#![allow(incomplete_features)]
#![feature(generic_const_exprs)]

use num_traits::{Float, NumCast};
use stroke::{
    ArrayPoint, BSpline, Bezier, OrientedBoundingBox, OrientedBoundingBoxDiagnostics,
    OrientedBoundingBoxError, PointIndex,
};

#[test]
fn default_budget_and_public_box_geometry() {
    let curve = Bezier::new(
        [
            ArrayPoint::from([0.0_f64, 0.0]),
            ArrayPoint::from([1.0, 3.0]),
            ArrayPoint::from([3.0, 3.0]),
            ArrayPoint::from([4.0, 0.0]),
        ]
        .map(|p| ArrayPoint::from([0.8 * p[0] - 0.6 * p[1], 0.6 * p[0] + 0.8 * p[1]])),
    )
    .unwrap();
    let bounds = curve.bounding_box_oriented_approx().unwrap();
    assert_eq!(
        bounds,
        curve.bounding_box_oriented_approx_with_budget(256).unwrap()
    );
    assert!(bounds.evaluations() <= 256);
    let diagnostics: Option<OrientedBoundingBoxDiagnostics<f64>> = bounds.diagnostics();
    assert!(diagnostics.is_some());
    assert!(bounds.contains(&ArrayPoint::from(bounds.center())));
    for upper in [[false, false], [false, true], [true, false], [true, true]] {
        assert!(bounds.contains(&ArrayPoint::from(bounds.corner(upper))));
    }
    let extents = bounds.half_extents();
    let center = bounds.center();
    let outside: [f64; 2] =
        core::array::from_fn(|i| center[i] + bounds.axes()[0][i] * (2.0 * extents[0] + 1.0));
    assert!(!bounds.contains(&ArrayPoint::from(outside)));
    let log_area = (2.0 * extents[0]).ln() + (2.0 * extents[1]).ln();
    assert!((log_area - bounds.log_volume()).abs() <= 16.0 * f64::EPSILON);
    for axis in bounds.axes() {
        assert!((axis[0] * axis[0] + axis[1] * axis[1] - 1.0).abs() <= 8.0 * f64::EPSILON);
    }
}

#[test]
fn constant_box_geometry_and_public_errors() {
    let curve = Bezier::new([ArrayPoint::from([1.0_f64, -2.0, 3.0])]).unwrap();
    let bounds = curve.bounding_box_oriented_approx().unwrap();
    assert_eq!(bounds.center(), [1.0, -2.0, 3.0]);
    assert_eq!(bounds.half_extents(), [0.0; 3]);
    assert_eq!(bounds.corner([true, false, true]), bounds.center());
    assert_eq!(bounds.evaluations(), 0);
    assert_eq!(bounds.log_volume(), f64::NEG_INFINITY);
    let error = curve
        .bounding_box_oriented_approx_with_budget(0)
        .unwrap_err();
    assert_eq!(error, OrientedBoundingBoxError::ZeroBudget);
    fn is_error(_: &dyn core::error::Error) {}
    is_error(&error);
}

fn validate_spline<
    T: Float,
    P: PointIndex<Scalar = T>,
    const K: usize,
    const C: usize,
    const D: usize,
>(
    curve: &BSpline<P, K, C, D>,
    budget: usize,
) -> OrientedBoundingBox<T, { P::DIM }>
where
    [(); D + 1]: Sized,
    [(); C + 64]: Sized,
    [(); P::DIM]: Sized,
{
    let bounds = curve
        .bounding_box_oriented_approx_with_budget(budget)
        .unwrap();
    assert!(bounds.evaluations() <= budget);
    let baseline = curve.bounding_box_oriented_approx_with_budget(1).unwrap();
    assert!(bounds.log_volume() <= baseline.log_volume());
    for (_, span) in curve.bezier_spans() {
        // Each span includes both endpoints, including the left-hand limit at
        // jumps that a uniform sample of BSpline::eval would miss.
        for i in 0..=64 {
            let t = <T as NumCast>::from(i).unwrap() / <T as NumCast>::from(64).unwrap();
            assert!(bounds.contains(&span.eval(t)), "span sample {i}/64 escaped");
        }
    }
    let (start, end) = curve.knot_domain();
    for i in 0..=256 {
        let u = <T as NumCast>::from(i).unwrap() / <T as NumCast>::from(256).unwrap();
        let t = if i == 0 {
            start
        } else if i == 256 {
            end
        } else {
            (start * (T::one() - u) + end * u).max(start).min(end)
        };
        assert!(
            bounds.contains(&curve.eval(t).unwrap()),
            "de Boor sample {i}/256 escaped"
        );
    }
    for i in 0..P::DIM {
        for j in 0..P::DIM {
            let product = (0..P::DIM).fold(T::zero(), |sum, k| {
                sum + bounds.axes()[i][k] * bounds.axes()[j][k]
            });
            let expected = if i == j { T::one() } else { T::zero() };
            assert!(
                (product - expected).abs()
                    <= <T as NumCast>::from(16 * P::DIM).unwrap() * T::epsilon()
            );
        }
    }
    bounds
}

#[test]
fn spline_default_budget_and_single_span_bezier_equivalence() {
    let controls = [[0.0_f64, 0.0], [1.0, 3.0], [3.0, -2.0], [4.0, 1.0]].map(ArrayPoint::from);
    let spline = BSpline::<_, 8, 4, 3>::new_clamped_uniform(controls).unwrap();
    let bounds = validate_spline(&spline, 256);
    assert_eq!(bounds, spline.bounding_box_oriented_approx().unwrap());
    let bezier = Bezier::new(controls)
        .unwrap()
        .bounding_box_oriented_approx()
        .unwrap();
    // The extraction margin differs, but the same cubic should fit equally well.
    assert!((bounds.log_volume() - bezier.log_volume()).abs() < 1e-9);
    for budget in [1, 2, 3, 17, 64] {
        let a = spline
            .bounding_box_oriented_approx_with_budget(budget)
            .unwrap();
        let b = spline
            .bounding_box_oriented_approx_with_budget(budget)
            .unwrap();
        assert_eq!(a, b);
        assert!(a.evaluations() <= budget);
    }
}

#[test]
fn unclamped_and_partly_clamped_splines_bound_the_active_curve() {
    let curve = BSpline::<_, 6, 3, 2>::new(
        [0.0_f64, 1.0, 2.0, 3.0, 4.0, 5.0],
        [[0.0, 100.0], [1.0, 0.0], [2.0, 100.0]].map(ArrayPoint::from),
    )
    .unwrap();
    let bounds = validate_spline(&curve, 128);
    assert!(!bounds.contains(curve.control_points().next().unwrap()));
    assert!(!bounds.contains(&ArrayPoint::from([1.0, 0.0])));
    let controls = [
        [0.0_f64, 0.0],
        [1.0, 3.0],
        [2.0, -2.0],
        [4.0, 2.0],
        [5.0, 0.0],
    ]
    .map(ArrayPoint::from);
    for knots in [
        [0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0],
        [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 5.0, 5.0],
    ] {
        validate_spline(&BSpline::<_, 8, 5, 2>::new(knots, controls).unwrap(), 128);
    }
}

#[test]
fn spline_jumps_include_every_one_sided_endpoint() {
    let controls = [[-10.0_f64, 0.0], [-9.0, 1.0], [9.0, -1.0], [10.0, 0.0]].map(ArrayPoint::from);
    let curve = BSpline::<_, 6, 4, 1>::new([0.0, 0.0, 0.5, 0.5, 1.0, 1.0], controls).unwrap();
    let bounds = validate_spline(&curve, 256);
    for point in controls {
        assert!(bounds.contains(&point));
    }
    assert_eq!(curve.eval(0.5).unwrap(), controls[2]);
    assert!(bounds.contains(&controls[1])); // The left-hand limit, which eval() skips at this jump.
    let cubic_controls = [
        [-10.0_f64, 0.0],
        [-8.0, 2.0],
        [-12.0, 3.0],
        [-9.0, 1.0],
        [9.0, -1.0],
        [8.0, -4.0],
        [12.0, 1.0],
        [10.0, 0.0],
    ]
    .map(ArrayPoint::from);
    let cubic = BSpline::<_, 12, 8, 3>::new(
        [0.0, 0.0, 0.0, 0.0, 0.5, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0, 1.0],
        cubic_controls,
    )
    .unwrap();
    let cubic_bounds = validate_spline(&cubic, 128);
    assert!(cubic_bounds.contains(&cubic_controls[3]));
    assert!(cubic_bounds.contains(&cubic_controls[4]));
    let steps = BSpline::<_, 4, 3, 0>::new(
        [0.0_f64, 1.0, 2.0, 3.0],
        [[-2.0, 1.0], [0.0, -1.0], [3.0, 2.0]].map(ArrayPoint::from),
    )
    .unwrap();
    validate_spline(&steps, 128);
}

#[test]
fn constant_splines_and_invalid_requests() {
    let point = ArrayPoint::from([2.0_f64, -3.0, 4.0]);
    let curve = BSpline::<_, 9, 5, 3>::new_clamped_uniform([point; 5]).unwrap();
    let bounds = validate_spline(&curve, 64);
    assert_eq!(bounds.center(), point.into_array());
    assert_eq!(bounds.half_extents(), [0.0; 3]);
    assert_eq!(bounds.evaluations(), 0);
    assert_eq!(bounds.log_volume(), f64::NEG_INFINITY);
    assert_eq!(
        curve
            .bounding_box_oriented_approx_with_budget(0)
            .unwrap_err(),
        OrientedBoundingBoxError::ZeroBudget
    );
    let empty =
        BSpline::<_, 2, 1, 0>::new([0.0_f64, 1.0], [ArrayPoint::<f64, 0>::from([])]).unwrap();
    assert_eq!(
        empty.bounding_box_oriented_approx().unwrap_err(),
        OrientedBoundingBoxError::ZeroDimension
    );
    let point_curve = BSpline::<_, 2, 1, 0>::new([0.0, 1.0], [point]).unwrap();
    validate_spline(&point_curve, 64);
    let unrepresentable = BSpline::<_, 4, 2, 1>::new(
        [0.0_f64, 0.0, 1.0, 1.0],
        [[f64::MAX, 0.0], [f64::MAX, f64::from_bits(1)]].map(ArrayPoint::from),
    )
    .unwrap();
    assert_eq!(
        unrepresentable.bounding_box_oriented_approx().unwrap_err(),
        OrientedBoundingBoxError::NonFiniteComputation
    );
    assert_eq!(
        unrepresentable
            .bounding_box_oriented_approx_with_budget(0)
            .unwrap_err(),
        OrientedBoundingBoxError::ZeroBudget
    );
}

macro_rules! spline_numeric_cases {
    ($name:ident, $float:ty, $small:expr, $large:expr, $offset:expr) => {
        #[test]
        fn $name() {
            for scale in [$small, 1.0, $large] {
                for offset in [0.0, $offset] {
                    let controls = [[0.0 as $float, 0.0], [1.0, 3.0], [3.0, -2.0], [4.0, 1.0]]
                        .map(|p| ArrayPoint::from(p.map(|x| x * scale + offset)));
                    let curve = BSpline::<_, 8, 4, 3>::new_clamped_uniform(controls).unwrap();
                    validate_spline(&curve, 64);
                }
            }
            for width in [
                <$float>::from_bits(2),
                <$float>::MIN_POSITIVE,
                $small,
                1.0,
                $large,
            ] {
                let curve = BSpline::<_, 4, 2, 1>::new(
                    [0.0 as $float, 0.0, width, width],
                    [[0.0, 0.0], [2.0, 3.0]].map(ArrayPoint::from),
                )
                .unwrap();
                validate_spline(&curve, 64);
            }
            let curve = BSpline::<_, 4, 2, 1>::new(
                [-<$float>::MAX, -<$float>::MAX, <$float>::MAX, <$float>::MAX],
                [[0.0 as $float, 0.0], [2.0, 3.0]].map(ArrayPoint::from),
            )
            .unwrap();
            validate_spline(&curve, 64);
            let curve = BSpline::<_, 4, 2, 1>::new(
                [0.0 as $float, 0.0, 1.0, 1.0],
                [[-<$float>::MAX, 0.0], [<$float>::MAX, 1.0]].map(ArrayPoint::from),
            )
            .unwrap();
            validate_spline(&curve, 64);
        }
    };
}

spline_numeric_cases!(
    spline_f64_scales_translations_and_knot_ranges,
    f64,
    1e-200_f64,
    1e200_f64,
    1e12_f64
);
spline_numeric_cases!(
    spline_f32_scales_translations_and_knot_ranges,
    f32,
    1e-20_f32,
    1e20_f32,
    1e6_f32
);

#[test]
fn splines_support_higher_dimensions_and_repeated_knots() {
    let controls = [
        [0.0_f64, 0.0, 0.0, 0.0],
        [1.0, 3.0, 2.0, 4.0],
        [4.0, -2.0, 1.0, -1.0],
        [5.0, 1.0, 4.0, 2.0],
        [3.0, -1.0, 3.0, 1.0],
    ]
    .map(ArrayPoint::from);
    let curve = BSpline::<_, 9, 5, 3>::new_clamped_uniform(controls).unwrap();
    validate_spline(&curve, 128);
    let repeated = BSpline::<_, 9, 6, 2>::new(
        [0.0_f64, 0.0, 0.0, 1.0, 1.0, 2.0, 3.0, 3.0, 3.0],
        [
            [0.0, 0.0],
            [1.0, 3.0],
            [2.0, -1.0],
            [3.0, 2.0],
            [4.0, -2.0],
            [5.0, 0.0],
        ]
        .map(ArrayPoint::from),
    )
    .unwrap();
    validate_spline(&repeated, 128);
    let quintic = BSpline::<_, 14, 8, 5>::new_clamped_uniform(
        [
            [0.0_f64, 0.0],
            [0.5, 4.0],
            [2.0, -3.0],
            [3.0, 4.0],
            [4.5, -2.0],
            [5.0, 1.0],
            [6.0, -1.0],
            [7.0, 0.0],
        ]
        .map(ArrayPoint::from),
    )
    .unwrap();
    validate_spline(&quintic, 128);
}

#[cfg(feature = "nalgebra")]
#[test]
fn spline_oriented_bounds_integrate_with_nalgebra() {
    let curve = BSpline::<_, 9, 5, 3>::new_clamped_uniform([
        nalgebra::SVector::<f64, 3>::new(0.0, 0.0, 0.0),
        nalgebra::SVector::<f64, 3>::new(1.0, 3.0, 2.0),
        nalgebra::SVector::<f64, 3>::new(4.0, -2.0, 1.0),
        nalgebra::SVector::<f64, 3>::new(5.0, 1.0, 4.0),
        nalgebra::SVector::<f64, 3>::new(3.0, -1.0, 3.0),
    ])
    .unwrap();
    validate_spline(&curve, 128);
}
