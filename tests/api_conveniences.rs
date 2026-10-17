//! Consumer-facing tests: public imports, trait bounds, and cross-representation behavior.

use stroke::{
    ArcLengthError, ArrayPoint, BSpline, BSplineError, BSplinePath, Bezier, BezierPath,
    BezierSegment, BoundingBox, CubicBezier, CurveConstructionError, CurveParameterError, FindRoot,
    LineSegment, Point, PointNorm, QuadraticBezier, RootFindingError,
};

type P = ArrayPoint<f64, 2>;
type S = ArrayPoint<f64, 1>;

// No Default, Debug, component indexing, or PartialEq: evaluation conveniences
// must remain usable with only the core Point abstraction.
#[derive(Clone, Copy)]
struct Minimal(f64);
impl core::ops::Add for Minimal {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        Self(self.0 + other.0)
    }
}
impl core::ops::Sub for Minimal {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        Self(self.0 - other.0)
    }
}
impl core::ops::Mul<f64> for Minimal {
    type Output = Self;
    fn mul(self, factor: f64) -> Self {
        Self(self.0 * factor)
    }
}
impl Point for Minimal {
    type Scalar = f64;
    const DIM: usize = 1;
    fn is_finite(&self) -> bool {
        self.0.is_finite()
    }
}

#[test]
fn basic_conveniences_do_not_strengthen_point_requirements() {
    let generic = Bezier::new([Minimal(0.0), Minimal(1.0), Minimal(3.0)]).unwrap();
    near(generic.derivative_at(0.5).0, 3.0);
    near(
        generic.subcurve(0.25..0.75).unwrap().eval(0.5).0,
        generic.eval(0.5).0,
    );
    let specialized = QuadraticBezier::from(generic);
    near(specialized.derivative_at(0.5).0, 3.0);
    let spline =
        BSpline::<Minimal, 6, 3, 2>::new_clamped_uniform(generic.into_control_points()).unwrap();
    near(spline.derivative_at(0.5).unwrap().0, 3.0);
    let mapped = spline.try_map_control_points(|p| p * 2.0).unwrap();
    near(
        mapped.eval(0.5).unwrap().0,
        spline.eval(0.5).unwrap().0 * 2.0,
    );
}

fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-11 * expected.abs().max(1.0),
        "{actual} != {expected}"
    );
}

fn point_near(actual: P, expected: P) {
    for axis in 0..2 {
        near(actual[axis], expected[axis]);
    }
}

#[test]
fn points_convert_and_borrow_without_copy_bounds() {
    let point = ArrayPoint::from([String::from("x"), String::from("y")]);
    assert!(core::ptr::eq(point.as_array(), point.as_ref()));
    assert_eq!(
        (&point).into_iter().map(String::as_str).collect::<Vec<_>>(),
        ["x", "y"]
    );
    let array: [String; 2] = point.into();
    assert_eq!(array, ["x", "y"]);
    let empty = ArrayPoint::<u8, 0>::from([]);
    assert_eq!((&empty).into_iter().count(), 0);
    assert_eq!(empty.into_array(), []);
}

#[test]
fn public_errors_integrate_with_standard_error_handling() {
    fn is_error<E: core::error::Error>() {}
    is_error::<ArcLengthError>();
    is_error::<BSplineError>();
    is_error::<CurveConstructionError>();
    is_error::<CurveParameterError>();
    is_error::<RootFindingError>();
    fn construct() -> Result<(), Box<dyn core::error::Error>> {
        LineSegment::new(S::new([0.0]), S::new([1.0]))?;
        BSpline::<S, 4, 2, 1>::new_clamped_uniform([S::new([0.0]), S::new([1.0])])?;
        Ok(())
    }
    construct().unwrap();
    struct Axis;
    impl FindRoot<S> for Axis {
        fn parameter_domain(&self) -> (f64, f64) {
            (0.0, 1.0)
        }
        fn axis_value(&self, t: f64, _: usize) -> Result<f64, RootFindingError> {
            Ok(t)
        }
        fn axis_derivative(&self, _: f64, _: usize) -> Result<f64, RootFindingError> {
            Ok(1.0)
        }
    }
    near(Axis.root_newton_axis(0.4, 0, 0.0, None, None).unwrap(), 0.4);
}

#[test]
fn equivalent_beziers_convert_losslessly_and_borrow_controls() {
    let points = [
        P::new([0.0, 1.0]),
        P::new([2.0, -1.0]),
        P::new([3.0, 2.0]),
        P::new([5.0, 0.0]),
    ];
    let generic = Bezier::new(points).unwrap();
    let cubic = CubicBezier::from(generic);
    assert_eq!(Bezier::from(cubic), generic);
    assert_eq!(cubic.into_control_points(), points);
    assert_eq!(cubic.control_points().copied().collect::<Vec<_>>(), points);
    assert!(core::ptr::eq(
        generic.control_points().next().unwrap(),
        (&generic).into_iter().next().unwrap()
    ));
    let quadratic = QuadraticBezier::new(points[0], points[1], points[2]).unwrap();
    assert_eq!(QuadraticBezier::from(Bezier::from(quadratic)), quadratic);
    assert_eq!(
        quadratic.control_points().copied().collect::<Vec<_>>(),
        &points[..3]
    );
    let line = LineSegment::new(points[0], points[1]).unwrap();
    assert_eq!(LineSegment::from(Bezier::from(line)), line);
    assert_eq!(
        line.control_points().copied().collect::<Vec<_>>(),
        &points[..2]
    );
}

#[test]
fn control_replacement_is_atomic_and_mapping_is_validated() {
    let points = [S::new([0.0]), S::new([1.0]), S::new([2.0]), S::new([3.0])];
    let mut curve = Bezier::new(points).unwrap();
    let mut invalid = points;
    invalid[2] = S::new([f64::NAN]);
    assert_eq!(
        curve.set_control_points(invalid),
        Err(CurveConstructionError::NonFiniteControlPoint { index: 2 })
    );
    assert_eq!(curve.into_control_points(), points);
    assert_eq!(
        curve.try_map_control_points(|_| S::new([f64::INFINITY])),
        Err(CurveConstructionError::NonFiniteControlPoint { index: 0 })
    );
    let mapped = curve
        .try_map_control_points(|p| p * 2.0 + S::new([4.0]))
        .unwrap();
    near(mapped.eval(0.3)[0], curve.eval(0.3)[0] * 2.0 + 4.0);
    curve
        .set_control_points(mapped.into_control_points())
        .unwrap();
    assert_eq!(curve, mapped);
    let mut cubic = CubicBezier::new(points[0], points[1], points[2], points[3]).unwrap();
    let original = cubic;
    assert!(cubic.set_control_points(invalid).is_err());
    assert_eq!(cubic, original);
    let mut quadratic = QuadraticBezier::new(points[0], points[1], points[2]).unwrap();
    let original = quadratic;
    assert!(
        quadratic
            .set_control_points([points[0], points[1], invalid[2]])
            .is_err()
    );
    assert_eq!(quadratic, original);
    let mut line = LineSegment::new(points[0], points[1]).unwrap();
    let original = line;
    assert!(line.set_control_points([invalid[2], points[1]]).is_err());
    assert_eq!(line, original);
    let mut spline = BSpline::<S, 8, 4, 3>::new_clamped_uniform(points).unwrap();
    let original = spline;
    assert_eq!(
        spline.set_control_points(invalid),
        Err(BSplineError::NonFiniteControlPoint { index: 2 })
    );
    assert_eq!(spline, original);
    assert_eq!(
        spline.try_map_control_points(|_| S::new([f64::NEG_INFINITY])),
        Err(BSplineError::NonFiniteControlPoint { index: 0 })
    );
    spline
        .set_control_points(mapped.into_control_points())
        .unwrap();
    assert_eq!(
        spline.knots().copied().collect::<Vec<_>>(),
        original.knots().copied().collect::<Vec<_>>()
    );
    let segment = BezierSegment::from(line);
    assert!(
        segment
            .try_map_control_points(|_| S::new([f64::NAN]))
            .is_err()
    );
}

#[test]
fn uniform_knots_cover_zero_linear_and_cubic_degrees() {
    let zero = BSpline::<S, 4, 3, 0>::new_clamped_uniform([S::new([1.0]); 3]).unwrap();
    assert_eq!(
        zero.knots().copied().collect::<Vec<_>>(),
        [0.0, 1.0 / 3.0, 2.0 / 3.0, 1.0]
    );
    let linear = BSpline::<S, 5, 3, 1>::new_clamped_uniform([S::new([1.0]); 3]).unwrap();
    assert_eq!(
        linear.knots().copied().collect::<Vec<_>>(),
        [0.0, 0.0, 0.5, 1.0, 1.0]
    );
    let cubic = BSpline::<S, 10, 6, 3>::new_clamped_uniform([S::new([2.0]); 6]).unwrap();
    assert_eq!(
        cubic.knots().copied().collect::<Vec<_>>(),
        [0.0, 0.0, 0.0, 0.0, 1.0 / 3.0, 2.0 / 3.0, 1.0, 1.0, 1.0, 1.0]
    );
    assert_eq!(cubic.knot_domain(), (0.0, 1.0));
    assert_eq!(
        BSpline::<S, 4, 1, 1>::new_clamped_uniform([S::new([0.0])]),
        Err(BSplineError::TooFewControlPoints)
    );
    assert_eq!(
        BSpline::<S, 3, 2, 1>::new_clamped_uniform([S::new([0.0]); 2]),
        Err(BSplineError::TooFewKnots)
    );
    assert_eq!(
        BSpline::<S, 5, 2, 1>::new_clamped_uniform([S::new([0.0]); 2]),
        Err(BSplineError::TooManyKnots)
    );
    assert_eq!(
        BSpline::<S, 4, 2, 1>::new_clamped_uniform([S::new([0.0]), S::new([f64::NAN])]),
        Err(BSplineError::NonFiniteControlPoint { index: 1 })
    );
}

#[test]
fn paths_reuse_storage_and_preserve_rejected_segments() {
    let line = LineSegment::new(P::new([0.0, 0.0]), P::new([1.0, 0.0])).unwrap();
    let segment = BezierSegment::from(line);
    let mut path = BezierPath::<P, 1>::new();
    assert_eq!(path.capacity(), 1);
    assert_eq!(path.remaining_capacity(), 1);
    path.try_push(segment).unwrap();
    let copy = path.clone();
    assert_eq!(copy.start(), path.start());
    assert_eq!(path.remaining_capacity(), 0);
    assert_eq!(path.try_push(segment), Err(segment));
    assert_eq!(path.len(), 1);
    assert_eq!(path.pop(), Some(segment));
    assert_eq!(path.pop(), None);
    path.try_push(segment).unwrap();
    path.clear();
    assert!(path.is_empty());
    assert_eq!(path.remaining_capacity(), 1);
    assert!(!format!("{copy:?}").is_empty());
    let mut empty = BezierPath::<P, 0>::new();
    assert_eq!(empty.try_push(segment), Err(segment));
    assert_eq!(empty.pop(), None);
    let spline = BSpline::<P, 4, 2, 1>::new_clamped_uniform(line.into_control_points()).unwrap();
    let mut path = BSplinePath::<P, 4, 2, 1, 1>::new();
    assert_eq!(path.capacity(), 1);
    path.try_push(spline).unwrap();
    let copy = path.clone();
    assert_eq!(copy.start(), path.start());
    assert_eq!(path.try_push(spline), Err(spline));
    assert_eq!(path.pop(), Some(spline));
    assert_eq!(path.pop(), None);
    path.try_push(spline).unwrap();
    path.clear();
    assert_eq!(path.remaining_capacity(), 1);
    assert!(!format!("{copy:?}").is_empty());
    let mut empty = BSplinePath::<P, 4, 2, 1, 0>::new();
    assert_eq!(empty.try_push(spline), Err(spline));
}

#[test]
fn boxes_handle_overlap_touching_degeneracy_and_extreme_centers() {
    let a = BoundingBox::new([0.0, 0.0], [2.0, 4.0]).unwrap();
    let b = BoundingBox::new([2.0, 1.0], [3.0, 3.0]).unwrap();
    assert_eq!(a.min(), [0.0, 0.0]);
    assert_eq!(a.max(), [2.0, 4.0]);
    assert_eq!(a.center(), [1.0, 2.0]);
    assert_eq!(a.extent(), [2.0, 4.0]);
    assert!(a.intersects(&b));
    assert_eq!(a.intersection(&b), BoundingBox::new([2.0, 1.0], [2.0, 3.0]));
    assert_eq!(a.intersection(&b), b.intersection(&a));
    let disjoint = BoundingBox::new([3.0, 0.0], [4.0, 1.0]).unwrap();
    assert!(!a.intersects(&disjoint));
    assert_eq!(a.intersection(&disjoint), None);
    let huge = BoundingBox::new([-f64::MAX, f64::MAX], [f64::MAX, f64::MAX]).unwrap();
    assert_eq!(huge.center(), [0.0, f64::MAX]);
    assert_eq!(huge.extent(), [f64::INFINITY, 0.0]);
    let tiny = f64::from_bits(1);
    assert_eq!(
        BoundingBox::new([tiny], [tiny * 2.0]).unwrap().center(),
        [tiny * 2.0]
    );
    assert_eq!(BoundingBox::new([tiny], [tiny]).unwrap().center(), [tiny]);
    let zero = BoundingBox::<f32, 0>::new([], []).unwrap();
    assert_eq!(zero.intersection(&zero), Some(zero));
    assert_eq!(zero.center(), []);
}

#[test]
fn direct_derivatives_agree_with_analytic_beziers() {
    let curve =
        QuadraticBezier::new(P::new([0.0, 0.0]), P::new([1.0, 2.0]), P::new([3.0, 0.0])).unwrap();
    let generic = Bezier::from(curve);
    let line = LineSegment::new(P::new([0.0, 0.0]), P::new([3.0, 4.0])).unwrap();
    let cubic = CubicBezier::new(
        curve.start(),
        P::new([1.0, 2.0]),
        P::new([2.0, -1.0]),
        curve.end(),
    )
    .unwrap();
    for t in [0.0, 0.2, 0.5, 0.8, 1.0] {
        let expected = P::new([2.0 + 2.0 * t, 4.0 - 8.0 * t]);
        point_near(curve.derivative_at(t), expected);
        point_near(generic.derivative_at(t), expected);
        near(curve.speed(t), expected.norm());
        near(generic.speed(t), expected.norm());
        point_near(cubic.derivative_at(t), cubic.derivative().eval(t));
        near(cubic.speed(t), cubic.derivative_at(t).norm());
        assert_eq!(line.derivative_at(t), P::new([3.0, 4.0]));
        assert_eq!(line.speed(t), 5.0);
        point_near(BezierSegment::from(curve).derivative_at(t), expected);
    }
    let constant = Bezier::new([P::new([4.0, -7.0])]).unwrap();
    assert_eq!(constant.derivative_at(0.3), P::new([0.0, 0.0]));
    assert_eq!(constant.speed(0.3), 0.0);
}

#[test]
fn bspline_derivatives_match_extracted_spans_and_handle_zero_degree() {
    let points = [
        P::new([0.0, 0.0]),
        P::new([1.0, 2.0]),
        P::new([3.0, -1.0]),
        P::new([4.0, 3.0]),
        P::new([5.0, 1.0]),
        P::new([6.0, 0.0]),
    ];
    let clamped = BSpline::<P, 10, 6, 3>::new_clamped_uniform(points).unwrap();
    let unclamped =
        BSpline::<P, 10, 6, 3>::new([0.0, 0.5, 1.0, 2.0, 3.0, 3.0, 5.0, 6.0, 7.0, 9.0], points)
            .unwrap();
    for spline in [clamped, unclamped] {
        for (range, span) in spline.bezier_spans() {
            for local in [0.0, 0.25, 0.5, 0.75] {
                let t = range.start + (range.end - range.start) * local;
                let expected = span.derivative_at(local) * (1.0 / (range.end - range.start));
                point_near(spline.derivative_at(t).unwrap(), expected);
                point_near(
                    spline.derivative_at(t).unwrap(),
                    spline.derivative().eval(t).unwrap(),
                );
                near(spline.speed(t).unwrap(), expected.norm());
            }
        }
        assert_eq!(
            spline.derivative_at(f64::NAN),
            Err(BSplineError::KnotDomainViolation)
        );
        assert_eq!(
            spline.derivative_at(spline.knot_domain().1 + 1.0),
            Err(BSplineError::KnotDomainViolation)
        );
    }
    let zero =
        BSpline::<P, 4, 3, 0>::new_clamped_uniform([points[0], points[1], points[2]]).unwrap();
    for t in [0.0, 0.25, 0.5, 1.0] {
        assert_eq!(zero.derivative_at(t).unwrap(), P::new([0.0, 0.0]));
        assert_eq!(zero.speed(t).unwrap(), 0.0);
    }
    let jump = BSpline::<S, 6, 4, 1>::new(
        [0.0, 0.0, 0.5, 0.5, 1.0, 1.0],
        [S::new([0.0]), S::new([1.0]), S::new([10.0]), S::new([13.0])],
    )
    .unwrap();
    assert_eq!(jump.derivative_at(0.25).unwrap()[0], 2.0);
    assert_eq!(jump.derivative_at(0.5).unwrap()[0], 6.0);
    assert_eq!(jump.derivative_at(1.0).unwrap()[0], 6.0);
}

#[test]
fn path_derivatives_apply_both_parameter_scales() {
    let first = LineSegment::new(P::new([0.0, 0.0]), P::new([1.0, 0.0])).unwrap();
    let second = LineSegment::new(P::new([20.0, 0.0]), P::new([24.0, 0.0])).unwrap();
    let mut path = BezierPath::<P, 2>::new();
    let mut splines = BSplinePath::<P, 4, 2, 1, 2>::new();
    assert_eq!(path.derivative_at(0.0), None);
    assert_eq!(splines.speed(0.0), None);
    path.try_push(first.into()).unwrap();
    path.try_push(second.into()).unwrap();
    splines
        .try_push(BSpline::new([2.0, 2.0, 6.0, 6.0], first.into_control_points()).unwrap())
        .unwrap();
    splines
        .try_push(BSpline::new([-3.0, -3.0, 5.0, 5.0], second.into_control_points()).unwrap())
        .unwrap();
    for (t, expected) in [(0.0, 2.0), (0.25, 2.0), (0.5, 8.0), (0.75, 8.0), (1.0, 8.0)] {
        assert_eq!(path.derivative_at(t), Some(P::new([expected, 0.0])));
        assert_eq!(splines.derivative_at(t), Some(P::new([expected, 0.0])));
        assert_eq!(path.speed(t), Some(expected));
        assert_eq!(splines.speed(t), Some(expected));
    }
    assert_eq!(path.derivative_at(f64::NAN), None);
    assert_eq!(splines.derivative_at(f64::NAN), None);
}

#[test]
fn bezier_subcurves_reparameterize_and_validate_ranges() {
    let curve = Bezier::new([
        P::new([0.0, 0.0]),
        P::new([1.0, 2.0]),
        P::new([3.0, -1.0]),
        P::new([5.0, 0.0]),
    ])
    .unwrap();
    for range in [
        0.0..1.0,
        0.0..0.4,
        0.3..1.0,
        0.25..0.75,
        0.5..0.5,
        1.0..1.0,
        0.0..0.0,
    ] {
        let part = curve.subcurve(range.clone()).unwrap();
        let specialized = CubicBezier::from(curve).subcurve(range.clone()).unwrap();
        let segment = BezierSegment::from(CubicBezier::from(curve))
            .subcurve(range.clone())
            .unwrap();
        for t in [0.0, 0.2, 0.5, 0.8, 1.0] {
            let expected = curve.eval(range.start + (range.end - range.start) * t);
            point_near(part.eval(t), expected);
            point_near(specialized.eval(t), expected);
            point_near(segment.eval(t), expected);
            point_near(
                part.derivative_at(t),
                curve.derivative_at(range.start + (range.end - range.start) * t)
                    * (range.end - range.start),
            );
        }
    }
    assert_eq!(curve.subcurve(0.0..1.0).unwrap(), curve);
    for range in [
        -0.1..0.5,
        0.5..0.4,
        0.0..1.1,
        f64::NAN..1.0,
        0.0..f64::INFINITY,
    ] {
        assert_eq!(
            curve.subcurve(range),
            Err(CurveParameterError::InvalidRange)
        );
    }
    let constant = Bezier::new([P::new([4.0, 7.0])]).unwrap();
    assert_eq!(constant.subcurve(0.3..0.8).unwrap(), constant);
    let extreme = Bezier::new([S::new([-1e308]), S::new([0.0]), S::new([1e308])]).unwrap();
    let part = extreme.subcurve(0.25..0.75).unwrap();
    near(part.start()[0] / 1e308, -0.5);
    near(part.end()[0] / 1e308, 0.5);
    let quadratic =
        QuadraticBezier::new(P::new([0.0, 0.0]), P::new([1.0, 2.0]), P::new([3.0, 0.0])).unwrap();
    point_near(
        quadratic.subcurve(0.2..0.8).unwrap().eval(0.5),
        quadratic.eval(0.5),
    );
    let line = LineSegment::new(P::new([0.0, 0.0]), P::new([4.0, 0.0])).unwrap();
    assert_eq!(
        line.subcurve(0.25..0.75).unwrap().into_control_points(),
        [P::new([1.0, 0.0]), P::new([3.0, 0.0])]
    );
}

#[test]
fn splitting_by_length_handles_clamping_and_failures() {
    let line = LineSegment::new(P::new([0.0, 0.0]), P::new([4.0, 0.0])).unwrap();
    for (length, x) in [(-1.0, 0.0), (1.0, 1.0), (9.0, 4.0)] {
        let (a, b) = line.split_at_length(length, 1e-10).unwrap();
        assert_eq!(a.end(), P::new([x, 0.0]));
        assert_eq!(b.start(), a.end());
    }
    assert_eq!(
        line.split_at_length(f64::NAN, 1e-10),
        Err(ArcLengthError::InvalidParameter)
    );
    assert_eq!(
        line.split_at_length(1.0, 0.0),
        Err(ArcLengthError::InvalidAccuracy)
    );
    let cubic = CubicBezier::new(
        P::new([0.0, 0.0]),
        P::new([1.0, 2.0]),
        P::new([3.0, -1.0]),
        P::new([5.0, 0.0]),
    )
    .unwrap();
    let distance = cubic.arclen(1e-10).unwrap() * 0.3;
    let (a, b) = BezierSegment::from(cubic)
        .split_at_length(distance, 1e-10)
        .unwrap();
    point_near(a.end(), b.start());
    assert!((a.arclen(1e-11).unwrap() - distance).abs() < 1e-10);
    let quadratic = QuadraticBezier::new(cubic.start(), P::new([1.0, 2.0]), cubic.end()).unwrap();
    let (a, b) = quadratic.split_at_length(1.0, 1e-10).unwrap();
    point_near(a.end(), b.start());
    assert!((a.arclen(1e-11).unwrap() - 1.0).abs() < 1e-10);
}

#[test]
fn derivatives_handle_tiny_and_huge_native_domains() {
    let tiny = f64::from_bits(1);
    let spline =
        BSpline::<S, 4, 2, 1>::new([0.0, 0.0, tiny, tiny], [S::new([0.0]), S::new([tiny])])
            .unwrap();
    near(spline.derivative_at(0.0).unwrap()[0], 1.0);
    let huge = BSpline::<S, 4, 2, 1>::new(
        [-1e308, -1e308, 1e308, 1e308],
        [S::new([-1e308]), S::new([1e308])],
    )
    .unwrap();
    near(huge.derivative_at(0.0).unwrap()[0], 1.0);
    let overflow =
        BSpline::<S, 4, 2, 1>::new([0.0, 0.0, tiny, tiny], [S::new([0.0]), S::new([1.0])]).unwrap();
    assert_eq!(
        overflow.derivative_at(0.0),
        Err(BSplineError::NonFiniteComputation)
    );
}

macro_rules! scalar_conveniences {
    ($name:ident, $scalar:ty, $accuracy:expr) => {
        #[test]
        fn $name() {
            type P = ArrayPoint<$scalar, 2>;
            let points = [
                P::from([0.0, 0.0]),
                P::from([1.0, 2.0]),
                P::from([2.0, -1.0]),
                P::from([4.0, 0.0]),
            ];
            let generic = Bezier::new(points).unwrap();
            let cubic = CubicBezier::from(generic);
            let spline = BSpline::<P, 8, 4, 3>::new_clamped_uniform(points).unwrap();
            for t in [0.0, 0.1, 0.25, 0.5, 0.9, 1.0] {
                let derivative = cubic.derivative().eval(t);
                for actual in [
                    generic.derivative_at(t),
                    cubic.derivative_at(t),
                    spline.derivative_at(t).unwrap(),
                ] {
                    assert!((actual - derivative).norm() <= $accuracy);
                }
                assert!((generic.speed(t) - spline.speed(t).unwrap()).abs() <= $accuracy);
            }
            let part = generic.subcurve(0.2..0.8).unwrap();
            for t in [0.0, 0.2, 0.6, 1.0] {
                assert!((part.eval(t) - generic.eval(0.2 + 0.6 * t)).norm() <= $accuracy);
            }
            let mut path = BezierPath::<P, 3>::new();
            let mut splines = BSplinePath::<P, 4, 2, 1, 3>::new();
            for i in 0..3 {
                let line = LineSegment::new(
                    P::from([i as $scalar, 0.0]),
                    P::from([i as $scalar + 1.0, 0.0]),
                )
                .unwrap();
                path.try_push(line.into()).unwrap();
                splines
                    .try_push(
                        BSpline::new([2.0, 2.0, 7.0, 7.0], line.into_control_points()).unwrap(),
                    )
                    .unwrap();
            }
            for t in [0.0, 1.0 / 3.0, 0.5, 2.0 / 3.0, 1.0] {
                for derivative in [
                    path.derivative_at(t).unwrap(),
                    splines.derivative_at(t).unwrap(),
                ] {
                    assert!((derivative[0] - 3.0).abs() <= $accuracy);
                }
            }
            let large = <$scalar>::MAX;
            assert_eq!(
                BoundingBox::new([-large, large], [large, large])
                    .unwrap()
                    .center(),
                [0.0, large]
            );
            let tiny = <$scalar>::from_bits(1);
            let spline = BSpline::<ArrayPoint<$scalar, 1>, 4, 2, 1>::new(
                [0.0, 0.0, tiny, tiny],
                [ArrayPoint::from([0.0]), ArrayPoint::from([tiny])],
            )
            .unwrap();
            assert!((spline.derivative_at(0.0).unwrap()[0] - 1.0).abs() <= $accuracy);
        }
    };
}

scalar_conveniences!(conveniences_f32, f32, 5e-6);
scalar_conveniences!(conveniences_f64, f64, 2e-12);

#[cfg(feature = "nalgebra")]
#[test]
fn conveniences_support_nalgebra() {
    use nalgebra::SVector;
    type V = SVector<f64, 3>;
    let curve = CubicBezier::new(
        V::new(0.0, 0.0, 0.0),
        V::new(1.0, 2.0, 3.0),
        V::new(2.0, -1.0, 4.0),
        V::new(3.0, 0.0, 5.0),
    )
    .unwrap();
    let generic = Bezier::from(curve);
    assert_eq!(CubicBezier::from(generic), curve);
    let spline = BSpline::<V, 8, 4, 3>::new_clamped_uniform(curve.into_control_points()).unwrap();
    for t in [0.0, 0.25, 0.5, 0.75, 1.0] {
        assert!((spline.derivative_at(t).unwrap() - curve.derivative_at(t)).norm() < 1e-12);
        near(spline.speed(t).unwrap(), curve.speed(t));
    }
    let part = curve.subcurve(0.25..0.75).unwrap();
    assert!((part.eval(0.5) - curve.eval(0.5)).norm() < 1e-12);
    assert!(spline.bounding_box().intersects(&curve.bounding_box()));
}
