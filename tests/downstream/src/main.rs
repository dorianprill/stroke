//! Built from a temporary directory in CI, without inheriting repository Cargo config.
#![allow(incomplete_features)]
#![feature(generic_const_exprs)]

use stroke::{ArrayPoint, BSpline, Bezier, BezierPath, CubicBezier, RootFindingError};

fn main() -> Result<(), Box<dyn core::error::Error>> {
    type P = ArrayPoint<f64, 2>;
    let points = [
        P::from([0.0, 0.0]),
        P::from([1.0, 2.0]),
        P::from([2.0, 1.0]),
        P::from([3.0, 0.0]),
    ];
    let curve = Bezier::new(points)?;
    let cubic = CubicBezier::from(curve);
    let mut path = BezierPath::<P, 2>::new();
    path.try_push(cubic.into()).expect("available capacity");
    assert!(path.speed(0.5).unwrap().is_finite());
    assert_eq!(curve.subcurve(0.0..1.0)?, curve);
    assert!(curve.bounding_box().intersects(&cubic.control_box()));
    let spline = BSpline::<P, 8, 4, 3>::new_clamped_uniform(points)?;
    assert!(spline.speed(0.5)?.is_finite());
    let _: Option<RootFindingError> = None;
    #[cfg(feature = "nalgebra")]
    {
        use nalgebra::SVector;
        let points = points.map(|p| SVector::<f64, 2>::from_row_slice(p.as_array()));
        let curve = Bezier::new(points)?;
        assert!(curve.speed(0.5).is_finite());
        assert!(curve.bounding_box().contains(&curve.eval(0.5)));
    }
    Ok(())
}
