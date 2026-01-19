use stroke::{BezierPath, CubicBezier, LineSegment, PointN, QuadraticBezier};

fn main() {
    let mut path: BezierPath<PointN<f64, 2>, 8> = BezierPath::new();

    path.push_line(LineSegment::new(
        PointN::new([0.0, 0.0]),
        PointN::new([1.0, 0.0]),
    ));

    path.push_quadratic(QuadraticBezier::new(
        PointN::new([1.0, 0.0]),
        PointN::new([1.5, 0.8]),
        PointN::new([2.0, 0.0]),
    ));

    path.push_cubic(CubicBezier::new(
        PointN::new([2.0, 0.0]),
        PointN::new([2.5, -0.6]),
        PointN::new([3.0, 0.8]),
        PointN::new([3.5, 0.0]),
    ));

    let sample = path.eval(0.35).unwrap();
    println!("sample: {:?}", sample);

    let bounds = path.bounding_box().unwrap();
    println!("bounds: {:?}", bounds);
}
