use stroke::{ArrayPoint, BezierPath, CubicBezier, LineSegment, QuadraticBezier};

fn main() {
    let mut path: BezierPath<ArrayPoint<f64, 2>, 8> = BezierPath::new();

    path.push_line(
        LineSegment::new(ArrayPoint::new([0.0, 0.0]), ArrayPoint::new([1.0, 0.0])).unwrap(),
    );

    path.push_quadratic(
        QuadraticBezier::new(
            ArrayPoint::new([1.0, 0.0]),
            ArrayPoint::new([1.5, 0.8]),
            ArrayPoint::new([2.0, 0.0]),
        )
        .unwrap(),
    );

    path.push_cubic(
        CubicBezier::new(
            ArrayPoint::new([2.0, 0.0]),
            ArrayPoint::new([2.5, -0.6]),
            ArrayPoint::new([3.0, 0.8]),
            ArrayPoint::new([3.5, 0.0]),
        )
        .unwrap(),
    );

    let sample = path.eval(0.35).unwrap();
    println!("sample: {:?}", sample);

    let bounds = path.bounding_box().unwrap();
    println!("bounds: {:?}", bounds);
}
