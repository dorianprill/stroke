use stroke::{ArrayPoint, BSpline, BSplinePath};

fn main() {
    let knots: [f64; 4] = [0.0, 0.0, 1.0, 1.0];

    let seg1: BSpline<ArrayPoint<f64, 2>, 4, 2, 1> = BSpline::new(
        knots,
        [ArrayPoint::new([0.0, 0.0]), ArrayPoint::new([1.0, 0.0])],
    )
    .unwrap();

    let seg2: BSpline<ArrayPoint<f64, 2>, 4, 2, 1> = BSpline::new(
        knots,
        [ArrayPoint::new([1.0, 0.0]), ArrayPoint::new([1.0, 1.0])],
    )
    .unwrap();

    let mut path: BSplinePath<ArrayPoint<f64, 2>, 4, 2, 1, 4> = BSplinePath::new();
    path.push(seg1);
    path.push(seg2);

    let sample = path.eval(0.6).unwrap();
    println!("sample: {:?}", sample);

    let bounds = path.bounding_box().unwrap();
    println!("bounds: {:?}", bounds);
}
