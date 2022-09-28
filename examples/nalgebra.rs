// build this example with cargo test --features nalgebra

extern crate stroke;
use stroke::Bezier;

extern crate nalgebra;
use nalgebra::Point;

fn main () {
    let p: Point<f32, 2>;

    let points: [Point<f32,2>; 6] = [
            Point::new(0f32, 1.77f32),
            Point::new(1.1f32, -1f32),
            Point::new(4.3f32, 3f32),
            Point::new(3.2f32, -4f32),
            Point::new(7.3f32, 2.7f32),
            Point::new(8.9f32, 1.7f32),
    ];

    let curve: Bezier<Point<f32, 2>, 6> = Bezier::new(points);
    curve.eval(0.5);
}