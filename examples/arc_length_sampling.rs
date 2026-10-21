//! Place six samples at equal arc-length intervals along a cubic Bezier curve.
//! The printed gaps are straight-line distances, which can differ as the curve bends.

use stroke::{ArrayPoint, Bezier, PointNorm};

fn main() {
    let curve = Bezier::<ArrayPoint<f64, 2>, 4>::new([
        ArrayPoint::new([0.0, 0.0]),
        ArrayPoint::new([1.0, 2.0]),
        ArrayPoint::new([3.0, -1.0]),
        ArrayPoint::new([4.0, 0.0]),
    ])
    .unwrap();

    let steps = 6;
    let accuracy = 1e-10;
    let total = curve.arclen(accuracy).unwrap();
    println!("total length (approx): {:.5}", total);

    let mut prev: Option<ArrayPoint<f64, 2>> = None;
    for i in 0..steps {
        let s = total * (i as f64 / (steps - 1) as f64);
        let t = curve.t_at_length(s, accuracy).unwrap();
        let p = curve.point_at_length(s, accuracy).unwrap();
        let gap = prev.map(|q| (p - q).norm());
        match gap {
            Some(d) => println!("i={}  s={:.4}  t={:.4}  p={:?}  gap={:.4}", i, s, t, p, d),
            None => println!("i={}  s={:.4}  t={:.4}  p={:?}", i, s, t, p),
        }
        prev = Some(p);
    }
}
