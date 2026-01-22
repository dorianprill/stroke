use stroke::{Bezier, PointN, PointNorm};

fn main() {
    let curve = Bezier::<PointN<f64, 2>, 4>::new([
        PointN::new([0.0, 0.0]),
        PointN::new([1.0, 2.0]),
        PointN::new([3.0, -1.0]),
        PointN::new([4.0, 0.0]),
    ]);

    let steps = 6;
    let nsteps = 512;
    let total = curve.arclen(nsteps);
    println!("total length (approx): {:.5}", total);

    let mut prev: Option<PointN<f64, 2>> = None;
    for i in 0..steps {
        let s = total * (i as f64 / (steps - 1) as f64);
        let t = curve.t_at_length_approx(s, nsteps);
        let p = curve.point_at_length_approx(s, nsteps);
        let gap = prev.map(|q| (p - q).squared_norm().sqrt());
        match gap {
            Some(d) => println!("i={}  s={:.4}  t={:.4}  p={:?}  gap={:.4}", i, s, t, p, d),
            None => println!("i={}  s={:.4}  t={:.4}  p={:?}", i, s, t, p),
        }
        prev = Some(p);
    }
}
