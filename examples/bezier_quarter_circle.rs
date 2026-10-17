use stroke::{ArrayPoint, Bezier};

fn main() {
    // Classic cubic Bezier approximation of a quarter circle.
    let kappa = 4.0 / 3.0 * (2.0_f64.sqrt() - 1.0);
    let curve = Bezier::<ArrayPoint<f64, 2>, 4>::new([
        ArrayPoint::new([1.0, 0.0]),
        ArrayPoint::new([1.0, kappa]),
        ArrayPoint::new([kappa, 1.0]),
        ArrayPoint::new([0.0, 1.0]),
    ])
    .unwrap();

    let samples = 100;
    let mut max_err = 0.0;
    for i in 0..=samples {
        let t = i as f64 / samples as f64;
        let p = curve.eval(t);
        let r = (p[0] * p[0] + p[1] * p[1]).sqrt();
        let err = (r - 1.0).abs();
        if err > max_err {
            max_err = err;
        }
    }

    println!("kappa={:.9}", kappa);
    println!(
        "max radial error over {} samples: {:.6e}",
        samples + 1,
        max_err
    );
}
