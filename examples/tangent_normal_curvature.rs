use stroke::{Bezier, PointN};

fn main() {
    let curve = Bezier::<PointN<f64, 2>, 4>::new([
        PointN::new([-1.0, 0.0]),
        PointN::new([-0.5, 1.5]),
        PointN::new([0.5, -1.5]),
        PointN::new([1.0, 0.0]),
    ]);

    println!("t    p(x,y)         tan(x,y)       kappa   normal sketch");
    let samples = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0];
    for t in samples {
        let p = curve.eval(t);
        let tan = curve.tangent(t);
        let kappa = curve.curvature(t);
        let arrow = match curve.normal(t) {
            Some(n) => arrow_from_vec(n[0], n[1]),
            None => '?',
        };
        let sketch = hedgehog_line(p[0], arrow);
        println!(
            "{:>3.1}  ({:>6.3},{:>6.3})  ({:>6.3},{:>6.3})  {:>6.3}  {}",
            t, p[0], p[1], tan[0], tan[1], kappa, sketch
        );
    }
}

fn arrow_from_vec(x: f64, y: f64) -> char {
    if x.abs() >= y.abs() {
        if x >= 0.0 {
            '>'
        } else {
            '<'
        }
    } else if y >= 0.0 {
        '^'
    } else {
        'v'
    }
}

fn hedgehog_line(x: f64, arrow: char) -> String {
    let width = 24usize;
    let mut line = vec!['.'; width];
    let mut pos = ((x + 1.0) * 0.5 * (width as f64 - 1.0)).round() as isize;
    if pos < 0 {
        pos = 0;
    }
    if pos >= width as isize {
        pos = width as isize - 1;
    }
    line[pos as usize] = arrow;
    line.into_iter().collect()
}
