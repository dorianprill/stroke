use stroke::{BSpline, PointN};

const DEGREE: usize = 3;
const CONTROL_POINTS: usize = 8;
const KNOTS: usize = CONTROL_POINTS + DEGREE + 1;

fn main() {
    // Sample a 1D signal.
    let samples: [f64; CONTROL_POINTS] = core::array::from_fn(|i| {
        let x = i as f64 / (CONTROL_POINTS - 1) as f64;
        let base = (2.0 * std::f64::consts::PI * x).sin();
        let ripple = 0.2 * (12.0 * std::f64::consts::PI * x).sin();
        base + ripple
    });

    // Interpolate the samples with a cubic B-spline.
    // We solve a tridiagonal system to recover control points that force the
    // spline to pass through the sample values at uniform parameter steps.
    let control = interpolate_cubic_uniform(&samples);
    let control_points: [PointN<f64, 1>; CONTROL_POINTS] =
        core::array::from_fn(|i| PointN::new([control[i]]));
    // Use a clamped uniform knot vector so the spline starts/ends at the samples.
    let knots = open_uniform_knots();
    let spline: BSpline<PointN<f64, 1>, KNOTS, CONTROL_POINTS, DEGREE> =
        BSpline::new(knots, control_points).expect("valid knot vector");

    let (kmin, kmax) = spline.knot_domain();
    let width = 80usize;
    let height = 20usize;

    let mut curve_values = Vec::with_capacity(width);
    for i in 0..width {
        let t = kmin + (kmax - kmin) * (i as f64 / (width - 1) as f64);
        curve_values.push(spline.eval(t).unwrap()[0]);
    }

    let mut min = curve_values[0];
    let mut max = curve_values[0];
    for v in curve_values
        .iter()
        .copied()
        .chain(samples.iter().copied())
        .chain(control.iter().copied())
    {
        if v < min {
            min = v;
        }
        if v > max {
            max = v;
        }
    }
    if (max - min).abs() < 1e-12 {
        max = min + 1.0;
    }

    let mut grid = vec![vec![' '; width]; height];
    for (i, value) in curve_values.iter().copied().enumerate() {
        let row = value_to_row(value, min, max, height);
        grid[row][i] = '*';
    }
    // Overlay control points as letters.
    for (i, value) in control.iter().copied().enumerate() {
        let col =
            ((i as f64) * (width - 1) as f64 / (control.len() - 1) as f64).round() as usize;
        let row = value_to_row(value, min, max, height);
        let letter = if i < CONTROL_POINTS / 2 { 'a' } else { 'b' };
        grid[row][col] = letter;
    }

    println!("B-spline interpolation of a 1D signal (* = spline, a/b = control points)");
    println!("domain: [{:.2}, {:.2}], range: [{:.2}, {:.2}]", kmin, kmax, min, max);
    for row in grid {
        let line: String = row.into_iter().collect();
        println!("{}", line);
    }
}

fn value_to_row(value: f64, min: f64, max: f64, height: usize) -> usize {
    let t = (value - min) / (max - min);
    let y = (1.0 - t) * (height as f64 - 1.0);
    y.round().clamp(0.0, (height - 1) as f64) as usize
}

// Open uniform (clamped) knot vector for a cubic spline with uniform spacing.
fn open_uniform_knots() -> [f64; KNOTS] {
    let end = (CONTROL_POINTS - DEGREE) as f64;
    core::array::from_fn(|i| {
        if i <= DEGREE {
            0.0
        } else if i >= KNOTS - DEGREE - 1 {
            end
        } else {
            (i - DEGREE) as f64
        }
    })
}

// Solve for control points P so the cubic B-spline interpolates the samples.
// For uniform knots, interior samples satisfy: P_{i-1} + 4 P_i + P_{i+1} = 6 D_i.
// Endpoints are clamped: P_0 = D_0, P_{n-1} = D_{n-1}.
fn interpolate_cubic_uniform(samples: &[f64; CONTROL_POINTS]) -> [f64; CONTROL_POINTS] {
    let n = samples.len();
    let mut a = [[0.0; CONTROL_POINTS + 1]; CONTROL_POINTS];

    // Endpoints are clamped.
    a[0][0] = 1.0;
    a[0][CONTROL_POINTS] = samples[0];
    a[n - 1][n - 1] = 1.0;
    a[n - 1][CONTROL_POINTS] = samples[n - 1];

    // Interior: P_{i-1} + 4 P_i + P_{i+1} = 6 D_i.
    for i in 1..n - 1 {
        a[i][i - 1] = 1.0;
        a[i][i] = 4.0;
        a[i][i + 1] = 1.0;
        a[i][CONTROL_POINTS] = 6.0 * samples[i];
    }

    gaussian_elimination(&mut a)
}

// Solve the dense system in-place with Gauss-Jordan elimination.
fn gaussian_elimination(matrix: &mut [[f64; CONTROL_POINTS + 1]; CONTROL_POINTS]) -> [f64; CONTROL_POINTS] {
    for i in 0..CONTROL_POINTS {
        let mut pivot = i;
        for r in i + 1..CONTROL_POINTS {
            if matrix[r][i].abs() > matrix[pivot][i].abs() {
                pivot = r;
            }
        }
        // If the pivot is tiny, the system is ill-conditioned; keep going with the current row.
        if matrix[pivot][i].abs() < 1e-12 {
            continue;
        }
        if pivot != i {
            matrix.swap(i, pivot);
        }
        let inv = 1.0 / matrix[i][i];
        for c in i..=CONTROL_POINTS {
            matrix[i][c] *= inv;
        }
        for r in 0..CONTROL_POINTS {
            if r == i {
                continue;
            }
            let factor = matrix[r][i];
            for c in i..=CONTROL_POINTS {
                matrix[r][c] -= factor * matrix[i][c];
            }
        }
    }

    core::array::from_fn(|i| matrix[i][CONTROL_POINTS])
}
