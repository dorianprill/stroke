use stroke::{BSpline, PointN};

const DEGREE: usize = 3;
const CONTROL_POINTS: usize = 8;
const KNOTS: usize = CONTROL_POINTS + DEGREE + 1;

fn main() -> Result<(), Box<dyn core::error::Error>> {
    // Sample a 1D signal.
    let samples: [f64; CONTROL_POINTS] = core::array::from_fn(|i| {
        let x = i as f64 / (CONTROL_POINTS - 1) as f64;
        let base = (2.0 * core::f64::consts::PI * x).sin();
        let ripple = 0.2 * (12.0 * core::f64::consts::PI * x).sin();
        base + ripple
    });

    // Interpolate the samples with a cubic B-spline.
    // We solve a basis-function collocation system to recover control points that force the
    // spline to pass through the sample values at uniform parameter steps.
    let control = interpolate_cubic_uniform(&samples)?;
    let control_points: [PointN<f64, 1>; CONTROL_POINTS] =
        core::array::from_fn(|i| PointN::new([control[i]]));
    // Use a clamped uniform knot vector so the spline starts/ends at the samples.
    let spline: BSpline<PointN<f64, 1>, KNOTS, CONTROL_POINTS, DEGREE> =
        BSpline::new_clamped_uniform(control_points)?;

    // Verify the example actually interpolates every sample, not just the endpoints.
    for (i, expected) in samples.iter().enumerate() {
        let t = i as f64 / (CONTROL_POINTS - 1) as f64;
        assert!((spline.eval(t)?[0] - expected).abs() < 1e-12);
    }

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
    // Overlay controls at their Greville abscissae (not at uniform sample times).
    for (i, value) in control.iter().copied().enumerate() {
        let t = spline.knots().skip(i + 1).take(DEGREE).sum::<f64>() / DEGREE as f64;
        let col = ((t - kmin) / (kmax - kmin) * (width - 1) as f64).round() as usize;
        let row = value_to_row(value, min, max, height);
        let letter = if i < CONTROL_POINTS / 2 { 'a' } else { 'b' };
        grid[row][col] = letter;
    }

    println!("B-spline interpolation of a 1D signal (* = spline, a/b = control points)");
    println!(
        "domain: [{:.2}, {:.2}], range: [{:.2}, {:.2}]",
        kmin, kmax, min, max
    );
    for row in grid {
        let line: String = row.into_iter().collect();
        println!("{}", line);
    }
    Ok(())
}

fn value_to_row(value: f64, min: f64, max: f64, height: usize) -> usize {
    let t = (value - min) / (max - min);
    let y = (1.0 - t) * (height as f64 - 1.0);
    y.round().clamp(0.0, (height - 1) as f64) as usize
}

// Solve for control points P so the cubic B-spline interpolates the samples.
// Build the actual basis matrix: the 1:4:1 cardinal-spline stencil does not
// describe a clamped spline sampled at arbitrary uniform parameter steps.
fn interpolate_cubic_uniform(
    samples: &[f64; CONTROL_POINTS],
) -> Result<[f64; CONTROL_POINTS], Box<dyn core::error::Error>> {
    let basis_curve =
        BSpline::<PointN<f64, 1>, KNOTS, CONTROL_POINTS, DEGREE>::new_clamped_uniform(
            [PointN::new([0.0]); CONTROL_POINTS],
        )?;
    let mut a = [[0.0; CONTROL_POINTS + 1]; CONTROL_POINTS];
    for (i, row) in a.iter_mut().enumerate() {
        let t = i as f64 / (CONTROL_POINTS - 1) as f64;
        let span = basis_curve.knot_span(t)?;
        let basis = basis_curve.basis_functions(span, t)?;
        row[span - DEGREE..=span].copy_from_slice(&basis);
        row[CONTROL_POINTS] = samples[i];
    }
    Ok(gaussian_elimination(&mut a)?)
}

// Solve the dense system in-place with Gauss-Jordan elimination.
fn gaussian_elimination(
    matrix: &mut [[f64; CONTROL_POINTS + 1]; CONTROL_POINTS],
) -> Result<[f64; CONTROL_POINTS], &'static str> {
    if matrix.iter().flatten().any(|value| !value.is_finite()) {
        return Err("interpolation system must be finite");
    }
    let scale = matrix
        .iter()
        .flat_map(|row| &row[..CONTROL_POINTS])
        .fold(0.0_f64, |m, value| m.max(value.abs()));
    let threshold = scale * f64::EPSILON * CONTROL_POINTS as f64;
    for i in 0..CONTROL_POINTS {
        let mut pivot = i;
        for r in i + 1..CONTROL_POINTS {
            if matrix[r][i].abs() > matrix[pivot][i].abs() {
                pivot = r;
            }
        }
        if matrix[pivot][i].abs() <= threshold {
            return Err("interpolation system is singular or ill-conditioned");
        }
        if pivot != i {
            matrix.swap(i, pivot);
        }
        let inv = 1.0 / matrix[i][i];
        for value in &mut matrix[i][i..] {
            *value *= inv;
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

    let result = core::array::from_fn(|i| matrix[i][CONTROL_POINTS]);
    if result.iter().all(|value| value.is_finite()) {
        Ok(result)
    } else {
        Err("interpolation arithmetic overflowed")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reconstructs_samples_with_clamped_basis() {
        for samples in [
            [1.0; CONTROL_POINTS],
            [0.0, 1.0, -0.5, 2.0, 0.25, -1.0, 0.5, 0.0],
        ] {
            let control = interpolate_cubic_uniform(&samples).unwrap();
            let curve =
                BSpline::<PointN<f64, 1>, KNOTS, CONTROL_POINTS, DEGREE>::new_clamped_uniform(
                    control.map(|value| PointN::new([value])),
                )
                .unwrap();
            for (i, expected) in samples.iter().enumerate() {
                assert!(
                    (curve.eval(i as f64 / (CONTROL_POINTS - 1) as f64).unwrap()[0] - expected)
                        .abs()
                        < 1e-12
                );
            }
        }
    }

    #[test]
    fn rejects_invalid_interpolation_systems() {
        assert!(gaussian_elimination(&mut [[0.0; CONTROL_POINTS + 1]; CONTROL_POINTS]).is_err());
        assert!(interpolate_cubic_uniform(&[f64::NAN; CONTROL_POINTS]).is_err());
    }
}
