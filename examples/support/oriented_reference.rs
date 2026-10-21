//! Polynomial extrema and numerical reference volumes for tests and the study.
//! Use modest coordinates. This reference code is outside the library.
use core::array;

pub type Frame<const D: usize> = [[f64; D]; D];

pub fn identity<const D: usize>() -> Frame<D> {
    array::from_fn(|i| array::from_fn(|j| if i == j { 1.0 } else { 0.0 }))
}
pub fn dot<const D: usize>(a: &[f64; D], b: &[f64; D]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}
pub fn rotated<const D: usize>(mut axes: Frame<D>, i: usize, j: usize, angle: f64) -> Frame<D> {
    let (s, c) = angle.sin_cos();
    let a = axes[i];
    let b = axes[j];
    axes[i] = array::from_fn(|d| c * a[d] + s * b[d]);
    axes[j] = array::from_fn(|d| -s * a[d] + c * b[d]);
    axes
}
pub fn random(state: &mut u64) -> f64 {
    *state = state
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    ((*state >> 11) as f64 / ((1_u64 << 53) as f64)) * 2.0 - 1.0
}
fn quadratic(a: f64, b: f64, c: f64) -> Vec<f64> {
    if a == 0.0 {
        return if b == 0.0 { vec![] } else { vec![-c / b] };
    }
    let d = b * b - 4.0 * a * c;
    if d < 0.0 {
        return vec![];
    }
    let q = -0.5 * (b + d.sqrt().copysign(b));
    if q == 0.0 {
        vec![-b / (2.0 * a)]
    } else {
        vec![q / a, c / q]
    }
}
fn horner(p: [f64; 4], t: f64) -> f64 {
    ((p[3] * t + p[2]) * t + p[1]) * t + p[0]
}
fn cubic_roots(p: [f64; 4]) -> Vec<f64> {
    if p[3] == 0.0 {
        return quadratic(p[2], p[1], p[0]);
    }
    // Partition at derivative roots. Each open interval is monotone; isolate
    // its root by bisection. Critical points themselves include repeated roots.
    let mut partitions = vec![0.0, 1.0];
    partitions.extend(
        quadratic(3.0 * p[3], 2.0 * p[2], p[1])
            .into_iter()
            .filter(|&t| t > 0.0 && t < 1.0),
    );
    partitions.sort_by(f64::total_cmp);
    let mut roots = partitions.clone();
    for interval in partitions.windows(2) {
        let mut lo = interval[0];
        let mut hi = interval[1];
        let sign = horner(p, lo).is_sign_negative();
        if sign == horner(p, hi).is_sign_negative() {
            continue;
        }
        for _ in 0..56 {
            let mid = (lo + hi) * 0.5;
            if horner(p, mid).is_sign_negative() == sign {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        roots.push((lo + hi) * 0.5);
    }
    roots
}
pub fn range<const N: usize>(p: [f64; N]) -> (f64, f64) {
    let mut lo = p[0].min(p[N - 1]);
    let mut hi = p[0].max(p[N - 1]);
    let roots = match N {
        1 | 2 => vec![],
        3 => quadratic(0.0, 2.0 * (p[0] - 2.0 * p[1] + p[2]), 2.0 * (p[1] - p[0])),
        4 => quadratic(
            3.0 * (-p[0] + 3.0 * p[1] - 3.0 * p[2] + p[3]),
            6.0 * (p[0] - 2.0 * p[1] + p[2]),
            3.0 * (p[1] - p[0]),
        ),
        5 => cubic_roots([
            4.0 * (p[1] - p[0]),
            12.0 * (p[0] - 2.0 * p[1] + p[2]),
            12.0 * (-p[0] + 3.0 * p[1] - 3.0 * p[2] + p[3]),
            4.0 * (p[0] - 4.0 * p[1] + 6.0 * p[2] - 4.0 * p[3] + p[4]),
        ]),
        _ => panic!("reference supports degrees zero through four"),
    };
    for t in roots.into_iter().filter(|&t| t > 0.0 && t < 1.0) {
        let mut q = p;
        for level in 1..N {
            for i in 0..N - level {
                q[i] = (1.0 - t) * q[i] + t * q[i + 1];
            }
        }
        lo = lo.min(q[0]);
        hi = hi.max(q[0]);
    }
    (lo, hi)
}
pub fn log_volume<const D: usize, const N: usize>(controls: &[[f64; D]; N], axes: Frame<D>) -> f64 {
    axes.map(|axis| {
        let (lo, hi) = range(controls.map(|p| dot(&axis, &p)));
        (hi - lo).ln()
    })
    .iter()
    .sum()
}
fn frame<const D: usize>(angles: &[f64]) -> Frame<D> {
    let mut axes = identity();
    let mut index = 0;
    for i in 0..D {
        for j in i + 1..D {
            axes = rotated(axes, i, j, angles[index]);
            index += 1;
        }
    }
    axes
}
// Nelder-Mead reference search over joint rotation angles. Only tests and
// the study use this optimizer. The library fitter has its own rotation search.
fn simplex(f: impl Fn(&[f64]) -> f64, initial: &[f64], budget: usize) -> f64 {
    let m = initial.len();
    let mut vertices = vec![initial.to_vec(); m + 1];
    for i in 0..m {
        vertices[i + 1][i] += 0.15;
    }
    let mut scores: Vec<_> = vertices.iter().map(|v| f(v)).collect();
    let mut evaluations = m + 1;
    while evaluations + 2 + m <= budget {
        let mut order: Vec<_> = (0..=m).collect();
        order.sort_by(|&a, &b| scores[a].total_cmp(&scores[b]));
        vertices = order.iter().map(|&i| vertices[i].clone()).collect();
        scores = order.iter().map(|&i| scores[i]).collect();
        let spread = vertices[1..]
            .iter()
            .flat_map(|v| v.iter().zip(&vertices[0]).map(|(a, b)| (a - b).abs()))
            .fold(0.0, f64::max);
        if spread < 1e-10 {
            break;
        }
        let center: Vec<_> = (0..m)
            .map(|d| vertices[..m].iter().map(|v| v[d]).sum::<f64>() / m as f64)
            .collect();
        let reflection: Vec<_> = center
            .iter()
            .zip(&vertices[m])
            .map(|(c, w)| 2.0 * c - w)
            .collect();
        let fr = f(&reflection);
        evaluations += 1;
        if fr < scores[0] {
            let expansion: Vec<_> = center
                .iter()
                .zip(&reflection)
                .map(|(c, r)| c + 2.0 * (r - c))
                .collect();
            let fe = f(&expansion);
            evaluations += 1;
            if fe < fr {
                vertices[m] = expansion;
                scores[m] = fe;
            } else {
                vertices[m] = reflection;
                scores[m] = fr;
            }
        } else if fr < scores[m - 1] {
            vertices[m] = reflection;
            scores[m] = fr;
        } else {
            let outside = fr < scores[m];
            let target = if outside { &reflection } else { &vertices[m] };
            let contraction: Vec<_> = center
                .iter()
                .zip(target)
                .map(|(c, t)| c + 0.5 * (t - c))
                .collect();
            let fc = f(&contraction);
            evaluations += 1;
            if fc < (if outside { fr } else { scores[m] }) {
                vertices[m] = contraction;
                scores[m] = fc;
            } else {
                for i in 1..=m {
                    vertices[i] = (0..m)
                        .map(|d| vertices[0][d] + 0.5 * (vertices[i][d] - vertices[0][d]))
                        .collect();
                    scores[i] = f(&vertices[i]);
                    evaluations += 1;
                }
            }
        }
    }
    scores.into_iter().fold(f64::INFINITY, f64::min)
}
pub fn reference<const D: usize, const N: usize>(controls: &[[f64; D]; N]) -> f64 {
    reference_for(|axes| log_volume(controls, axes))
}

pub fn reference_for<const D: usize>(log_volume: impl Fn(Frame<D>) -> f64) -> f64 {
    let m = D * (D - 1) / 2;
    let objective = |x: &[f64]| log_volume(frame(x));
    if D == 1 {
        return objective(&[]);
    }
    if D == 2 {
        let step = core::f64::consts::FRAC_PI_2 / 4096.0;
        let mut angle = 0.0;
        let mut best = objective(&[angle]);
        for i in 1..4096 {
            let a = i as f64 * step;
            let score = objective(&[a]);
            if score < best {
                best = score;
                angle = a;
            }
        }
        let mut lo = angle - step;
        let mut hi = angle + step;
        for _ in 0..64 {
            let a = lo + (hi - lo) * 0.3819660112501051;
            let b = hi - (hi - lo) * 0.3819660112501051;
            let fa = objective(&[a]);
            let fb = objective(&[b]);
            best = best.min(fa).min(fb);
            if fa < fb {
                hi = b;
            } else {
                lo = a;
            }
        }
        return best;
    }
    let mut starts = vec![vec![0.0; m]];
    if D == 3 {
        let mut grid = Vec::new();
        for a in 0..18 {
            for b in 0..18 {
                for c in 0..18 {
                    let x = vec![a as f64, b as f64, c as f64]
                        .into_iter()
                        .map(|x| x * core::f64::consts::PI / 18.0)
                        .collect::<Vec<_>>();
                    grid.push((objective(&x), x));
                }
            }
        }
        grid.sort_by(|a, b| a.0.total_cmp(&b.0));
        starts.extend(grid.into_iter().take(32).map(|(_, x)| x));
    }
    let mut state = 0xf00d_cafe;
    starts.extend((0..64).map(|_| {
        (0..m)
            .map(|_| random(&mut state) * core::f64::consts::PI)
            .collect::<Vec<_>>()
    }));
    starts
        .iter()
        .map(|x| simplex(objective, x, if D > 3 { 4096 } else { 2048 }))
        .fold(f64::INFINITY, f64::min)
}

// Include the fitted frame among the reference search's starting points.
// Reference extrema and optimization are independent, but these starts come
// from the library fitter.
pub fn reference_at<const D: usize, const N: usize>(
    controls: &[[f64; D]; N],
    base: Frame<D>,
) -> f64 {
    reference_at_for(|axes| log_volume(controls, axes), base)
}

pub fn reference_at_for<const D: usize>(
    log_volume: impl Fn(Frame<D>) -> f64,
    base: Frame<D>,
) -> f64 {
    let m = D * (D - 1) / 2;
    let objective = |angles: &[f64]| {
        let mut axes = base;
        let mut k = 0;
        for i in 0..D {
            for j in i + 1..D {
                axes = rotated(axes, i, j, angles[k]);
                k += 1;
            }
        }
        log_volume(axes)
    };
    let initial = vec![0.0; m];
    if m == 0 {
        return objective(&initial);
    }
    simplex(objective, &initial, 4096).min(objective(&initial))
}
