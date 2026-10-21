#![allow(incomplete_features)]
#![feature(generic_const_exprs)]

//! Measure oriented-box accuracy and runtime for Bezier and B-spline curves.
//! Run in release mode. `--extended` adds cases to the fixed 2D, 3D and 4D corpus.
//! Timings include model construction. Polynomial extrema check every fitted face.
//! A separate optimizer searches rotations for best-known reference volumes,
//! including starts from the fitted frames. It does not certify global minima.
//! Spline references reuse span extraction, but compute extrema independently.

#[path = "support/oriented_reference.rs"]
mod reference;

use core::{array, fmt::Write, hint::black_box};
use reference::{dot, identity, random, rotated};
use std::time::Instant;
use stroke::{ArrayPoint, BSpline, Bezier, OrientedBoundingBox, OrientedBoundingBoxError, Point};

#[derive(Default)]
struct Summary {
    ratios: [Vec<f64>; 2],
    times: [f64; 2],
}

impl Summary {
    fn report(&mut self, kind: &str, dimension: usize, degree: usize) {
        for (b, budget) in [256, 1024].into_iter().enumerate() {
            let r = &mut self.ratios[b];
            r.sort_by(f64::total_cmp);
            let n = r.len();
            println!(
                "{kind} {dimension}D degree {degree} / budget {budget}: n={n}, median excess={:.6}%, p95={:.6}%, worst={:.6}%, mean total time={:.1}us",
                100.0 * (r[n / 2] - 1.0),
                100.0 * (r[(n * 95 / 100).min(n - 1)] - 1.0),
                100.0 * (r[n - 1] - 1.0),
                self.times[b] / n as f64,
            );
        }
    }
}

struct Case<'a, const D: usize, const N: usize> {
    kind: &'static str,
    knots: &'static str,
    index: usize,
    aspect: f64,
    pieces: &'a [[[f64; D]; N]],
}

fn measure<const D: usize, const N: usize, const BOX_D: usize>(
    case: Case<'_, D, N>,
    fit: impl Fn(usize) -> Result<OrientedBoundingBox<f64, BOX_D>, OrientedBoundingBoxError>,
    summary: &mut Summary,
    csv: &mut String,
) -> Result<(), Box<dyn core::error::Error>> {
    assert_eq!(D, BOX_D);
    let log_volume = |axes: reference::Frame<D>| {
        if case.pieces.len() == 1 {
            return reference::log_volume(&case.pieces[0], axes);
        }
        axes.iter()
            .map(|axis| {
                let (lo, hi) =
                    case.pieces
                        .iter()
                        .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), p| {
                            let (a, b) = reference::range(p.map(|p| dot(axis, &p)));
                            (lo.min(a), hi.max(b))
                        });
                (hi - lo).ln()
            })
            .sum()
    };
    let mut numerical = if case.pieces.len() == 1 {
        reference::reference(&case.pieces[0])
    } else {
        reference::reference_for(log_volume)
    };
    let mut records = Vec::new();
    let mut previous = f64::INFINITY;
    for (b, budget) in [256, 1024].into_iter().enumerate() {
        let start = Instant::now();
        let bounds = fit(budget)?;
        let micros = start.elapsed().as_secs_f64() * 1e6;
        assert!(bounds.evaluations() <= budget);
        assert!(
            bounds.log_volume() <= previous,
            "more work worsened the incumbent"
        );
        previous = bounds.log_volume();
        let axes = array::from_fn(|i| array::from_fn(|j| bounds.axes()[i][j]));
        numerical = numerical.min(if case.pieces.len() == 1 {
            reference::reference_at(&case.pieces[0], axes)
        } else {
            reference::reference_at_for(log_volume, axes)
        });
        let center = bounds.center();
        let half = bounds.half_extents();
        for piece in case.pieces {
            for (i, axis) in axes.iter().enumerate() {
                let (lo, hi) = reference::range(
                    piece.map(|p| (0..D).map(|d| axis[d] * (p[d] - center[d])).sum()),
                );
                assert!(
                    lo >= -half[i] && hi <= half[i],
                    "{} {D}D case {}: face {i} escaped",
                    case.kind,
                    case.index
                );
            }
        }
        records.push((b, budget, micros, bounds));
    }
    for (b, budget, micros, bounds) in records {
        let ratio = (bounds.log_volume() - numerical).exp();
        summary.ratios[b].push(ratio);
        summary.times[b] += micros;
        let info = bounds.diagnostics();
        writeln!(
            csv,
            "{},{D},{},{},{},{},{budget},{},{micros:.3},{ratio:.12},{:.17},{numerical:.17},{},{},{},{:.12}",
            case.kind,
            N - 1,
            case.index,
            case.knots,
            case.aspect,
            bounds.evaluations(),
            bounds.log_volume(),
            info.map_or(0, |d| d.representative_points),
            info.is_some_and(|d| d.representative_limit_reached),
            info.map_or(0, |d| d.point_set_evaluations),
            info.map_or(0.0, |d| d.model_log_volume_gap)
        )?;
    }
    Ok(())
}

fn controls<const D: usize, const N: usize>(state: &mut u64, aspect: f64) -> [[f64; D]; N] {
    let mut axes = identity();
    for i in 0..D {
        for j in i + 1..D {
            axes = rotated(axes, i, j, random(state) * core::f64::consts::PI);
        }
    }
    array::from_fn(|_| {
        let p = array::from_fn(|d| random(state) * if d == 0 { 1.0 } else { aspect });
        axes.map(|axis| dot(&axis, &p))
    })
}

fn beziers<const D: usize, const N: usize>(
    count: usize,
    csv: &mut String,
) -> Result<(), Box<dyn core::error::Error>>
where
    [(); <ArrayPoint<f64, D> as Point>::DIM]: Sized,
{
    let mut state = 0x0bb_5eed + D as u64;
    let mut summary = Summary::default();
    for index in 0..count {
        let aspect = [1.0, 0.3, 0.03, 0.003][index % 4];
        let points = controls::<D, N>(&mut state, aspect);
        let curve = Bezier::new(points.map(ArrayPoint::from))?;
        measure(
            Case {
                kind: "bezier",
                knots: "none",
                index,
                aspect,
                pieces: &[points],
            },
            |budget| black_box(&curve).bounding_box_oriented_approx_with_budget(budget),
            &mut summary,
            csv,
        )?;
    }
    summary.report("bezier", D, N - 1);
    Ok(())
}

fn bsplines<const D: usize, const K: usize, const C: usize, const DEG: usize>(
    count: usize,
    csv: &mut String,
) -> Result<(), Box<dyn core::error::Error>>
where
    [(); <ArrayPoint<f64, D> as Point>::DIM]: Sized,
    [(); DEG + 1]: Sized,
    [(); C + 64]: Sized,
{
    let mut state = 0xb51_5eed + D as u64;
    let mut summary = Summary::default();
    for index in 0..count {
        let aspect = [1.0, 0.3, 0.03, 0.003][(index / 4) % 4];
        let points = controls::<D, C>(&mut state, aspect).map(ArrayPoint::from);
        let uniform = BSpline::<_, K, C, DEG>::new_clamped_uniform(points)?;
        let mut knots = array::from_fn(|i| uniform.knots().as_slice()[i]);
        let kind = match index % 4 {
            0 => "clamped",
            1 => {
                knots = array::from_fn(|i| i as f64);
                "unclamped"
            }
            2 => {
                knots[DEG + 2] = knots[DEG + 1];
                "repeated"
            }
            _ => {
                assert_eq!(C % (DEG + 1), 0);
                let blocks = C / (DEG + 1);
                knots = array::from_fn(|i| (i / (DEG + 1)) as f64 / blocks as f64);
                "discontinuous"
            }
        };
        let curve = BSpline::<_, K, C, DEG>::new(knots, points)?;
        let pieces: Vec<[[f64; D]; DEG + 1]> = curve
            .bezier_spans()
            .map(|(_, span)| span.into_control_points().map(ArrayPoint::into_array))
            .collect();
        measure(
            Case {
                kind: "bspline",
                knots: kind,
                index,
                aspect,
                pieces: &pieces,
            },
            |budget| black_box(&curve).bounding_box_oriented_approx_with_budget(budget),
            &mut summary,
            csv,
        )?;
    }
    summary.report("bspline", D, DEG);
    Ok(())
}

fn main() -> Result<(), Box<dyn core::error::Error>> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    let extended = match args.as_slice() {
        [] => false,
        [arg] if arg == "--extended" => true,
        _ => return Err("usage: oriented_box_study [--extended]".into()),
    };
    let multiplier = if extended { 2 } else { 1 };
    let mut csv = String::from(
        "curve,dimension,degree,case,knot_kind,aspect,budget,full_evaluations,total_us,volume_over_reference,log_volume,reference_log_volume,representative_points,representative_limit_reached,point_set_evaluations,model_log_volume_gap\n",
    );
    beziers::<2, 4>(64 * multiplier, &mut csv)?;
    beziers::<3, 4>(64 * multiplier, &mut csv)?;
    // A cubic lies in an affine subspace of dimension at most three. Use
    // quartics in 4D so a zero-volume optimum cannot distort volume ratios.
    beziers::<4, 5>(16 * multiplier, &mut csv)?;
    let spline_count = if extended { 64 } else { 16 };
    bsplines::<2, 12, 8, 3>(spline_count, &mut csv)?;
    bsplines::<3, 12, 8, 3>(spline_count, &mut csv)?;
    bsplines::<4, 15, 10, 4>(spline_count, &mut csv)?;
    std::fs::create_dir_all("target")?;
    let path = if extended {
        "target/oriented_box_study_extended.csv"
    } else {
        "target/oriented_box_study.csv"
    };
    std::fs::write(path, csv)?;
    println!("Wrote {path}; references are numerical, not certified global optima.");
    Ok(())
}
