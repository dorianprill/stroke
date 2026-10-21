//! Plot cubic Bernstein and B-spline weights over their parameter domain.
//! Run `cargo run --example kuva_basis` to regenerate `basis.svg`.

use core::array;
use kuva::prelude::{Figure, Layout, LegendPosition, LinePlot, Plot, SvgBackend};
use stroke::{ArrayPoint, BSpline, BSplineError, Bezier, CurveConstructionError};

type Point = ArrayPoint<f64, 1>;
type Series = Vec<(f64, f64)>;

const OUTPUT: &str = "basis.svg";
const STEPS: usize = 400;
const DEGREE: usize = 3;
const KNOTS: [f64; 11] = [0.0, 0.0, 0.0, 0.0, 0.25, 0.5, 0.75, 1.0, 1.0, 1.0, 1.0];
const COLORS: [&str; 7] = [
    "#0072b2", "#e69f00", "#009e73", "#7c3aed", "#d55e00", "#56b4e9", "#cc397b",
];

fn bezier() -> Result<Bezier<Point, 4>, CurveConstructionError> {
    Bezier::new(array::from_fn(|i| Point::from([i as f64])))
}

fn spline() -> Result<BSpline<Point, 11, 7, DEGREE>, BSplineError> {
    BSpline::new(KNOTS, array::from_fn(|i| Point::from([i as f64])))
}

fn sample<const N: usize, E>(
    mut evaluate: impl FnMut(f64) -> Result<[f64; N], E>,
) -> Result<[Series; N], E> {
    let mut series = array::from_fn(|_| Vec::with_capacity(STEPS + 1));
    for step in 0..=STEPS {
        let t = step as f64 / STEPS as f64;
        let weights = evaluate(t)?;
        for (series, weight) in series.iter_mut().zip(weights) {
            series.push((t, weight));
        }
    }
    Ok(series)
}

fn spline_weights(curve: &BSpline<Point, 11, 7, DEGREE>, t: f64) -> Result<[f64; 7], BSplineError> {
    let span = curve.knot_span(t)?;
    let active = curve.basis_functions(span, t)?;
    // The API returns only the four active weights, in global control-point order.
    let mut weights = [0.0; 7];
    weights[span - DEGREE..=span].copy_from_slice(&active);
    Ok(weights)
}

fn panel<const N: usize>(
    title: &str,
    subtitle: &str,
    prefix: &str,
    series: [Series; N],
    interior_knots: &[f64],
) -> (Vec<Plot>, Layout) {
    let mut plots: Vec<Plot> = interior_knots
        .iter()
        .map(|&t| {
            LinePlot::new()
                .with_data([(t, 0.0), (t, 1.0)])
                .with_color("#a0a8b3")
                .with_stroke_width(1.0)
                .with_dashed()
                .into()
        })
        .collect();
    for (i, data) in series.into_iter().enumerate() {
        plots.push(
            LinePlot::new()
                .with_data(data)
                .with_color(COLORS[i])
                .with_stroke_width(2.8)
                .with_legend(format!("{prefix}{i}"))
                .into(),
        );
    }
    let layout = Layout::auto_from_plots(&plots)
        .with_title(title)
        .with_title_size(20)
        .with_subtitle(subtitle)
        .with_subtitle_size(14)
        .with_x_label("Parameter t")
        .with_y_label("Basis weight")
        .with_x_axis_min(-0.015)
        .with_x_axis_max(1.015)
        .with_y_axis_min(-0.025)
        .with_y_axis_max(1.05)
        .with_x_tick_step(0.25)
        .with_y_tick_step(0.25)
        .with_show_grid(true)
        .with_legend_position(LegendPosition::OutsideBottomColumns)
        .with_legend_col_limit(7)
        .with_legend_box(false);
    (plots, layout)
}

fn render() -> Result<String, Box<dyn core::error::Error>> {
    let bezier = bezier()?;
    let spline = spline()?;
    let panels = vec![
        panel(
            "Cubic Bézier basis",
            "4 Bernstein weights",
            "B",
            sample(|t| bezier.basis_functions(t))?,
            &[],
        ),
        panel(
            "Cubic B-spline basis",
            "7 weights · clamped uniform knots",
            "N",
            sample(|t| spline_weights(&spline, t))?,
            &KNOTS[4..7],
        ),
    ];
    let (plots, layouts) = panels.into_iter().unzip();
    let scene = Figure::new(1, 2)
        .with_plots(plots)
        .with_layouts(layouts)
        .with_cell_size(560.0, 360.0)
        .with_spacing(24.0)
        .with_padding(24.0)
        .render();
    Ok(SvgBackend.render_scene(&scene))
}

fn main() -> Result<(), Box<dyn core::error::Error>> {
    std::fs::write(OUTPUT, render()?)?;
    println!(
        "Wrote {OUTPUT} (1x2 panels, {} samples per basis weight)",
        STEPS + 1
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn check_partition<const N: usize>(series: &[Series; N]) {
        assert!(series.iter().all(|series| series.len() == STEPS + 1));
        for step in 0..=STEPS {
            let t = step as f64 / STEPS as f64;
            let mut sum = 0.0;
            for series in series {
                assert_eq!(series[step].0, t);
                let weight = series[step].1;
                assert!(weight.is_finite() && (0.0..=1.0).contains(&weight));
                sum += weight;
            }
            assert!((sum - 1.0).abs() <= 32.0 * f64::EPSILON);
        }
        for (i, series) in series.iter().enumerate() {
            assert_eq!(series[0].1, f64::from(i == 0));
            assert_eq!(series[STEPS].1, f64::from(i == N - 1));
        }
    }

    #[test]
    fn plotted_weights_are_nonnegative_partitions_of_unity() {
        let bezier = bezier().unwrap();
        let spline = spline().unwrap();
        check_partition(&sample(|t| bezier.basis_functions(t)).unwrap());
        check_partition(&sample(|t| spline_weights(&spline, t)).unwrap());
    }

    #[test]
    fn spline_series_keep_global_indices_and_local_support() {
        let spline = spline().unwrap();
        let series = sample(|t| spline_weights(&spline, t)).unwrap();
        for (i, series) in series.iter().enumerate() {
            for &(t, weight) in series {
                if t < KNOTS[i] || t > KNOTS[i + DEGREE + 1] {
                    assert_eq!(weight, 0.0);
                }
            }
        }
        let middle = spline_weights(&spline, 0.5).unwrap();
        for (actual, expected) in
            middle
                .into_iter()
                .zip([0.0, 0.0, 1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0, 0.0, 0.0])
        {
            assert!((actual - expected).abs() <= 8.0 * f64::EPSILON);
        }
    }

    #[test]
    fn figure_renders_both_bases_and_every_legend_entry() {
        let svg = render().unwrap();
        assert!(svg.starts_with("<svg"));
        assert!(svg.contains("Cubic Bézier basis"));
        assert!(svg.contains("Cubic B-spline basis"));
        for label in [
            "B0", "B1", "B2", "B3", "N0", "N1", "N2", "N3", "N4", "N5", "N6",
        ] {
            assert!(svg.contains(&format!(">{label}</text>")));
        }
        assert!(!svg.contains("NaN"));
        assert!(!svg.contains("Infinity"));
    }
}
