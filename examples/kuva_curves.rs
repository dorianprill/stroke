#![allow(incomplete_features)]
#![feature(generic_const_exprs)]

//! Plot four curves with their control points and axis-aligned boxes using Kuva.
//! Solid green boxes bound the curves. Dashed gray boxes bound the controls.
//!
//! Run `cargo run --example kuva_curves` to write `curves.svg` in the current
//! directory. No system font libraries are needed.

use kuva::prelude::{
    Figure, Layout, LegendEntry, LegendPosition, LegendShape, LinePlot, MarkerShape, Plot,
    ScatterPlot, SvgBackend,
};
use stroke::{ArrayPoint, BSpline, Bezier};

type Point = ArrayPoint<f64, 2>;
type Bounds = [(f64, f64); 2];

const OUTPUT: &str = "curves.svg";
const STEPS: usize = 1000;
const CURVE_COLOR: &str = "#d44336";
const CONTROL_COLOR: &str = "#3267b1";
const BOUNDS_COLOR: &str = "#26945d";
const CONTROL_BOUNDS_COLOR: &str = "#777777";
const KNOT_COLOR: &str = "#bc6b00";

fn coordinates(point: Point) -> (f64, f64) {
    (point[0], point[1])
}

fn box_plot(bounds: Bounds) -> LinePlot {
    let [(xmin, xmax), (ymin, ymax)] = bounds;
    LinePlot::new().with_data([
        (xmin, ymin),
        (xmin, ymax),
        (xmax, ymax),
        (xmax, ymin),
        (xmin, ymin),
    ])
}

fn panel(
    title: &str,
    degree: usize,
    control_points: &[Point],
    samples: Vec<Point>,
    bounds: Bounds,
    control_bounds: Bounds,
) -> (Vec<Plot>, Layout) {
    // Allow a few representable steps of roundoff at these coordinate scales.
    let tolerance = 128.0
        * f64::EPSILON
        * control_bounds.iter().fold(1.0_f64, |scale, &(lo, hi)| {
            scale.max(lo.abs()).max(hi.abs())
        });
    for axis in 0..2 {
        assert!(bounds[axis].0 >= control_bounds[axis].0 - tolerance);
        assert!(bounds[axis].1 <= control_bounds[axis].1 + tolerance);
        for point in control_points {
            assert!(point[axis] >= control_bounds[axis].0 - tolerance);
            assert!(point[axis] <= control_bounds[axis].1 + tolerance);
        }
    }
    for point in &samples {
        for axis in 0..2 {
            assert!(point[axis].is_finite());
            assert!(point[axis] >= bounds[axis].0 - tolerance);
            assert!(point[axis] <= bounds[axis].1 + tolerance);
        }
    }

    let plots: Vec<Plot> = vec![
        box_plot(control_bounds)
            .with_color(CONTROL_BOUNDS_COLOR)
            .with_stroke_width(1.5)
            .with_dashed()
            .with_legend("Control box · control_box()")
            .into(),
        box_plot(bounds)
            .with_color(BOUNDS_COLOR)
            .with_stroke_width(2.0)
            .with_legend("Tight box · bounding_box()")
            .into(),
        LinePlot::new()
            .with_data(control_points.iter().copied().map(coordinates))
            .with_color(CONTROL_COLOR)
            .with_stroke_width(1.0)
            .with_dashed()
            .with_legend("Control polygon")
            .into(),
        LinePlot::new()
            .with_data(samples.into_iter().map(coordinates))
            .with_color(CURVE_COLOR)
            .with_stroke_width(2.5)
            .with_legend("Curve")
            .into(),
        ScatterPlot::new()
            .with_data(control_points.iter().copied().map(coordinates))
            .with_color(CONTROL_COLOR)
            .with_size(4.5)
            .with_labels(
                control_points
                    .iter()
                    .enumerate()
                    .map(|(i, p)| format!("P{i} ({}, {})", p[0], p[1])),
            )
            .with_repel_labels()
            .with_legend("Control points")
            .into(),
    ];

    // Include every control point, even when it lies outside the curve's box.
    // Leave room for labels and preserve equal scales on the geometric axes.
    let extent = control_points.iter().fold(bounds, |mut extent, point| {
        for axis in 0..2 {
            extent[axis].0 = extent[axis].0.min(point[axis]);
            extent[axis].1 = extent[axis].1.max(point[axis]);
        }
        extent
    });
    let padding = 0.2 * (extent[0].1 - extent[0].0).max(extent[1].1 - extent[1].0);
    let mut layout = Layout::auto_from_plots(&plots)
        .with_title(title)
        .with_title_size(18)
        .with_subtitle(format!("Degree {degree}"))
        .with_subtitle_size(14)
        .with_x_label("x")
        .with_y_label("y")
        .with_x_axis_min(extent[0].0 - padding)
        .with_x_axis_max(extent[0].1 + padding)
        .with_y_axis_min(extent[1].0 - padding)
        .with_y_axis_max(extent[1].1 + padding)
        .with_equal_aspect()
        .with_show_grid(true);
    // Common series use the figure legend; spline panels add their own knot legends.
    layout.show_legend = false;

    (plots, layout)
}

fn with_knot_legend(
    layout: Layout,
    knots: impl IntoIterator<Item = f64>,
    domain: (f64, f64),
) -> Layout {
    let knots = knots
        .into_iter()
        .map(|knot| knot.to_string())
        .collect::<Vec<_>>()
        .join(", ");
    layout
        .with_legend_title(format!("Knots: [{knots}]"))
        .with_legend_entries(vec![LegendEntry {
            label: format!("Parameter domain: [{}, {}]", domain.0, domain.1),
            color: CURVE_COLOR.into(),
            shape: LegendShape::Line,
            dasharray: None,
        }])
        .with_legend_position(LegendPosition::OutsideBottomCenter)
        .with_legend_width(360.0)
}

fn main() -> Result<(), Box<dyn core::error::Error>> {
    // All Bezier degrees use Bezier<P, N>, including the quadratic and cubic.
    let quadratic = Bezier::new([
        Point::new([0.0, 0.0]),
        Point::new([2.0, 4.0]),
        Point::new([4.0, 0.0]),
    ])
    .map_err(|error| error.to_string())?;
    let cubic = Bezier::new([
        Point::new([0.0, 1.77]),
        Point::new([1.1, -1.0]),
        Point::new([4.3, 3.0]),
        Point::new([3.2, -4.0]),
    ])
    .map_err(|error| error.to_string())?;
    // A degree-5 spline with six copies of each endpoint knot is exactly the
    // original quintic Bezier when it uses the same six control points.
    let clamped = BSpline::<Point, 12, 6, 5>::new(
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
        [
            Point::new([0.0, 0.0]),
            Point::new([0.5, 4.0]),
            Point::new([2.0, -3.0]),
            Point::new([3.0, 4.0]),
            Point::new([4.5, -2.0]),
            Point::new([5.0, 1.0]),
        ],
    )
    .map_err(|error| error.to_string())?;
    // Strictly increasing knots leave the spline unclamped. Its endpoints
    // need not coincide with the first or last control point.
    let spline = BSpline::<Point, 10, 6, 3>::new(
        [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0],
        [
            Point::new([0.0, 0.0]),
            Point::new([1.0, 3.0]),
            Point::new([2.0, -2.0]),
            Point::new([3.0, 3.0]),
            Point::new([4.0, -1.0]),
            Point::new([5.0, 1.0]),
        ],
    )
    .map_err(|error| error.to_string())?;

    // Include both endpoints. B-splines are sampled over their knot domain,
    // which need not be [0, 1]. Propagate any evaluation errors to the caller.
    let (start, end) = clamped.knot_domain();
    let clamped_samples = (0..=STEPS)
        .map(|i| clamped.eval(start + (end - start) * i as f64 / STEPS as f64))
        .collect::<Result<Vec<_>, _>>()
        .map_err(|error| error.to_string())?;
    let (clamped_plots, clamped_layout) = panel(
        "BSpline<ArrayPoint<f64, 2>, 12, 6, 5>",
        5,
        &clamped.control_points().copied().collect::<Vec<_>>(),
        clamped_samples,
        clamped.bounding_box().into_axes(),
        clamped.control_box().into_axes(),
    );
    let clamped_layout = with_knot_legend(clamped_layout, clamped.knots().copied(), (start, end))
        .with_subtitle("Degree 5 · clamped");
    let (start, end) = spline.knot_domain();
    let spline_samples = (0..=STEPS)
        .map(|i| spline.eval(start + (end - start) * i as f64 / STEPS as f64))
        .collect::<Result<Vec<_>, _>>()
        .map_err(|error| error.to_string())?;
    let spline_endpoints = [
        coordinates(spline_samples[0]),
        coordinates(spline_samples[STEPS]),
    ];
    let (mut spline_plots, spline_layout) = panel(
        "BSpline<ArrayPoint<f64, 2>, 10, 6, 3>",
        3,
        &spline.control_points().copied().collect::<Vec<_>>(),
        spline_samples,
        spline.bounding_box().into_axes(),
        spline.control_box().into_axes(),
    );
    // Mark the actual endpoints separately from the blue control points.
    spline_plots.push(
        ScatterPlot::new()
            .with_data(spline_endpoints)
            .with_color(CURVE_COLOR)
            .with_marker(MarkerShape::Diamond)
            .with_size(4.5)
            .with_labels(["Start", "End"])
            .with_repel_labels()
            .with_legend("Unclamped endpoints")
            .into(),
    );
    // Knots are parameter values. Only interior knots in the valid domain
    // get their own markers; domain endpoints are already marked above.
    let interior_knots: Vec<_> = spline
        .knots()
        .copied()
        .filter(|&knot| knot > start && knot < end)
        .collect();
    let knot_points = interior_knots
        .iter()
        .map(|&knot| spline.eval(knot).map(coordinates))
        .collect::<Result<Vec<_>, _>>()
        .map_err(|error| error.to_string())?;
    spline_plots.push(
        ScatterPlot::new()
            .with_data(knot_points)
            .with_color(KNOT_COLOR)
            .with_marker(MarkerShape::Square)
            .with_size(4.5)
            .with_labels(interior_knots.iter().map(|knot| format!("t = {knot}")))
            .with_repel_labels()
            .with_legend("Interior knots")
            .into(),
    );
    let spline_layout = with_knot_legend(spline_layout, spline.knots().copied(), (start, end))
        .with_subtitle("Degree 3 · unclamped");
    let panels = vec![
        panel(
            "Bezier<ArrayPoint<f64, 2>, 3>",
            2,
            &quadratic.into_control_points(),
            (0..=STEPS)
                .map(|i| quadratic.eval(i as f64 / STEPS as f64))
                .collect(),
            quadratic.bounding_box().into_axes(),
            quadratic.control_box().into_axes(),
        ),
        panel(
            "Bezier<ArrayPoint<f64, 2>, 4>",
            3,
            &cubic.into_control_points(),
            (0..=STEPS)
                .map(|i| cubic.eval(i as f64 / STEPS as f64))
                .collect(),
            cubic.bounding_box().into_axes(),
            cubic.control_box().into_axes(),
        ),
        (clamped_plots, clamped_layout),
        (spline_plots, spline_layout),
    ];
    let (plots, layouts) = panels.into_iter().unzip();
    let scene = Figure::new(2, 2)
        .with_title("Const-generic Bezier and B-spline curves")
        .with_plots(plots)
        .with_layouts(layouts)
        .with_cell_size(560.0, 440.0)
        .with_spacing(24.0)
        .with_padding(24.0)
        .with_labels()
        .with_shared_legend_bottom()
        .with_keep_panel_legends()
        .render();

    std::fs::write(OUTPUT, SvgBackend.render_scene(&scene))?;
    println!(
        "Wrote {OUTPUT} (2x2 panels, {} samples per curve)",
        STEPS + 1
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn quadratic_panel_plots_distinct_tight_and_control_boxes() {
        let curve = Bezier::new([
            Point::new([0.0, 0.0]),
            Point::new([2.0, 4.0]),
            Point::new([4.0, 0.0]),
        ])
        .unwrap();
        let tight = curve.bounding_box().into_axes();
        let control = curve.control_box().into_axes();
        assert_eq!(tight, [(0.0, 4.0), (0.0, 2.0)]);
        assert_eq!(control, [(0.0, 4.0), (0.0, 4.0)]);
        let (plots, _) = panel(
            "Quadratic",
            2,
            &curve.into_control_points(),
            (0..=STEPS)
                .map(|i| curve.eval(i as f64 / STEPS as f64))
                .collect(),
            tight,
            control,
        );
        for (plot, bounds, color, dashed, label) in [
            (
                &plots[0],
                control,
                CONTROL_BOUNDS_COLOR,
                true,
                "Control box · control_box()",
            ),
            (
                &plots[1],
                tight,
                BOUNDS_COLOR,
                false,
                "Tight box · bounding_box()",
            ),
        ] {
            let Plot::Line(line) = plot else {
                panic!("expected a box outline");
            };
            let expected = box_plot(bounds);
            assert_eq!(
                line.data.iter().map(|p| (p.x, p.y)).collect::<Vec<_>>(),
                expected.data.iter().map(|p| (p.x, p.y)).collect::<Vec<_>>()
            );
            assert_eq!(line.color, color);
            assert_eq!(line.line_style.dasharray().is_some(), dashed);
            assert_eq!(line.legend_label.as_deref(), Some(label));
        }
    }

    #[test]
    fn box_outline_preserves_degenerate_bounds() {
        let line = box_plot([(2.0, 2.0), (-1.0, 3.0)]);
        assert_eq!(
            line.data.iter().map(|p| (p.x, p.y)).collect::<Vec<_>>(),
            [
                (2.0, -1.0),
                (2.0, 3.0),
                (2.0, 3.0),
                (2.0, -1.0),
                (2.0, -1.0)
            ]
        );
    }
}
