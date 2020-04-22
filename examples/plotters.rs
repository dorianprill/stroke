extern crate plotters;
use plotters::prelude::*;

extern crate brezel;
use brezel::cubic_bezier::CubicBezier;
use brezel::point2::Point2;
use brezel::point::Point;


fn main() -> Result<(), Box<dyn std::error::Error>> {

    let cp = vec![(0f64,  1.77f64),
                 (1.1f64, -1f64),
                 (4.3f64, 3f64),
                 (3.2f64, -4f64)];
    
    let ch = vec![(0f64,  1.77f64),
                    (1.1f64, -1f64),
                    (3.2f64, -4f64),
                    (4.3f64, 3f64),
                    (0f64,  1.77f64)];

    // TODO point2 needs a constructor for tuples
    // TODO beziers need a convex_hull() fn that returns tuples
    let bezier = CubicBezier::new( 
                Point2::new(0f64,  1.77f64),
                Point2::new(1.1f64, -1f64),
                Point2::new(4.3f64,3f64),
                Point2::new(3.2f64, -4f64),
            );

    let ((xmin, ymin), (xmax, ymax)) = bezier.bounding_box::<f64>();

    let nsteps: usize =  1000;                            
    let mut bezier_path: Vec<(f64, f64)> = Vec::with_capacity(nsteps);       
    for t in 0..nsteps {
        let t = t as f64 * 1f64/(nsteps as f64);
        let p = bezier.eval_casteljau(t);
        bezier_path.push((p.x(), p.y()));
    }

    let root = BitMapBackend::new("cubic_bezier_bounding_box.png", (640, 480)).into_drawing_area();
    root.fill(&WHITE)?;

    // setup the chart
    let mut chart = ChartBuilder::on(&root)
        .caption("Cubic Bezier Curve", ("sans-serif", 21).into_font())
        .margin(5)
        .x_label_area_size(30)
        .y_label_area_size(30)
        .build_ranged((xmin-2.0)..xmin+6.0, ymin-1.0..ymax+3.0)?; // make graph a bit bigger than bounding box 

    chart.configure_mesh().draw()?;

    // draw the control points
    chart
        .draw_series(PointSeries::of_element(
            cp.clone(),
            5,
            &BLUE,
            &|coord, size, style| {
                EmptyElement::at(coord)
                    + Circle::new((0, 0), size, style)
                    + Text::new(
                        format!("{:?}", coord),
                        (0, 15),
                        ("sans-serif", 15).into_font(),
                    )
            },
        ))?
        .label("Control Points")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &BLUE));

    // draw the actual bezier curve
    chart
        .draw_series(LineSeries::new(
            bezier_path,
            &RED,
        ))?
        .label("B(t)")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &RED));

    // draw the bounding box
    chart.draw_series(
        AreaSeries::new(
            vec![(xmin, ymin),
                (xmin, ymax),
                (xmax, ymax),
                (xmax, ymin),
                (xmin, ymin)],
            0.0,
            &GREEN.mix(0.05),
        ).border_style(&GREEN)
    )?
    .label("Bounding Box")
    .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &GREEN));

    // draw the control polygon
    chart
        .draw_series(
            AreaSeries::new(
            ch.clone(),
            0.0,
            &BLUE.mix(0.0),
            ).border_style(&BLUE)
        )?
        .label("Control Polygon")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &BLUE));

    chart
        .configure_series_labels()
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .draw()?;

    Ok(())
}
        