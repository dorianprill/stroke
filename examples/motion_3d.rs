//! Sample a cubic 3D tool trajectory at 0.5 coordinate units per second.
//! Convert elapsed time to distance, then invert arc length to find the curve parameter.
//! CSV output gives position, velocity and curvature at 21 equally spaced times.
//! Assume instantaneous starts and stops, with no acceleration or jerk limits.

use stroke::{ArrayPoint, BSpline};

fn main() -> Result<(), Box<dyn core::error::Error>> {
    type P = ArrayPoint<f64, 3>;
    let curve = BSpline::<P, 10, 6, 3>::new_clamped_uniform([
        P::from([0.0, 0.0, 0.0]),
        P::from([1.0, 0.0, 0.5]),
        P::from([2.0, 2.0, 1.0]),
        P::from([3.0, 2.0, 2.0]),
        P::from([4.0, 0.0, 2.5]),
        P::from([5.0, 0.0, 3.0]),
    ])?;
    let accuracy = 1e-9;
    let feed_rate = 0.5; // coordinate units per second
    let length = curve.arclen(accuracy)?;
    let duration = length / feed_rate;

    println!("seconds,x,y,z,vx,vy,vz,curvature");
    for i in 0..=20 {
        let seconds = duration * (i as f64 / 20.0);
        let distance = seconds * feed_rate;
        let t = curve.t_at_length(distance, accuracy)?;
        let point = curve.point_at_length(distance, accuracy)?;
        // With s = feed_rate * seconds, velocity = (dP/ds) * feed_rate.
        // This curve has nonzero speed throughout, so dP/ds is its unit tangent.
        let velocity = curve.tangent(t)? * feed_rate;
        println!(
            "{seconds:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6}",
            point[0],
            point[1],
            point[2],
            velocity[0],
            velocity[1],
            velocity[2],
            curve.curvature(t)?
        );
    }
    Ok(())
}
