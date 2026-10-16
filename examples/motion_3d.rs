//! Constant-path-speed position sampling for a 3D tool or robot trajectory.
//! This is geometry, not an acceleration/jerk-limited or collision-safe motion planner.

use stroke::{BSpline, PointN};

fn main() -> Result<(), Box<dyn core::error::Error>> {
    type P = PointN<f64, 3>;
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
        // Native parameter speed is not physical speed. The unit tangent times
        // the chosen feed rate gives the physical velocity on this regular curve.
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
