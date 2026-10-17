//! Sample a scalar cubic B-spline envelope over five seconds at 0.25 s intervals.
//! Fourfold endpoint knots fix the initial and final values at zero. The simple
//! knot at 2 s joins two polynomial pieces with continuous first/second derivatives.
//! CSV reports the envelope value and its native-time derivative in units/s.

use stroke::{ArrayPoint, BSpline};

fn main() -> Result<(), Box<dyn core::error::Error>> {
    let envelope = BSpline::<ArrayPoint<f64, 1>, 9, 5, 3>::new(
        [0.0, 0.0, 0.0, 0.0, 2.0, 5.0, 5.0, 5.0, 5.0],
        [0.0, 0.25, 1.2, 0.6, 0.0].map(|value| ArrayPoint::from([value])),
    )?;
    println!("seconds,value,rate_per_second");
    for i in 0..=20 {
        let seconds = i as f64 * 0.25;
        let value = envelope.eval(seconds)?[0];
        let rate = envelope.derivative_at(seconds)?[0];
        println!("{seconds:.2},{value:.6},{rate:.6}");
    }
    Ok(())
}
