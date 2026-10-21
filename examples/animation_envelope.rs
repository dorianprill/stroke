//! Sample a scalar cubic B-spline over five seconds at 0.25-second intervals.
//! Four copies of each endpoint knot make both endpoint values zero. A single
//! interior knot at 2 seconds joins pieces with continuous first and second derivatives.
//! The parameter is time in seconds. CSV output gives the value and its rate per second.

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
