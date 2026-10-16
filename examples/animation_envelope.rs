//! A scalar animation/control envelope whose native parameter is time in seconds.
//! Interior control values shape the envelope; they are not interpolated keyframes.

use stroke::{BSpline, PointN};

fn main() -> Result<(), Box<dyn core::error::Error>> {
    let envelope = BSpline::<PointN<f64, 1>, 9, 5, 3>::new(
        [0.0, 0.0, 0.0, 0.0, 2.0, 5.0, 5.0, 5.0, 5.0],
        [0.0, 0.25, 1.2, 0.6, 0.0].map(|value| PointN::from([value])),
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
