#[cfg(feature = "nalgebra")]
fn main() {
    use nalgebra::SVector;
    use stroke::Bezier;

    let curve = Bezier::<SVector<f32, 2>, 3>::new([
        SVector::<f32, 2>::new(0.0, 0.0),
        SVector::<f32, 2>::new(1.0, 0.0),
        SVector::<f32, 2>::new(1.0, 1.0),
    ]);

    let mid = curve.eval(0.5);
    let bounds = curve.bounding_box();
    println!("mid={:?}, bounds={:?}", mid, bounds);
}

#[cfg(not(feature = "nalgebra"))]
fn main() {
    eprintln!("Enable the nalgebra feature: cargo run --example nalgebra_basic --features nalgebra");
}
