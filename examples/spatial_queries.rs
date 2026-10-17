//! CAD/GIS-style candidate filtering with axis-aligned boxes.
//! Box overlap is only a broad-phase test, not proof of curve intersection.

use stroke::{ArrayPoint, BoundingBox, QuadraticBezier};

fn main() -> Result<(), Box<dyn core::error::Error>> {
    let curve = QuadraticBezier::new(
        ArrayPoint::from([0.0_f64, 0.0]),
        ArrayPoint::from([2.0, 4.0]),
        ArrayPoint::from([4.0, 0.0]),
    )?;
    let regions = [
        BoundingBox::new([1.0, 3.0], [3.0, 3.5]).unwrap(),
        BoundingBox::new([1.0, 0.5], [3.0, 1.0]).unwrap(),
    ];
    let tight = curve.bounding_box();
    let control = curve.control_box();
    for (i, region) in regions.iter().enumerate() {
        println!(
            "region {i}: control candidate={}, tight candidate={}, overlap={:?}",
            control.intersects(region),
            tight.intersects(region),
            tight.intersection(region)
        );
    }
    assert!(control.intersects(&regions[0]));
    assert!(!tight.intersects(&regions[0]));
    assert!(tight.intersects(&regions[1]));
    let trimmed = curve.subcurve(0.25..0.75)?;
    println!(
        "middle-half box: {:?}, center: {:?}, extent: {:?}",
        trimmed.bounding_box(),
        tight.center(),
        tight.extent()
    );
    Ok(())
}
