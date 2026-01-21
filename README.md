# stroke  

![Rust](https://github.com/dorianprill/stroke/actions/workflows/rust.yml/badge.svg)
![crates.io](https://img.shields.io/crates/v/stroke.svg)

![A Cubic Bézier Curve with Bounding Box and Convex Hull rendered by plotters.rs](https://raw.githubusercontent.com/dorianprill/stroke-rs/main/cubic_bezier_bounding_box.png)  

made with [plotters.rs](https://github.com/38/plotters)  

A zero-allocation library providing const-generic implementations of Bézier curves, B-Spline curves and specialized implementations of up to cubic Bézier curves in N-dimensional euclidean space. It is intended for general/embedded/wasm use supporting #![no_std] environments written in 100% safe Rust with minimal dependencies.  

The library makes heavy use of const-generics and `generic_const_exprs`, so the nightly compiler is required.  
It comes with a const-generic N-dimensional Point type so you can use the library without any other dependencies.  
`PointN<T, N>` uses your chosen scalar `T` (typically `f32`/`f64`) and keeps the dimension at compile time.  
Should you want to integrate with types provided by another library, implement the small `Point` trait. Optional extension traits (`PointIndex`, `PointDot`, `PointNorm`) unlock component access and geometric helpers when needed.  

### Integration features

- `nalgebra`: implements `Point` for `nalgebra::SVector<T, D>` (no_std + libm).
  Requires `T: nalgebra::RealField + num_traits::Float` (typically `f32` or `f64`).

Example:
```rust
use nalgebra::SVector;
use stroke::Bezier;

let curve = Bezier::<SVector<f32, 2>, 3>::new([
    SVector::<f32, 2>::new(0.0, 0.0),
    SVector::<f32, 2>::new(1.0, 0.0),
    SVector::<f32, 2>::new(1.0, 1.0),
]);

let mid = curve.eval(0.5);
```

Enable with:
```
stroke = { version = "0.2.0", features = ["nalgebra"] }
nalgebra = { version = "0.32", default-features = false, features = ["libm"] }
```

Right now, the generic versions don't implement all methods that the specialized versions do (as the algorithms get a bit more complicated) but should reach parity eventually.

## Goals

- Provide generic Bézier curves and B-Splines. Due to their frequent usage, further provide lines, quadratic and cubic Bézier curves
- Support no-std for all targets
- Extensive unit testing and code coverage
- Integration tests for other generic math libraries (TBD - maybe optimath, aljabar, micromath, nalgebra) since Point types are replicated in many libraries

## Non-Goals

- Focus on use for rendering or highest performance (no GPU)

## Features

These are the main supported features. Some further utility methods are exposed where useful.  

### Quadratic and Cubic Bézier Curves

- [x] evaluation (De Casteljau, direct)
- [x] split
- [x] derivative
- [x] arc length (linear approx.)
- [ ] arc length (Legendre-Gauss)
- [ ] curvature/radius (Frenet-Serret Frame)
- [x] bounding box
- [ ] tight box

### Generic Bézier Curves

- [x] evaluation (De Casteljau)
- [x] split
- [x] derivative
- [x] arc length (linear approx.)
- [ ] arc length (Legendre-Gauss)
- [ ] curvature/radius (Frenet-Serret Frame)
- [x] bounding box
- [ ] tight box

### Generic B-Splines

- [x] evaluation (De Boor)
- [x] split
- [x] derivative
- [x] arc length (linear approx.; returns Result on eval errors)
- [ ] arc length (Legendre-Gauss)
- [ ] curvature/radius (Frenet-Serret Frame)
- [x] bounding box
- [ ] tight box

### Approximation Notes

- `arclen` and `distance_to_point*` use sampling-based approximations; increase steps for accuracy.
- `BSpline::arclen` returns `Result` because evaluation is fallible outside the knot domain.

## Examples

- `bezier_path`: mixed Bezier path segments; run `cargo run --example bezier_path`.
- `bspline_path`: basic B-spline path evaluation; run `cargo run --example bspline_path`.
- `plotters_cubic_bezier`: render a cubic Bezier with plotters; run `cargo run --example plotters_cubic_bezier`.
- `nalgebra_basic`: Bezier with nalgebra SVector; run `cargo run --example nalgebra_basic --features nalgebra`.
  
## Related  

If you're looking for a published crate for rendering with gpu support you should check out [Lyon](https://github.com/nical/lyon) from which I draw some inspiration, it's really good. It features lines, quadratic and cubic Béziers, in addition to arcs, triangles, and other geometric objects but no general Bézier curves or B-Splines. It also does seem to support wasm.  

Also, there's [Twinklebear/bspline](https://github.com/Twinklebear/bspline) which is a very clean and useful library for just bsplines. However, it depends on std and its simple Interpolate trait defines no way to access the individual dimensions of the points and hence implements no derivatives in the library.  

This clear online book [A Primer on Bézier Curves](https://pomax.github.io/bezierinfo/) helped me with a lot of the algorithms involved.
