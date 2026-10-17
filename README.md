# stroke  

![Rust](https://github.com/dorianprill/stroke/actions/workflows/rust.yml/badge.svg)
![crates.io](https://img.shields.io/crates/v/stroke.svg)

A zero-allocation library providing const-generic implementations of Bézier curves, B-Spline curves and specialized implementations of up to cubic Bézier curves in N-dimensional euclidean space. It is intended for general/embedded/wasm use supporting #![no_std] environments written in 100% safe Rust with minimal dependencies.

![Generic Bézier and B-spline curves with bounding boxes and labeled control points](curves.svg)

Made with [Kuva](https://github.com/Psy-Fer/kuva). Regenerate with `cargo run --example kuva_curves`.

The library makes heavy use of const-generics and `generic_const_exprs`, so the nightly compiler is required.  
This repo uses `nightly` via `rust-toolchain`.  
Consumers using const-expression-based APIs (such as bounding boxes) must also
enable `#![feature(generic_const_exprs)]` at their crate root on nightly.
`ArrayPoint` supplies an array-backed implementation of `Point` for use without an external math library.
`ArrayPoint<T, N>` uses your chosen scalar `T` (typically `f32`/`f64`) and keeps the dimension at compile time.
Point scalar types must implement `num_traits::Float`; the library therefore currently targets floating-point numeric representations such as `f32` and `f64`.

Should you want to integrate with types provided by another library, implement the small `Point` trait. Optional extension traits (`PointIndex`, `PointDot`, `PointNorm`) unlock component access and geometric helpers when needed.  

Right now, the generic versions don't implement all methods that the specialized versions do (as the algorithms get a bit more complicated) but should reach parity eventually.

## Integration features

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
]).unwrap();

let mid = curve.eval(0.5);
```

Enable with:

```toml
stroke = { version = "0.4", features = ["nalgebra"] }
nalgebra = { version = "0.34", default-features = false, features = ["libm"] }
```

## Goals

- Provide generic Bézier curves and B-Splines. Due to their frequent usage, further provide lines, quadratic and cubic Bézier curves
- Support no-std for all targets
- Extensive unit testing and code coverage
- Integration tests for other generic math libraries (currently only nalgebra)

## Non-Goals

- Focus on use for rendering or highest performance (no GPU)

## Features

These are the main supported features. Some further utility methods are exposed where useful.  

### Quadratic and Cubic Bézier Curves

- [x] evaluation (De Casteljau)
- [x] split
- [x] derivative
- [x] tangent (unit)
- [x] normal (principal)
- [x] curvature (magnitude)
- [x] arc length (linear approx.)
- [x] arc length (adaptive Gauss-Legendre)
- [ ] Frenet-Serret frame / torsion
- [x] bounding box
- [x] tight axis-aligned bounding box and conservative control-point bounds

### Generic Bézier Curves

- [x] evaluation (De Casteljau)
- [x] split
- [x] derivative
- [x] tangent (unit)
- [x] normal (principal)
- [x] curvature (magnitude)
- [x] arc length (linear approx.)
- [x] arc length (adaptive Gauss-Legendre)
- [ ] Frenet-Serret frame / torsion
- [x] bounding box
- [x] tight axis-aligned bounding box and conservative control-point bounds

### Generic B-Splines

- [x] evaluation (De Boor)
- [x] Bézier span extraction (`bezier_spans`, without allocation)
- [x] derivative
- [x] basis functions (and derivatives)
- [x] tangent (unit)
- [x] normal (principal)
- [x] curvature (magnitude)
- [x] arc length (polyline approximation)
- [x] arc length (adaptive Gauss-Legendre)
- [ ] Frenet-Serret frame / torsion
- [x] bounding box
- [x] tight axis-aligned bounding box and conservative control-point bounds

### General use

`arclen(accuracy)`, `arclen_between(start, end, accuracy)`, and inverse-length
queries use adaptive Gauss-Legendre quadrature. Choose a positive finite
`accuracy` in coordinate units; subdivision follows an error estimate and is
limited by scalar roundoff. These queries return errors for invalid input,
non-finite arithmetic, or exhausted convergence budgets.
`arclen_polyline(segments)` uses fixed subdivisions per curve or B-spline span.

```rust
use stroke::{ArrayPoint, Bezier};

let curve = Bezier::new([
    ArrayPoint::from([0.0_f64, 0.0]),
    ArrayPoint::from([1.0, 2.0]),
    ArrayPoint::from([2.0, 0.0]),
]).unwrap();
let length = curve.arclen(1e-10).unwrap();
let halfway = curve.point_at_length(length * 0.5, 1e-10).unwrap();
let bounds = curve.bounding_box();
```

Bézier parameters lie in `[0, 1]`; paths divide that interval uniformly among
segments. B-splines use `[knots[D], knots[C]]`. `derivative_at(t)` and `speed(t)`
use these parameter units. At a boundary, evaluation selects the piece on the
right, except at the final endpoint. Gaps contribute no arc length.

Inverse lookups clamp finite distances to the endpoints; zero-length curves
select their start. Empty paths have length zero, `Ok(None)` for inverse
lookups, and `None` for boxes. `bounding_box()` gives tight bounds to floating-point
resolution; the cheaper `control_box()` encloses the control points and curve.

Clamped uniform B-splines interpolate their endpoint controls; interior
controls shape the curve. Use `bezier_spans()` to extract their polynomial
pieces. Mapping control points transforms the whole curve geometrically when
the map is affine.

Curved `distance_to_point*` queries use sampling and local refinement, with no
global accuracy guarantee. Non-finite query points produce NaN.
Newton root queries use `eps` as an absolute coordinate residual tolerance;
a failed local search can miss existing roots. At stationary points,
tangent and curvature return zero and normal returns `None`.

## Examples

- `bezier_path`: mixed Bezier path segments; run `cargo run --example bezier_path`.
- `bspline_path`: basic B-spline path evaluation; run `cargo run --example bspline_path`.
- `bspline_signal_1d`: 1D signal interpolation with a cubic B-spline; run `cargo run --example bspline_signal_1d`.
- `bezier_quarter_circle`: cubic Bezier quarter-circle approximation; run `cargo run --example bezier_quarter_circle`.
- `arc_length_sampling`: equal-arc-length sampling along a Bezier curve; run `cargo run --example arc_length_sampling`.
- `motion_3d`: sample a 3D tool trajectory at a fixed feed rate; CSV includes time, position, velocity, and curvature. Run `cargo run --example motion_3d`.
- `spatial_queries`: CAD/GIS broad-phase filtering with tight versus control boxes, overlap queries, and trimmed curve bounds. Run `cargo run --example spatial_queries`.
- `animation_envelope`: scalar B-spline value and rate over a native time domain, producing CSV output. Run `cargo run --example animation_envelope`.
- `tangent_normal_curvature`: sample tangent/normal/curvature; run `cargo run --example tangent_normal_curvature`.
- `kuva_curves`: render a 2×2 SVG figure with generic quadratic and cubic `Bezier` curves, a clamped quintic `BSpline`, and an unclamped cubic `BSpline`, with curve bounding boxes, control polygons, named control points, knot-vector legends, and marked interior knots and unclamped endpoints; run `cargo run --example kuva_curves` to write `curves.svg`. Uses Kuva's SVG backend with optional features disabled; no system font library or font installation is required.
- `nalgebra_basic`: Bezier with nalgebra SVector; run `cargo run --example nalgebra_basic --features nalgebra`.
  
## Development

Install the nightly tools once with `rustup component add rustfmt clippy`.
Run `cargo lint` for all targets with all features, or `cargo lint-default`
for all targets without default features. Both treat warnings as errors;
CI runs both commands. `.cargo/config.toml` also selects the compatible
nightly trait solver for local builds, Clippy, and doctests.
CI checks the library with and without nalgebra on WASM and bare-metal ARM,
verifies the release package, and denies warnings in all-feature rustdoc.
Run `cargo test --examples --all-features` for example assertions. CI also
builds a separate downstream consumer from outside the repository, without
inheriting this configuration. The API pass was tested with
`rustc 1.101.0-nightly (7df81381d 2026-10-14)`; stable Rust is not supported.

Enable automatic checks before each commit in a fresh clone:

```sh
git config --local core.hooksPath .githooks
```

The hook checks formatting and runs both lint commands. It requires no extra
hook manager. To disable it, run `git config --local --unset core.hooksPath`.

Lint levels are in `Cargo.toml`, with lint options in `clippy.toml`.
The policy keeps Clippy's default checks, adds Cargo metadata checks,
portable `core`/`alloc` API checks, documentation formatting, redundant-clone
checks, and checks for debug macros, unfinished code, and process exits.
Unsafe code is forbidden for every target. Clippy reads the declared MSRV
from `package.rust-version` and checks test code against it as well.

We enable extra lints individually rather than enabling all of `pedantic`,
`nursery`, or `restriction`. Exact floating-point knot comparisons, indexed
arrays, scalar casts, and short mathematical variable names are intentional
in these algorithms. Unwraps are also intentional for scalar conversions
and test assertions; they are not banned globally. Any lint exception should
be scoped to the relevant item and explain why it applies.

## Related  

If you're looking for a published crate for rendering with gpu support you should check out [Lyon](https://github.com/nical/lyon) from which I draw some inspiration, it's really good. It features lines, quadratic and cubic Béziers, in addition to arcs, triangles, and other geometric objects but no general Bézier curves or B-Splines. It also does seem to support wasm.  

Also, there's [Twinklebear/bspline](https://github.com/Twinklebear/bspline) which is a very clean and useful library for just bsplines. However, it depends on std and its simple Interpolate trait defines no way to access the individual dimensions of the points and hence implements no derivatives in the library.  

[Kurbo](https://github.com/linebender/kurbo) provides 2D curves and vector paths using `f64`, with a focus on accuracy and efficient high-precision queries.

This clear online book [A Primer on Bézier Curves](https://pomax.github.io/bezierinfo/) helped me with a lot of the algorithms involved.
