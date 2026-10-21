# stroke  

![Rust](https://github.com/dorianprill/stroke/actions/workflows/rust.yml/badge.svg)
![crates.io](https://img.shields.io/crates/v/stroke.svg)

A zero-allocation library providing const-generic implementations of Bézier curves, B-Spline curves and specialized implementations of up to cubic Bézier curves in N-dimensional euclidean space. It is intended for general/embedded/wasm use supporting #![no_std] environments written in 100% safe Rust with minimal dependencies.

![Generic Bézier and B-spline curves with tight and control boxes and labeled control points](https://raw.githubusercontent.com/dorianprill/stroke/main/curves.svg)

Solid green boxes show `bounding_box()`; dashed gray boxes show `control_box()`. Made with [Kuva](https://github.com/Psy-Fer/kuva). Regenerate with `cargo run --example kuva_curves`.

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

☑ implemented; ☐ not implemented; — not applicable. Checkmarks refer to public APIs.

| Feature | Quadratic Bézier | Cubic Bézier | Generic Bézier | Generic B-Splines |
|---|:---:|:---:|:---:|:---:|
| evaluation | ☑ | ☑ | ☑ | ☑ |
| split | ☑ | ☑ | ☑ | ☐ |
| derivative | ☑ | ☑ | ☑ | ☑ |
| tangent (unit) | ☑ | ☑ | ☑ | ☑ |
| normal (principal) | ☑ | ☑ | ☑ | ☑ |
| curvature (magnitude) | ☑ | ☑ | ☑ | ☑ |
| arc length (polyline approximation) | ☑ | ☑ | ☑ | ☑ |
| arc length (adaptive Gauss-Legendre) | ☑ | ☑ | ☑ | ☑ |
| Frenet-Serret frame / torsion | ☐ | ☐ | ☐ | ☐ |
| bounding box | ☑ | ☑ | ☑ | ☑ |
| tight axis-aligned bounding box and conservative control-point bounds | ☑ | ☑ | ☑ | ☑ |
| approximate oriented bounding box (heuristic orientation search) | ☑ | ☑ | ☑ | ☑ |
| basis functions (and derivatives) | ☑ | ☑ | ☑ | ☑ |

Béziers use de Casteljau evaluation; B-splines use de Boor. Polyline arc length
uses `segments` intervals per Bézier curve or per nonzero B-spline knot span.

## Basis functions

Each color shows one control point's weight in curve evaluation. Both sets of
weights sum to one on `[0, 1]`; dashed lines mark the B-spline's interior knots.

![Cubic Bézier and B-spline basis weights in two side-by-side panels](https://raw.githubusercontent.com/dorianprill/stroke/main/basis.svg)

For a cubic Bézier with controls $P_0,\ldots,P_3$, the Bernstein weights are

$$
\begin{aligned}
B_0^3(t)&=(1-t)^3, & B_1^3(t)&=3t(1-t)^2,\\
B_2^3(t)&=3t^2(1-t), & B_3^3(t)&=t^3.
\end{aligned}
$$

The curve point is $C(t)=\sum_{i=0}^{3}B_i^3(t)P_i$.

A cubic B-spline uses $C(t)=\sum_i N_{i,3}(t)P_i$, with at most four nonzero
weights at each parameter. Given knots $u_i$, start with $N_{i,0}(t)=1$ on
$[u_i,u_{i+1})$ and zero elsewhere, then apply the Cox-de Boor recurrence for
$p=1,2,3$:

$$
N_{i,p}(t)=\frac{t-u_i}{u_{i+p}-u_i}N_{i,p-1}(t)
+\frac{u_{i+p+1}-t}{u_{i+p+1}-u_{i+1}}N_{i+1,p-1}(t).
$$

A term with a zero denominator contributes zero. At the final domain endpoint,
use the left-hand limit. The plot uses clamped knots
`[0, 0, 0, 0, 0.25, 0.5, 0.75, 1, 1, 1, 1]`.

Béziers and line segments expose Bernstein weights through `basis_functions(t)`
and `basis_functions_with_derivatives::<R>(t)`. Values follow control-point order;
derivative row `r` is with respect to `t`, and orders above the degree are zero.
Finite parameters outside `[0, 1]` extrapolate. Non-finite inputs or results
return `BezierBasisError`. Use `eval(t)` when only a curve point is needed.

Regenerate the figure with `cargo run --example kuva_basis`.

## General use

Constructors return `Result` and reject non-finite controls. B-splines also
validate their knots. Control-point setters leave the curve unchanged on error.

```rust
#![allow(incomplete_features)]
#![feature(generic_const_exprs)]
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

`arclen(accuracy)`, `arclen_between(start, end, accuracy)` and inverse-length
queries use adaptive Gauss-Legendre integration for curved pieces. Choose a
finite, positive `accuracy` in coordinate units. It is an absolute error target
based on an estimate, limited by floating-point roundoff. Invalid inputs,
non-finite arithmetic or exhausted work limits return errors.
`arclen_polyline(segments)` uses fixed parameter steps per curve or B-spline span.

Bézier parameters lie in `[0, 1]`. Paths divide that interval equally among
stored segments. B-splines use `[knots[D], knots[C]]`. Derivatives and speeds
use these parameter units. At a boundary, evaluation selects the piece to the
right, except at the final endpoint. Gaps add no arc length.

Inverse-length queries clamp finite distances to the endpoints. Zero-length
curves select their start. Empty paths have length zero, return `Ok(None)`
for inverse lookups and `None` for boxes.

`bounding_box()` gives tight axis-aligned curve bounds to floating-point
resolution. `control_box()` bounds the controls and curve, usually with less work.
`bounding_box_oriented_approx()` searches for a smaller box with rotated axes,
using up to 256 full-curve trials. Its `_with_budget(n)` variant limits trials.
Neither the runtime nor volume error is fixed by that budget. Larger budgets
retain earlier improvements. The search may stop before reaching even a local
minimum and has no fixed-factor approximation guarantee. `diagnostics()` reports
model work and error estimates. B-splines include every nonzero span and both limits at jumps.
Their stack storage grows with control count and degree.

Clamped uniform B-splines pass through their endpoint controls. Interior controls
shape the curve. Use `bezier_spans()` to extract its polynomial pieces.
Affine control-point maps, such as translations, also map the whole curve.

Curved `distance_to_point*` queries sample and refine near the closest sample.
They can miss a closer point and have no error bound. Non-finite query points
produce NaN. Newton root searches use `eps` to bound coordinate error, and can
fail even when a root exists. At a stationary point, tangent and curvature
return zero and normal returns `None`.

## Examples

Run `cargo run --example NAME`. The study and nalgebra examples need the extra
options shown below.

- `bezier_path` joins line, quadratic and cubic segments.
- `bspline_path` evaluates a path of B-spline segments.
- `bspline_signal_1d` fits a cubic B-spline through scalar samples.
- `bezier_quarter_circle` approximates a quarter circle with a cubic curve.
- `arc_length_sampling` places samples at equal arc-length intervals.
- `motion_3d` samples a tool trajectory at a fixed feed rate. Its CSV contains
  time, position, velocity and curvature.
- `spatial_queries` filters curve candidates with bounding-box overlap and
  bounds trimmed portions of a curve.
- `animation_envelope` evaluates a scalar B-spline and its rate of change over
  a time domain, with CSV output.
- `tangent_normal_curvature` samples those geometric quantities.
- `kuva_basis` writes `basis.svg`, comparing cubic Bernstein and clamped
  B-spline basis weights, with interior knot guides.
- `kuva_curves` writes `curves.svg`, with curve boxes, control boxes, labeled
  controls and knot markers. It needs no system font libraries.
- `oriented_box_study` checks oriented bounds against polynomial extrema and
  numerical reference searches in 2D, 3D and 4D. Run it with
  `cargo run --release --example oriented_box_study`. Add `-- --extended`
  for a larger corpus. It writes timings and diagnostics to CSV. Reference
  orientations are best-known numerical results, not certified optima.
- `nalgebra_basic` uses `SVector` points. Run it with
  `cargo run --example nalgebra_basic --features nalgebra`.

## Development

Install the tools with `rustup component add rustfmt clippy`.

- `cargo lint` checks all targets with all features.
- `cargo lint-default` checks all targets without default features.
- `cargo test --all-features` runs the tests and doc examples.
- `cargo test --examples --all-features` runs the example assertions.

Both lint commands deny warnings. CI also checks WASM and bare-metal ARM builds
with and without nalgebra, verifies the release package and denies rustdoc
warnings. A separate consumer test runs outside the repository's Cargo
configuration. Local builds use the nightly trait solver selected in
`.cargo/config.toml`.

To check formatting and run both lints before each commit:

```sh
git config --local core.hooksPath .githooks
```

The hook needs no extra manager. Disable it with
`git config --local --unset core.hooksPath`.

Lint settings are in `Cargo.toml` and `clippy.toml`. Unsafe code is forbidden in
the library, tests and examples. Clippy checks the declared minimum Rust version
in `package.rust-version`, including test code.

Extra lints are enabled individually. Exact knot comparisons, indexed arrays,
scalar casts and short mathematical names are intentional. Scalar conversions
and test assertions may use `unwrap()`. Keep lint exceptions local and explain
why they are needed.

## Related  

If you're looking for a published crate for rendering with gpu support you should check out [Lyon](https://github.com/nical/lyon) from which I draw some inspiration, it's really good. It features lines, quadratic and cubic Béziers, in addition to arcs, triangles, and other geometric objects but no general Bézier curves or B-Splines. It also does seem to support wasm.  

Also, there's [Twinklebear/bspline](https://github.com/Twinklebear/bspline) which is a very clean and useful library for just bsplines. However, it depends on std and its simple Interpolate trait defines no way to access the individual dimensions of the points and hence implements no derivatives in the library.  

[Kurbo](https://github.com/linebender/kurbo) provides 2D curves and vector paths using `f64`, with a focus on accuracy and efficient high-precision queries.

This clear online book [A Primer on Bézier Curves](https://pomax.github.io/bezierinfo/) helped me with a lot of the algorithms involved.
