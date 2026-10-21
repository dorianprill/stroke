# Changelog

All notable changes to this project will be documented in this file.

## [Unreleased]

### Added

- Added `kuva_basis` and a second README figure comparing cubic Bernstein and clamped B-spline basis weights, with interior knot guides.
- Added allocation-free Bernstein `basis_functions(t)` and `basis_functions_with_derivatives::<R>(t)` queries for generic, quadratic and cubic Béziers and line segments. Finite parameters may extrapolate; derivative orders above the degree return zero. `BezierBasisError` reports non-finite parameters or unrepresentable results. Existing curve evaluation and B-spline APIs are unchanged.

## [0.4.0] - 2026-10-17

### Breaking changes
- Renamed `PointN` and `point_generic` to `ArrayPoint` and `array_point` for easier distinction between implementation and trait; the `Point` trait is unchanged.
- Removed specialized Bézier `dd()` methods; use `derivative_at(t)[axis]`. Generic Bézier `derivative()` now requires at least two control points; constant curves retain zero-valued `derivative_at()`.
- B-spline basis queries reject spans that do not own the parameter with `InvalidSpan`; non-finite results return `NonFiniteComputation`.
- Coincidence predicates measure distance from the first control point; non-finite coincidence/linearity tolerances return false.
- Added root-search input/arithmetic error variants and defined root accuracy as an absolute coordinate residual tolerance; tiny parameter steps no longer establish success.
- Changed Bézier `control_points()` to borrowed iteration; `into_control_points()` returns the owned array. Removed mutable generic Bézier iteration in favor of checked control-point replacement/mapping.
- Added `Point::is_finite`; custom point implementations must report whether all coordinates are finite.
- Made line and Bezier `new` constructors fallible so they reject non-finite control points; generic Bezier construction also rejects an empty control-point array.
- Replaced `arclen(nsteps)` with `arclen(accuracy) -> Result<_, ArcLengthError>` across curves, segments, and paths; the old sampling calculation is now `arclen_polyline(segments)`.
- Length-to-parameter and length-to-point methods now require an accuracy target and return `Result`; paths additionally retain `Option` for empty lookups. Removed the old `_approx` inverse-length methods.
- Bounding-box methods now return `BoundingBox<Scalar, DIM>` instead of raw axis arrays; use `into_axes()` for the previous representation.
- B-spline construction rejects empty active knot domains and knot multiplicities greater than degree + 1. Partly clamped splines use the standard `[knots[D], knots[C]]` active domain.

### Added
- Array conversions and borrowed iteration for `ArrayPoint`, lossless specialized/generic Bézier conversions, and conditional debugging/comparison traits for generic curves and cloning/debugging for paths.
- Validated control-point mapping and atomic replacement, clamped-uniform B-spline construction, and direct `derivative_at`/`speed` helpers, including path chain-rule scaling.
- Checked Bézier `subcurve` extraction and distance-based splitting; fixed-capacity path reuse helpers and a segment-preserving `try_push` API.
- Bounding-box coordinate, center, extent, overlap, and intersection helpers, including touching and zero-width boxes.
- 3D motion, spatial-query, and time-domain animation examples, and external-consumer smoke coverage with/without nalgebra.
- Adaptive Gauss-Legendre arc length, consistent partial-length APIs, and bracketed inverse-length lookup, shared across all curve types and paths.
- Allocation-free B-spline `bezier_spans()` extraction using the polar de Boor algorithm.
- Tight generic Bézier, B-spline, and path boxes; explicit `control_box()` methods and a `BoundingBox` type with union and containment helpers.
- Stable `PointNorm::norm()` accumulation for very small and very large representable magnitudes.

### Changed
- Specialized Bézier `eval()` and `eval_casteljau()` share endpoint-preserving de Casteljau interpolation. Array-backed point arithmetic uses array mapping/iteration.
- Condensed README usage guidance and documented the time/parameter units in motion and animation examples.
- Added a project Clippy policy, Cargo lint aliases, and a pre-commit hook; CI checks both feature configurations with warnings treated as errors.
- Replaced the Plotters example with a Kuva SVG figure showing two generic Bezier curves and clamped/unclamped generic B-splines, their bounding boxes, labeled control points, knot-vector legends, and marked interior knots and unclamped endpoints.
- Removed Plotters and the CI fontconfig installation; the plotting example uses Kuva without optional features.
- Configured CI to use the compatible nightly trait solver for `generic_const_exprs`.
- Cleaned up existing curve-code Clippy findings and added strict Clippy checks to CI.

### Fixed
- Prevent empty generic derivative curves. Scaled differences and shared geometry helpers stabilize tangents, normals, and curvature across coordinate scales, including overflowing derivative magnitudes.
- Reuse robust interpolation for specialized evaluation and splitting; avoid NaNs from finite extreme controls.
- Compare geometric distances without squaring in coincidence/linearity predicates, preserving exact zero-tolerance distinctions and avoiding overflow/underflow.
- Validate B-spline span ownership without index overflow and carry separate exponents through basis/derivative recurrences to preserve wide and tiny knot intervals.
- Validate Newton inputs and report success only for a finite satisfactory residual at the returned domain-clamped parameter. Share the solver across curve types and support direct derivative evaluation for constant pieces.
- Preserve every better coarse/refinement distance candidate, including exact endpoint and interior-sample hits, while retaining the heuristic search contract. Use stable Euclidean norms and overflow-safe B-spline sampling parameters.
- Scale and center line projections to avoid squared-length overflow/underflow and loss of small central offsets on huge segments; retain midpoint rounding tails and distinguish zero-length segments structurally.
- Verify release packages, warning-free all-feature docs, and library portability with/without nalgebra on WASM and bare-metal ARM in CI.
- Export public root-finding errors and common segment/B-spline error types at the crate root; all public errors implement `core::error::Error`.
- Correct the signal interpolation example using actual clamped basis functions and report singular/non-finite solver failures instead of skipping pivots.
- Align the README nalgebra dependency with the adapter's 0.34 dependency.
- Evaluate unclamped B-splines at the final domain parameter using de Boor's algorithm instead of incorrectly returning the last control point.
- Reject non-finite B-spline knots and control points, reporting the first invalid control-point index.
- Treat B-spline knot multiplicity and zero-width spans structurally, while using scale-aware tolerances only for algebraic cancellation.
- Preserve discontinuities when measuring paths and B-splines: gaps no longer acquire implicit connecting-line lengths.
- Replace derivative sign sampling and coordinate-scale-dependent root tolerances in bounds with Bernstein subdivision.
- Normalize analytic quadratic-root coefficients to avoid discriminant overflow/underflow, and use the cancellation-resistant quadratic formula.
- Default B-splines now have a valid clamped knot domain instead of an all-zero knot vector.
- Preserve the selected path segment at discontinuities using shared representable parameter boundaries and direct segment-local point lookup.
- Refine generic bounds in local coordinates, avoiding origin-dependent tolerances while preserving the original curve endpoints.
- Compute partial line lengths from scaled coordinate displacements, including representable portions of lines whose full length overflows. Linear Bézier and B-spline spans share this exact line formula.
- Preserve exact interpolation endpoints and avoid subtracting opposite finite extremes when their difference overflows.

## [0.3.0] - 2026-01-22

### Breaking changes
- Refactored Point trait bounds and curve trait requirements; downstream implementations may need updates.
- Removed the old Spline trait.
- BSpline arc length and distance helpers are now fallible and return Result when evaluation can fail outside the knot domain.

### Added
- Generic Bezier and BSpline root finding via the FindRoot helper.
- BSpline basis functions (with derivatives) and bounding box support.
- BezierPath and BSplinePath types for multi-segment paths.
- nalgebra feature adapter for SVector integration.
- New examples: B-spline 1D interpolation, generic Bezier quarter-circle, arc-length sampling, and tangent/normal/curvature sampling.

### Changed
- Upgraded to Rust 2024 edition and pinned a nightly toolchain for const-generic features.
- Improved arc-length and distance approximation routines across curve types.
- Plotters example renamed/excluded from tests; example documentation expanded.

### Fixed
- Corrected De Casteljau and line math edge cases.
- Fixed B-spline derivative edge cases, knot span search, and out-of-bounds handling.
- Addressed clippy lint warnings and rustdoc issues; expanded CI coverage.

## [0.2.0] - 2024-09-27

### Breaking changes
- Major API rework of Point trait bounds and curve generics.
- BSpline const generic switched from "order" to "degree" and related type parameters updated.
- Removed unfinished RationalBezier type.
- Added BSplineError / KnotVectorKind, changing BSpline construction and validation behavior.

### Added
- distance_to_point helpers for Bezier, BSpline, QuadraticBezier, and CubicBezier.
- Bezier arc length approximation.
- BSpline derivative support.
- Root solver module groundwork (root_newton_raphson).
- Control point iterators for Bezier types.

### Changed
- Updated to Rust 2021 edition and refreshed documentation/readme.
- Enabled tinyvec const-generic support to improve PointN usability.

### Fixed
- Corrected generic Bezier derivative implementation.
- Resolved BSpline distance_to_point test failures and related edge cases.
- Addressed clippy warnings and architecture-specific constants.

See: https://github.com/dorianprill/stroke/releases/tag/v0.2.0

## [0.1.0] - 2022-02-24

### Added
- Initial no_std release with const-generic PointN and Point trait.
- Specialized line, quadratic, and cubic Bezier types with eval, split, derivative, arc length, and bounding box helpers.
- Generic Bezier and BSpline types with eval, split, derivative, and arc length support.
- Plotters-based example image and basic documentation.

See: https://github.com/dorianprill/stroke/releases/tag/v0.1.0
