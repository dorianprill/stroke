# Changelog

## [0.4.1] - Unreleased

### Added

- Added `kuva_basis` and a second README figure comparing cubic Bernstein and clamped B-spline basis weights, with interior knot guides.
- Added allocation-free Bernstein `basis_functions(t)` and `basis_functions_with_derivatives::<R>(t)` queries for generic, quadratic and cubic Béziers and line segments. Finite parameters may extrapolate; derivative orders above the degree return zero. `BezierBasisError` reports non-finite parameters or unrepresentable results. Existing curve evaluation and B-spline APIs are unchanged.
- Added `bounding_box_oriented_approx()` and `_with_budget(n)` to generic, quadratic and cubic Bézier curves and B-splines. Both return `OrientedBoundingBox` with containment and coordinate methods, or `OrientedBoundingBoxError`. Specialized Béziers delegate to the generic implementation. The shared search uses subdivision and plane rotations without allocation or new dependencies. B-splines include every nonzero span and both limits at jumps.
- Added `OrientedBoundingBoxDiagnostics` and `diagnostics()` to report model work, error estimates and subdivision limits. Larger budgets retain the best earlier box. Documented the method, references and stack usage. The search has no global-optimality or approximation-ratio guarantee, and its floating-point bounds are not certified.
- Added `oriented_box_study` for Bézier and B-spline accuracy and timings in 2D, 3D and 4D. It writes CSV diagnostics and checks faces against polynomial extrema. Tests also cover known minimum boxes, jumps, extreme scales and higher dimensions. Reference orientations are numerical results, not certified optima.

### Changed

- Excluded the generated `curves.svg` from the crates.io package; the README loads the image from GitHub.
- Show tight curve boxes and dashed control boxes together in the plotting example and README figure.
- Rewrote the README, API comments and example explanations. Clarified parameter units, trait requirements and numerical limits.

## [0.4.0] - 2026-10-17

### Breaking changes

- Renamed `PointN` and `point_generic` to `ArrayPoint` and `array_point`. The `Point` trait is unchanged.
- Removed specialized Bézier `dd()` methods; use `derivative_at(t)[axis]`. Generic Bézier `derivative()` now requires at least two control points; constant curves retain zero-valued `derivative_at()`.
- B-spline basis queries return `InvalidSpan` when the supplied span does not contain the parameter. Non-finite results return `NonFiniteComputation`.
- Coincidence predicates measure distance from the first control point; non-finite coincidence/linearity tolerances return false.
- Added root-search input and arithmetic errors. Root accuracy now bounds coordinate error at the returned parameter. A tiny parameter step alone no longer establishes success.
- Changed Bézier `control_points()` to borrowed iteration. `into_control_points()` returns the owned array. Replaced mutable generic Bézier iteration with validated control-point replacement and mapping.
- Added `Point::is_finite`; custom point implementations must report whether all coordinates are finite.
- Made line and Bezier `new` constructors fallible so they reject non-finite control points; generic Bezier construction also rejects an empty control-point array.
- Replaced `arclen(nsteps)` with `arclen(accuracy) -> Result<_, ArcLengthError>` across curves, segments, and paths; the old sampling calculation is now `arclen_polyline(segments)`.
- Length-to-parameter and length-to-point methods now require an accuracy target and return `Result`. Paths also use `Option` for empty lookups. Removed the old `_approx` inverse-length methods.
- Bounding-box methods now return `BoundingBox<Scalar, DIM>` instead of raw axis arrays; use `into_axes()` for the previous representation.
- B-spline construction rejects empty active knot domains and knot multiplicities greater than degree + 1. Partly clamped splines use the standard `[knots[D], knots[C]]` active domain.

### Added

- Added array conversions and borrowed iteration for `ArrayPoint`, plus lossless conversions between specialized and generic Béziers. Generic curves gained conditional debug and comparison traits. Paths gained clone and debug implementations.
- Added validated control-point mapping and replacement that leaves curves unchanged on error. Added clamped-uniform B-spline construction and direct `derivative_at` and `speed` methods. Path derivatives include the chain-rule scale.
- Added validated Bézier `subcurve` extraction and distance-based splitting. Paths gained capacity and reuse methods. `try_push` returns the supplied segment if the path is full.
- Bounding-box coordinate, center, extent, overlap, and intersection helpers, including touching and zero-width boxes.
- Added 3D motion, spatial-query and time-domain animation examples. Added external-consumer build tests with and without nalgebra.
- Adaptive Gauss-Legendre arc length, consistent partial-length APIs, and bracketed inverse-length lookup, shared across all curve types and paths.
- Allocation-free B-spline `bezier_spans()` extraction using the polar de Boor algorithm.
- Added tight boxes for generic Béziers, B-splines and paths. Added `control_box()` methods and a `BoundingBox` type with union and containment queries.
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
- Reuse endpoint-preserving interpolation for specialized evaluation and splitting to avoid NaNs from finite extreme controls.
- Compare geometric distances without squaring in coincidence/linearity predicates, preserving exact zero-tolerance distinctions and avoiding overflow/underflow.
- Validate B-spline span indices without overflow. Keep separate binary exponents during basis and derivative calculations to handle wide and tiny knot intervals.
- Validate Newton inputs and check coordinate error at the returned, domain-clamped parameter before reporting success. Share the solver across curve types and evaluate derivatives directly for constant pieces.
- Keep the best distance found during sampling and refinement, including exact endpoint and interior-sample hits. Use scaled Euclidean lengths and overflow-safe B-spline sampling parameters. The search remains heuristic.
- Scale and center line projections to avoid squared-length overflow, underflow and loss of small offsets on huge segments. Keep midpoint rounding errors and detect zero-length segments by exact equality.
- Verify release packages, warning-free all-feature docs and WASM/bare-metal ARM builds in CI, with and without nalgebra.
- Export public root-finding errors and common segment/B-spline error types at the crate root; all public errors implement `core::error::Error`.
- Correct the signal interpolation example using actual clamped basis functions and report singular/non-finite solver failures instead of skipping pivots.
- Align the README nalgebra dependency with the adapter's 0.34 dependency.
- Evaluate unclamped B-splines at the final domain parameter using de Boor's algorithm instead of incorrectly returning the last control point.
- Reject non-finite B-spline knots and control points, reporting the first invalid control-point index.
- Use exact comparisons for B-spline knot multiplicity and zero-width spans. Reserve scale-aware tolerances for algebraic cancellation.
- Preserve discontinuities when measuring paths and B-splines: gaps no longer acquire implicit connecting-line lengths.
- Replace derivative sign sampling and coordinate-scale-dependent root tolerances in bounds with Bernstein subdivision.
- Normalize analytic quadratic-root coefficients to avoid discriminant overflow/underflow, and use the cancellation-resistant quadratic formula.
- Default B-splines now have a valid clamped knot domain instead of an all-zero knot vector.
- Use the same representable path boundaries for evaluation and inverse-length lookup. Direct point lookup keeps the selected segment at jumps.
- Refine generic bounds in local coordinates, avoiding origin-dependent tolerances while preserving the original curve endpoints.
- Compute partial line lengths from scaled coordinate displacements, including representable portions of lines whose full length overflows. Linear Bézier and B-spline spans share this exact line formula.
- Preserve exact interpolation endpoints and avoid subtracting opposite finite extremes when their difference overflows.

## [0.3.0] - 2026-01-22

### Breaking changes

- Refactored Point trait bounds and curve trait requirements; downstream implementations may need updates.
- Removed the old Spline trait.
- `BSpline` arc-length and distance methods now return `Result` because evaluation can fail outside the knot domain.

### Added

- Generic Bezier and BSpline root finding via the FindRoot helper.
- Added `BSpline` basis functions, their derivatives and bounding boxes.
- BezierPath and BSplinePath types for multi-segment paths.
- nalgebra feature adapter for SVector integration.
- Added examples for scalar B-spline interpolation, a generic Bezier quarter circle, arc-length sampling and differential geometry.

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
- Added the first `root_newton_raphson` solver.
- Control point iterators for Bezier types.

### Changed

- Updated to Rust 2021 edition and refreshed documentation/readme.
- Enabled tinyvec const-generic support to improve PointN usability.

### Fixed

- Corrected generic Bezier derivative implementation.
- Resolved BSpline distance_to_point test failures and related edge cases.
- Addressed clippy warnings and architecture-specific constants.

See the [0.2.0 release](https://github.com/dorianprill/stroke/releases/tag/v0.2.0).

## [0.1.0] - 2022-02-24

### Added

- Initial no_std release with const-generic PointN and Point trait.
- Specialized line, quadratic, and cubic Bezier types with eval, split, derivative, arc length, and bounding box helpers.
- Generic Bezier and BSpline types with eval, split, derivative, and arc length support.
- Plotters-based example image and basic documentation.

See the [0.1.0 release](https://github.com/dorianprill/stroke/releases/tag/v0.1.0).
