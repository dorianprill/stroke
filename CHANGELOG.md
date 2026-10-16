# Changelog

All notable changes to this project will be documented in this file.

## [Unreleased]

### Breaking changes
- Added `Point::is_finite`; custom point implementations must report whether all coordinates are finite.
- Made line and Bezier `new` constructors fallible so they reject non-finite control points; generic Bezier construction also rejects an empty control-point array.

### Changed
- Added a project Clippy policy, Cargo lint aliases, and a pre-commit hook; CI checks both feature configurations with warnings treated as errors.
- Replaced the Plotters example with a Kuva SVG figure showing two generic Bezier curves and clamped/unclamped generic B-splines, their bounding boxes, labeled control points, knot-vector legends, and marked interior knots and unclamped endpoints.
- Removed Plotters and the CI fontconfig installation; the plotting example uses Kuva without optional features.
- Configured CI to use the compatible nightly trait solver for `generic_const_exprs`.
- Cleaned up existing curve-code Clippy findings and added strict Clippy checks to CI.

### Fixed
- Evaluate unclamped B-splines at the final domain parameter using de Boor's algorithm instead of incorrectly returning the last control point.
- Reject non-finite B-spline knots and control points, reporting the first invalid control-point index.
- Treat B-spline knot multiplicity and zero-width spans structurally, while using scale-aware tolerances only for algebraic cancellation.

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
