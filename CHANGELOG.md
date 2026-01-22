# Changelog

All notable changes to this project will be documented in this file.

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
