# brezel

WIP minimalistic implementation of (up to) cubic bezier curves (on generic point types with a minimal exposed trait) for general use

## Goals

- [x] Support no-std for both 32 and 64 bit targets
- [ ] Provide quadratic beziers
- [x] Provide cubic beziers
- [ ] Integrate well with other libraries like micromath, nalgebra etc. since Point2/3 types are reproduced in many libraries. Achieve this by using generics as much as possible and expose only a minimal trait to the outside
- [ ] Good test coverage for both unit and integration test
Also for me to truly learn rusts' advanced concepts as well as touch up on bernstein basis/bezier curves, maybe b-splines/NURBS.

## Non-Goals

- Threading or parallelization  
- Focus on use for rendering
- Pushing it to crates.io and maintaining it
  