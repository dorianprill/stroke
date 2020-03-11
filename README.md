# brezel

WIP minimalistic implementation of (up to) cubic bezier curves (on generic point types with a minimal exposed trait) for general use

## Goals

- [ ] Support no-std for both 32 and 64 bit targets
- [x] Provide (up to) cubic beziers in 2d space
- [ ] Provide (up to) cubic beziers in 3d space (Maybe even make the dimension generic?)
- [ ] Integrate well with other libraries like micromath, nalgebra etc. since Point2/3 types are reproduced in many libraries. Achieve this by using generics as much a spossible and expose a minimal trait to the outside

Also for me to truly learn rusts' advanced concepts as well as touch up on bernstein basis/bezier curves, maybe b-splines/NURBS.

## Non-Goals

- No fancy muli threading or parallelization  
- No focus on use for rendering
- pushing it to crate.io and supporting it
  