# brezel  
![Rust](https://github.com/dorianprill/brezel/workflows/Rust/badge.svg)  
WIP minimalistic implementation of up to cubic bezier curves (on generic 2D-point types with a minimal exposed trait) for general use on all platforms including embedded and wasm through #![no_std]

Brezel is the German word for "pretzel" - I'm sure you appreciate the parallels to those delicious curvy snacks

![A Cubic Bezier Curve with Bounding Box and Convex Hull rendered by plotters.rs](https://github.com/dorianprill/brezel/blob/master/cubic_bezier_bounding_box.png)  

made with [plotters.rs](https://github.com/38/plotters)  

## Goals

- [x] Support no-std for both 32 and 64 bit targets
- [x] Provide lines
- [x] Provide quadratic beziers
- [x] Provide cubic beziers
- [x] Where applicable: evaluate, split, arc length, curvature/radius, bounding box, derivative functions and shortcuts for evaluating them
- [ ] tight box, the curve's normal
- [ ] Integrate well with other libraries like micromath, nalgebra etc. since Point2/3 types are reproduced in many libraries. Achieve this by using generics as much as possible and expose only a minimal trait to the outside
- [ ] Good test coverage for both unit and integration test
Also for me to truly learn rusts' advanced concepts as well as touch up on bernstein basis/bezier curves, maybe b-splines/NURBS.

## Non-Goals

- Threading or parallelization  
- Focus on use for rendering
- Pushing it to crates.io and maintaining it

## Related 
If you're looking for a published crate with gpu support you should check out [Lyon](https://github.com/nical/lyon) from which I draw some inspiration, it's really good.
  
