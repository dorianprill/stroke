# stroke-rs  
![Rust](https://github.com/dorianprill/brezel/workflows/Rust/badge.svg)  
WIP minimalistic implementation of up to cubic bezier curves (on generic point types with a minimal exposed trait) for general use on all platforms including embedded and wasm through #![no_std]

![A Cubic Bezier Curve with Bounding Box and Convex Hull rendered by plotters.rs](https://github.com/dorianprill/brezel/blob/master/cubic_bezier_bounding_box.png)  

made with [plotters.rs](https://github.com/38/plotters)  

## Goals

- [x] Support no-std for both 32 and 64 bit targets
- [x] Provide lines
- [x] Provide quadratic beziers
- [x] Provide cubic beziers
- [x] Where applicable: evaluate, split, arc length, curvature/radius, bounding box, derivative functions and shortcuts for evaluating them
- [ ] Where applicable: tight box, the curve's normal
- [x] Provide general B-Splines (wip using const generics)
- [ ] Integrate well with other libraries like (possibly) micromath, nalgebra etc. since Point2/3 types are reproduced in many libraries. Achieve this by using generics as much as possible and expose only a minimal Point-trait to the outside. Only possible with libraries that also do not make a distinction between a point and its position vector
- [ ] Good test coverage for both unit and integration test
My personal goal is to truly learn rusts' advanced concepts as well as touch up on bernstein basis/bezier curves, B-Splines, maybe even NURBS.

## Non-Goals

- Focus on use for rendering or highest performance

## Related 
If you're looking for a published crate with gpu support you should check out [Lyon](https://github.com/nical/lyon) from which I draw some inspiration, it's really good.  
Also, there's [Twinklebear/bspline](https://github.com/Twinklebear/bspline) which is a very clean and useful library for just bsplines.  
This clear online book [A Primer on Bézier Curves](https://pomax.github.io/bezierinfo/) helped me with a lot of the algorithms involved.
