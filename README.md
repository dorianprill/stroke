# stroke-rs  

![Rust](https://github.com/dorianprill/brezel/workflows/Rust/badge.svg)  
General implementations of Bézier curves and B-Splines and specialized implementations of up to cubic Bézier curves (on generic point types [WIP] with a minimal exposed trait) for general use on all platforms including embedded and wasm through #![no_std]

![A Cubic Bézier Curve with Bounding Box and Convex Hull rendered by plotters.rs](https://raw.githubusercontent.com/dorianprill/stroke-rs/main/cubic_bezier_bounding_box.png)  

made with [plotters.rs](https://github.com/38/plotters)  

## Goals

- [x] Support no-std for both 32 and 64 bit targets
- [x] Provide lines
- [x] Provide quadratic Béziers
- [x] Provide cubic Béziers
- [x] Where applicable: evaluate, split, arc length, curvature/radius, bounding box, derivative functions and shortcuts for evaluating them
- [ ] Where applicable: tight box, the curve's normal
- [x] Provide general Bézier curves (WIP using const generics)
- [x] Provide general B-Splines (WIP using const generics)
- [ ] Integrate well with other libraries like (possibly) micromath, nalgebra etc. since Point2/3 types are reproduced in many libraries. Achieve this by using generics as much as possible and expose only a minimal Point-trait to the outside. Only possible with libraries that also do not make a distinction between a point and its position vector
- [ ] Good test coverage for both unit and integration test
My personal goal is to truly learn rusts' advanced concepts as well as touch up on bernstein basis/Bézier curves, B-Splines, maybe even NURBS.

## Non-Goals

- Focus on use for rendering or highest performance

## Related  

If you're looking for a published crate for rendering with gpu support you should check out [Lyon](https://github.com/nical/lyon) from which I draw some inspiration, it's really good. It features lines, quadratic and cubic Béziers, in addition to arcs, triangles, and other geometric objects but no general Bézier curves or B-Splines. It also does seem to support wasm.  

Also, there's [Twinklebear/bspline](https://github.com/Twinklebear/bspline) which is a very clean and useful library for just bsplines. However, it depends on std and its simple Interpolate trait defines no way to access the individual dimensions of the points and hence implements no derivatives in the library.  

This clear online book [A Primer on Bézier Curves](https://pomax.github.io/Bézierinfo/) helped me with a lot of the algorithms involved.
