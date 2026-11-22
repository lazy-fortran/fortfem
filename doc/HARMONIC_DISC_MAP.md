# Harmonic map onto the disc and pullback O-grids

Module `fortfem_harmonic_disc` (`src/discretisation/harmonic_disc.f90`),
test `test_harmonic_disc`.

## Purpose

A block-structured spline mesh of a plane domain bounded by a smooth simple
closed curve, valid without tuning for nonconvex and strongly shaped curves.
Ruled or transfinite block maps between a curved wall and a straight core
fold on concave walls; the construction here does not depend on the shape.

## Construction

G = (ξ, η) solves Δξ = Δη = 0 in the domain, with the curve mapped onto the
unit circle by normalised arc length from the first sample. The boundary
map is a homeomorphism onto the boundary of a convex set, so G is a
diffeomorphism of the domain onto the disc (Radó, Kneser, Choquet) whatever
the shape of the domain. The mesh is the pullback of the five-block O-grid
of the unit disc (`layout_ogrid`): the block maps are the tensor Greville
interpolants of G⁻¹ ∘ (exact ruled disc block map), so the mesh wall
interpolates the curve.

G needs no mesh. With z = x₁ + i x₂, ξ = Re v for a holomorphic v, taken as
the Cauchy integral of a real density μ on the curve (the double-layer
potential). On the curve, (1/2 + Re C) μ = boundary value, C the Cauchy
principal-value operator; its kernel is smooth for a smooth curve, and the
trapezoidal rule on samples that are equally spaced in a smooth periodic
parameter converges spectrally. The interior limit of v at the nodes follows
from μ after subtracting the singularity, and inside the domain v and v′ are
evaluated from these boundary values by the barycentric form of Cauchy's
formula, which keeps its accuracy up to the curve. G⁻¹ is Newton's method
from a polar table built by marching inward from the curve.

## Reported quantities (error control)

- `harmonic_disc_t%residual`: largest error of the boundary values of G at
  the nodes. It measures whether the samples resolve the curve; double the
  samples until it is small (1e-9 or less in the tests).
- `spline_jacobian_certificate`: the Jacobian determinant of a block map is,
  per cell, a polynomial of degree 2p − 1 per direction. Positive Bernstein
  coefficients on every cell prove det > 0 everywhere, not only at
  quadrature points. The condition is sufficient, not necessary; failed
  cells are counted, not subdivided.
- `harmonic_disc_wall_error`: largest distance of the mesh wall from the
  curve and the relative error of the enclosed area.
- `preserve_area` of `harmonic_disc_pull_back`: the wall data are moved
  along the curve normal by one constant (secant iteration) so that the mesh
  encloses the area of the curve. Volume integrals over the mesh then carry
  no first-order error from the wall fit.

## Limits

- A harmonic map spreads cells away from convex tips. The wall distance of
  the test shapes falls by 3.3 to 12.6 per halving between 8 and 16 cells
  per block side, below the spline order where a tip is not resolved. No
  grading of the boundary correspondence is implemented.
- Dense O(n³) solve and O(n) evaluation per point for n samples (256 to 512
  in use); meant for section meshes, not for large n.
- The first sample fixes the disc angle 0. A family of curves gives a smooth
  family of meshes when that point moves smoothly.

## Sources

- R. Kress, Linear Integral Equations (Springer): double-layer potential and
  the Nyström method with the trapezoidal rule on closed smooth curves.
- J. Helsing and R. Ojala, On the evaluation of layer potentials close to
  their sources, J. Comput. Phys. 227, 2899 (2008): boundary values of the
  Cauchy integral and evaluation near the boundary.
- R. Babin, F. Hindenlang, O. Maj and R. Köberl, Construction of an
  invertible mapping to boundary conforming coordinates for arbitrarily
  shaped toroidal domains, arXiv:2411.04683: harmonic map by boundary
  integrals, then inverted; cites Radó, Kneser and Choquet for
  invertibility. Idea only; their code was not read.
