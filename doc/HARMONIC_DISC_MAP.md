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
unit circle by normalised arc length from the first sample (or by a
curvature-weighted measure, "Graded boundary correspondence" below). The
boundary map is a homeomorphism onto the boundary of a convex set, so G is a
diffeomorphism of the domain onto the disc (Radó, Kneser, Choquet) whatever
the shape of the domain. The mesh is the pullback of the five-block O-grid
of the unit disc (`layout_ogrid`): the block maps are the tensor Greville
interpolants of G⁻¹ ∘ (exact ruled disc block map).

The wall of the mesh is one periodic spline around the whole curve. Its
breaks are the cell breaks of the wall sides of the four outer blocks (as
disc angles), it interpolates the curve at its Greville angles, and the wall
data of a block are its values at the block's Greville points. On a block
the periodic spline lies in the spline space of that block side, so the
block interpolation reproduces it: the wall is C^(p−1) at the block corners
as inside a block. Interpolating the curve block by block (the first
version) made the two wall splines that meet at a block corner end with
different tangents. On a coarse mesh of a strongly shaped section that is a
kink of several degrees (up to 11.6° on the planes of a quasi-axisymmetric
stellarator section at 4 cells per block side, degree 3), a corner of the
domain that the curve does not have; a consumer found it as a local spike of
the current density of an equilibrium solved on the mesh (MILCA,
docs/flux-map.md, "Wall layer of the force-balance residual").
`test_harmonic_disc` checks the continuity of the wall tangent (every
degree) and of the wall curvature (degree 3 and above) at the block corners.

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

## Graded boundary correspondence (optional)

`harmonic_disc_build(..., grading=c)` replaces the arc length l in the
boundary correspondence by the measure

    dω = (1 + (c κ L/2π)²)^(1/3) dl,

κ the curvature and L the length of the curve (κL/2π = 1 on a circle, which
the grading leaves unchanged). The disc angle of a curve point is 2π ω/ω(L),
so equal arcs of the circle, and with them the wall cells of a pulled-back
mesh, hold equal amounts of ω: more cells where the curve bends. The weight
is a smooth function of the curve (no |κ| with its kink at inflection
points), so the trapezoidal rule keeps its spectral accuracy; it grows like
κ^(2/3) at a tip. The exponent is a choice between the densities κ^(1/2)
(equal chord error) and κ (equal turning angle per cell); it is not derived
from an error estimate of the splines.
Any such correspondence is a homeomorphism onto the circle, so the theorem
of Radó, Kneser and Choquet and the certificate apply unchanged. c = 0 (the
default) is the arc length, bit for bit.

Measured (`test_harmonic_disc`, 4 cells per block side, degree 2, c = 1,
256 samples, area preserved): the wall distance falls from 5.45e-2 to
1.70e-2 on the bean and from 3.30e-2 to 1.12e-2 on the rounded triangle;
the Bernstein ratio rises from 0.082 to 0.184 and from 0.103 to 0.199; the
boundary residual stays at 1e-12. The correspondence then depends on the
curvature: for a family of curves whose tips move (the sections of a
stellarator along the torus) the meshes slide along the wall with the
family parameter and need more resolution in that parameter (a consumer's
measurement: MILCA, docs/flux-map.md, "Mesh grading towards tips").

## Reported quantities (error control)

- `harmonic_disc_t%residual`: largest error of the boundary values of G at
  the nodes. It measures whether the samples resolve the curve; double the
  samples until it is small (1e-9 or less in the tests).
- `spline_jacobian_certificate`: the Jacobian determinant of a block map is,
  per cell, a polynomial of degree 2p − 1 per direction. Positive Bernstein
  coefficients on every cell prove det > 0 everywhere, not only at
  quadrature points. The condition is sufficient, not necessary; failed
  cells are counted, not subdivided. The generated interpolation tables
  cover degrees 1 to 7.
- `harmonic_disc_wall_error`: largest distance of the mesh wall from the
  curve and the relative error of the enclosed area.
- `preserve_area` of `harmonic_disc_pull_back`: the wall data are moved
  along the curve normal by one constant (secant iteration) so that the mesh
  encloses the area of the curve. Volume integrals over the mesh then carry
  no first-order error from the wall fit.

## Limits

- A harmonic map spreads cells away from convex tips. The wall distance of
  the test shapes falls by 3.3 to 13.5 per halving between 8 and 16 cells
  per block side, below the spline order where a tip is not resolved. The
  graded boundary correspondence above places the wall cells at the tips;
  the interior cells follow the harmonic map.
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
