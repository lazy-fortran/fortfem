---
title: FCI support-volume geometry
---

# FCI support-volume geometry

`compute_fci_staggered_flux_box_volumes` is the geometry-side contract for the
positive staggered volumes used by the FCI support divergence.  It combines
the two field-line flux-expansion integrals with the plane-cell area and local
toroidal field:

\[
  \omega_\mu = (v^+_\mu+v^-_\mu)\,A_\mu\,B_{\varphi,\mu}.
\]

The routine deliberately accepts already traced factors.  It does not know
the equilibrium, mesh, or MPI layout, and it rejects non-finite, non-positive,
or mismatched inputs.  This keeps the volume construction reusable by the
matrix-free and assembled support-operator paths.

The value, JVP, and VJP products are emitted by the pinned FortSym generator
`gen_fci_support_volume_products`; the focused test compares the value and JVP
against independent product-rule oracles and checks the VJP dot-product
identity.

For unstructured planar cells,
`compute_fci_quadrilateral_cell_areas_2d` accepts
`cell_vertices(2,4,ncell)` in counter-clockwise boundary order and returns a
positive shoelace measure.  Repeated vertices, self-intersections, clockwise
orientation, degeneracy, non-finite coordinates, and shape mismatches are
reported before a measure is used.  The companion `_jvp` and `_vjp` routines
differentiate the fixed-topology map with respect to all vertex coordinates.
Their kernels are emitted by the pinned `gen_fci_quadrilateral_area_products`
FortSym generator.  The independent test uses a separately written shoelace
oracle, central differences, and the real dot-product identity; the gallery
fixture renders the polygons and their areas.

The generic `compute_fci_polygon_cell_areas_2d` contract now accepts
`cell_vertices(2,nvertex,ncell)` with any fixed `nvertex >= 3`.  It sums a
FortSym-generated per-edge Green contribution, validates repeated vertices
and non-adjacent edge intersections, and exposes matching `_jvp` and `_vjp`
actions.  Shared endpoint cotangents are accumulated across the two incident
edges.  The pentagon test and gallery fixture use independent shoelace
oracles, central differences, and a real dot-product identity.  Moving
connectivity and genuinely curved support-volume measures remain planned
extensions.  The curved quadrilateral map is a fixed-topology quadratic
Bezier-edge contract: its vertex ordering and finite control data are
validated, while callers remain responsible for keeping the curved boundary
non-self-intersecting.  Neither straight nor curved FCI topology is
differentiable across a connectivity or orientation event.

For edge `e`, with endpoints `p_e`, `p_{e+1}` and control point `c_e`, the
boundary is

\[
  r_e(t)=(1-t)^2p_e+2(1-t)t c_e+t^2p_{e+1},\qquad 0\leq t\leq1,
\]

and the generated measure evaluates the Green line integral
`1/2 integral (x dy - y dx)` exactly for these quadratic edges.  The
independent test evaluates the same line integral with three-point Gauss
quadrature, then checks central differences and the real VJP dot-product
identity.  The gallery samples the actual Bezier boundaries rather than
replacing them by a polygonal preview.

## Provenance

The factorization follows the support-operator construction described in the
local PARALLAX FCI documentation and by
[Stegmeir et al.](https://doi.org/10.1016/j.cpc.2016.12.014).  FortFEM contains
an independent algebraic contract and no PARALLAX source or benchmark data.
