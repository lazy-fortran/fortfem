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

Higher-order polygonal cells, moving connectivity, and genuinely curved
support-volume measures remain planned extensions.  The quadrilateral map is
the first reusable unstructured plane-cell contract and does not imply that a
cut or curved FCI topology is differentiable across an event.

## Provenance

The factorization follows the support-operator construction described in the
local PARALLAX FCI documentation and by
[Stegmeir et al.](https://doi.org/10.1016/j.cpc.2016.12.014).  FortFEM contains
an independent algebraic contract and no PARALLAX source or benchmark data.
