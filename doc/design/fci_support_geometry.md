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
identity.  Curved/unstructured plane-cell measures remain planned extensions.

## Provenance

The factorization follows the support-operator construction described in the
local PARALLAX FCI documentation and by
[Stegmeir et al.](https://doi.org/10.1016/j.cpc.2016.12.014).  FortFEM contains
an independent algebraic contract and no PARALLAX source or benchmark data.
