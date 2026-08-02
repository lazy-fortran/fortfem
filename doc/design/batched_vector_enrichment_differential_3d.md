---
title: Batched vector enrichment differential
---

# Batched vector enrichment differential

`evaluate_batched_vector_enrichment_differential_3d` applies one fixed
partition-of-unity enrichment to every vector basis function at every
quadrature point:

\[
 b_e=\psi b,\qquad
 \nabla\times b_e=\psi\,\nabla\times b+\nabla\psi\times b,\qquad
 \nabla\cdot b_e=\psi\,\nabla\cdot b+\nabla\psi\cdot b.
\]

The batched layout is `(component,basis,point)` for vector quantities and
`(basis,point)` for scalar quantities. It is the missing matrix-free
composition between shifted vector XFEM/XIGA basis values and FEEC curl/div
diagnostics. Value, product-rule JVP, and real VJP actions are public.

The primitive does not decide Piola covariant/contravariant maps, cut-cell
quadrature, continuity, global numbering, or whether a physical jump should
be reported as an exact-sequence defect. Those are caller-owned layers and can
compose this contract with the existing FortSym-derived Piola products.
