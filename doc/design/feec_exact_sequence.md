---
title: FEEC exact-sequence diagnostics
---

# FEEC exact-sequence diagnostics

`assemble_feec_exact_sequence` is the metric-independent diagnostic shared by
simplicial FEEC, tensor-product IGA, multipatch quotients, and periodic cell
complexes. Given coefficient-space maps

\[
G:V^0\to V^1,\qquad C:V^1\to V^2,\qquad D:V^2\to V^3,
\]

it returns the algebraic defects `C G` and `D C`. Metric/Hodge matrices and
material tensors are intentionally absent: exactness is a topological
property and must be checked independently of quadrature and constitutive
assembly.

The `_jvp` routine differentiates both matrix products with respect to all
three maps. The `_vjp` routine returns the real reverse products, so a
FortSym-generated or hand-composed incidence map can participate in geometry
and parameter sensitivity without a hidden tape. Integer orientation and
quotient choices remain fixed topology, consistent with the tree--cotree and
multipatch contracts.

The focused test uses a triangle-like oriented cycle and a second map whose
rows are linearly dependent, verifies both exact identities with an independent
oracle, perturbs all maps for the JVP, and checks the real dot-product identity
for the VJP. Shape mismatches are rejected before any matrix product.
