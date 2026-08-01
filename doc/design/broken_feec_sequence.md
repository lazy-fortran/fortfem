---
title: Broken FEEC sequence
---

# Broken FEEC sequence

`fortfem_broken_feec_sequence` embeds caller-owned cell-local de Rham maps in
independent cell-major blocks. For local maps

\[
 G_K: V^0_K\to V^1_K,\qquad
 C_K: V^1_K\to V^2_K,\qquad
 D_K: V^2_K\to V^3_K,
\]

the public assembly returns

\[
 G=\operatorname{diag}(G_K),\qquad
 C=\operatorname{diag}(C_K),\qquad
 D=\operatorname{diag}(D_K).
\]

The corresponding compositions are exposed as `C*G` and `D*C`. Thus a
cell-local exact sequence remains exact after embedding, while no continuity
or interface law is silently introduced between cells. This is the neutral
block used by DG, HDG, cut/XFEM, and patch-local IGA clients before they add
their own trace, penalty, mortar, or metric maps.

The JVP differentiates every local block and both products. The VJP includes
the direct block cotangents and the reverse products from both composition
diagnostics. Cell count and all local ranks are fixed topology; malformed or
non-finite maps are rejected before output is used.

The independent test checks cell-local exactness, absence of off-block
coupling, finite-difference JVP agreement, the real dot-product VJP identity,
and inconsistent-cell rejection.
