---
title: Regularized surface-current layer
---

# Regularized surface-current layer

`evaluate_regularized_surface_current_layer` supplies the geometry-neutral
Gaussian approximation

\[
 \delta_\epsilon(d)=\frac{\exp[-(d/\epsilon)^2]}
 {\sqrt{\pi}\,\epsilon},\qquad
 \mathbf J_\epsilon(d)=\delta_\epsilon(d)\,\mathbf K.
\]

Here `d` is a caller-owned signed normal distance, `K` contains any number of
tangential sheet-current components, and `epsilon` is a positive physical
thickness. Because the primitive consumes distances rather than a particular
coordinate map, slab, cylinder, sphere, torus, cut-cell, and IGA callers share
the same kernel. It assigns no Ampere sign and contains no material, HKT, or
plasma closure.

The scalar contraction and its exact JVP/VJP are generated together by
FortSym. The public wrapper applies them sample by sample and accumulates the
distance and global-thickness cotangents. All primal, tangent, and cotangent
inputs must be finite; the thickness must be strictly positive. Fixed
signed-distance topology is differentiable, while a change in closest surface
branch or interface topology remains an explicit event outside this contract.

`evaluate_regularized_surface_current_integral` reports

\[
 N_h=\sum_q w_q\delta_\epsilon(d_q),\qquad
 \mathbf I_h=\sum_q w_q\mathbf J_\epsilon(d_q).
\]

For a resolved normal quadrature, `N_h` converges to the exact Gaussian
normalization \(\int_{-\infty}^{\infty}\delta_\epsilon(d)\,dd=1\). This gives
clients a direct resolution diagnostic and verifies convergence to the
integrated sheet current without importing a geometry-specific solver.

The focused test uses an independent analytical Gaussian oracle, a centered
finite difference for all inputs including thickness, the real adjoint dot
product identity, and direct quadrature sums for both diagnostics.
