---
title: Generalized complex eigen-residual
---

# Generalized complex eigen-residual

FortFEM exposes a neutral residual contract for modal and linear-response
clients:

\[
  r(K,M,\lambda,u) = K u - \lambda M u.
\]

`assemble_generalized_eigen_residual` evaluates this expression without
choosing an eigensolver, normalization, coordinate system, or application
model.  `K` and `M` are caller-owned complex square matrices, so the same
contract can be used for a Fourier-FEM block, a reduced curl--curl system, a
resistive linearization, or an externally supplied oracle.

The tangent action is

\[
\begin{aligned}
Dr[\dot K,\dot M,\dot\lambda,\dot u] ={}&
 \dot K u + K\dot u - \dot\lambda M u\\
 &- \lambda\dot M u - \lambda M\dot u.
\end{aligned}
\]

The VJP uses the real part of the complex inner product.  For a residual
cotangent `r_bar`, it returns

\[
\begin{aligned}
\bar K_{ij} &= \bar r_i\,\overline{u_j},\\
\bar M_{ij} &= -\overline{\lambda}\,\bar r_i\,\overline{u_j},\\
\bar u &= K^H\bar r-\overline{\lambda}M^H\bar r,\\
\bar\lambda &= -\sum_i \bar r_i\,\overline{(Mu)_i}.
\end{aligned}
\]

The test uses an independent matrix-product oracle, a forward-difference
check for the JVP, and the real-complex adjoint identity for the VJP.  Shape
and finite-value validation returns a nonzero status instead of allowing an
invalid modal block into an external adapter.

This is intentionally a foundation contract.  FortFEM does not read GLISS,
GPEC, MARS-F, VMEC, or equilibrium files here; application adapters remain
responsible for their own metadata, normalization, and physical operators.
