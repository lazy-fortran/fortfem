---
title: Toroidal modal convolution
---

# Toroidal modal convolution

`apply_toroidal_modal_convolution` is a matrix-free retained-mode action for a
periodic scalar Green operator.  Given a fixed Fourier registry, it evaluates

\[
  y_k=\sum_\ell K_{k-\ell}x_\ell,
\]

including a term only when the difference mode is retained.  The kernel
coefficients, normalization, regularization, and zero-mode policy are
caller-owned.  Consequently the primitive is usable for NESTOR-like
regularized spectral Green actions without importing a toroidal equilibrium
or a particular surface convention.

JVP and VJP products cover both the modal kernel and source.  VJPs use
`Re(sum(conjg(bar)*dot))`, matching the other complex frequency-domain
contracts.  The implementation scans mode labels and deliberately does not
allocate an O(N²) difference map; this keeps auxiliary memory bounded for
large truncations.  A production client may replace the lookup with a sorted
or accelerator-backed map while preserving the same action contract.

The test compares the result with an independently written mode-label sum,
checks the product-rule derivative and verifies the complex adjoint identity.
