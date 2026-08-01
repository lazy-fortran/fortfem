---
title: Tensor diffusion matrix contraction
---

# Tensor diffusion matrix contraction

`assemble_tensor_diffusion_matrix` evaluates the fixed-topology quadrature
block

\[
A_{ij}=\sum_q w_q\,\nabla\varphi_i(q)^T K(q)\nabla\varphi_j(q).
\]

The basis-gradient array has shape `(quadrature, basis, 2)`, the tensor array
has shape `(quadrature, 2, 2)`, and every quadrature weight must be positive.
The tensor is caller-owned and may be symmetric or nonsymmetric. This keeps
the block neutral for anisotropic conductivity, resistivity, elasticity, and
tensor-pressure clients; constitutive laws and global sparse ownership remain
outside the routine.

The JVP differentiates gradients, tensor entries, and weights. The VJP returns
the corresponding real transpose products. The focused
`test_tensor_diffusion_matrix` fixture checks a hand-evaluated two-point
anisotropic matrix, a central-difference JVP, the real dot-product identity,
and invalid-weight rejection.

For a global FEM assembly, callers evaluate this block on each cell and apply
their own orientation, numbering, boundary, and constraint maps. The same
contract can therefore be used with fitted, cut, DG, or IGA quadrature.
