---
title: Tensor diffusion matrix contraction
---

# Tensor diffusion matrix contraction

`assemble_tensor_diffusion_matrix` evaluates the fixed-topology two-dimensional
quadrature block

\[
A_{ij}=\sum_q w_q\,\nabla\varphi_i(q)^T K(q)\nabla\varphi_j(q).
\]

The compatibility routine has basis-gradient shape `(quadrature, basis, 2)`
and tensor shape `(quadrature, 2, 2)`. The new
`assemble_tensor_diffusion_matrix_3d` family uses the corresponding `3`
dimension, while `assemble_tensor_diffusion_matrix_nd` accepts any positive
spatial dimension with basis-gradient shape `(quadrature, basis, dimension)`
and tensor shape `(quadrature, dimension, dimension)`. Every quadrature weight
must be positive. All three entry points expose matching JVP and VJP routines;
the legacy 2D contract is unchanged.
The tensor is caller-owned and may be symmetric or nonsymmetric. This keeps
the block neutral for anisotropic conductivity, resistivity, elasticity, and
tensor-pressure clients; constitutive laws and global sparse ownership remain
outside the routine.

The JVP differentiates gradients, tensor entries, and weights. The VJP returns
the corresponding real transpose products. The focused
The `test_tensor_diffusion_matrix` fixture checks a hand-evaluated two-point
anisotropic matrix, a central-difference JVP, the real dot-product identity,
and invalid-weight rejection. The independent
`test_tensor_diffusion_matrix_3d` fixture checks the 3D and general-dimension
contracts, a hand-evaluated 3D matrix, central-difference JVP, transpose
identity, and strict dimension rejection.

For a global FEM assembly, callers evaluate this block on each cell and apply
their own orientation, numbering, boundary, and constraint maps. The same
contract can therefore be used with fitted, cut, DG, or IGA quadrature.
