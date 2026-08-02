---
title: Complex coupled field residual
---

# Complex coupled field residual

`assemble_complex_coupled_field_residual` is the complex counterpart of the
real rectangular coupled residual.  For caller-owned blocks

\[
  r_f=A u-f, \qquad r_c=C u-g,
\]

`A` and `C` may have different row counts.  This is the small algebraic
boundary used to compose complex FEM, BEM, DtN, PML, Fourier, interface, or
wall blocks without prescribing a PDE or a global sparse storage layout.

The module also provides complete product-rule JVP and VJP routines.  VJPs
use the real-part complex inner product

\[
  \langle a,b\rangle=\operatorname{Re}\sum_i\overline{a_i}b_i,
\]

so a caller can combine the result with complex frequency-domain or linear
response operators without changing conventions.  The implementation only
uses caller-provided blocks and output arrays; it does not form a global
matrix or allocate a dense auxiliary system.

This layer is intentionally separate from
`assemble_linear_response_operator`, which owns the square
`K-omega^2 M+i omega R+V+W` block composition, and from the real
`assemble_coupled_field_residual`.  It is therefore suitable for a rectangular
field/constraint port such as an FEM volume paired with a BEM trace or a PML
field paired with wall/interface constraints.

The independent test checks the matrix-action oracle, a complex central/forward
JVP, the real-part adjoint identity, and rejection of incompatible dimensions.
