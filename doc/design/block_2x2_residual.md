# Two-field 2x2 block residual

`fortfem_block_2x2_residual` is the smallest neutral composition layer for
two independent state fields.  For caller-owned rectangular blocks it evaluates

\[
 r_1=A_{11}x+A_{12}y-f_1,\qquad
 r_2=A_{21}x+A_{22}y-f_2.
\]

The four blocks can be local FEM, BEM, DtN, PML, interface, tensor, or Fourier
actions.  The module does not materialize a global sparse matrix, select a
Schur complement, impose a gauge, or encode application physics.  Those
choices remain with the caller, while this contract makes composition and
derivatives uniform across mixed fields.

The public API provides the primal residual and product-rule real JVP/VJP:

```fortran
call assemble_block_2x2_residual( &
    block_11, block_12, block_21, block_22, state_1, state_2, rhs_1, rhs_2, &
    residual_1, residual_2, status)
call assemble_block_2x2_residual_jvp( &
    block_11, block_12, block_21, block_22, state_1, state_2, rhs_1, rhs_2, &
    block_11_dot, block_12_dot, block_21_dot, block_22_dot, state_1_dot, &
    state_2_dot, rhs_1_dot, rhs_2_dot, residual_1_dot, residual_2_dot, status)
call assemble_block_2x2_residual_vjp( &
    block_11, block_12, block_21, block_22, state_1, state_2, rhs_1, rhs_2, &
    residual_1_bar, residual_2_bar, block_11_bar, block_12_bar, block_21_bar, &
    block_22_bar, state_1_bar, state_2_bar, rhs_1_bar, rhs_2_bar, status)
```

For a fixed topology, the VJP is checked against

\[
 \langle\bar r_1,\dot r_1\rangle+\langle\bar r_2,\dot r_2\rangle
 =\sum_{i,j}\langle\bar A_{ij},\dot A_{ij}\rangle
  +\langle\bar x,\dot x\rangle+\langle\bar y,\dot y\rangle
  +\langle\bar f_1,\dot f_1\rangle+\langle\bar f_2,\dot f_2\rangle.
\]

`test_block_2x2_residual` uses an independently written matrix expression,
finite-differences every block/state/data direction, verifies this adjoint
identity, and checks incompatible block dimensions.  A future N-field or
sparse block graph can build on this contract without changing existing
single-field APIs.
