# Coupled field residual

`fortfem_coupled_field_residual` is the neutral composition layer for a
caller-owned collection of finite-element, boundary, Fourier, tensor, or
interface blocks.  Given a field operator (A), a state (u), field data
(f), a constraint operator (C), and target data (g), it evaluates

\[
  r_f = A u - f, \qquad r_c = C u - g.
\]

The field and constraint operators may be rectangular and may have different
row counts.  FortFEM does not decide how their rows are assembled or how
global sparse ownership is stored; clients can use this contract after
composing compatible H1, H(curl), H(div), DG, IGA, BEM, DtN, PML, tensor, or
surface-current blocks.

The public API provides the primal residual and the exact real JVP/VJP:

```fortran
call assemble_coupled_field_residual( &
    field_operator, state, source, constraint_operator, constraint_target, &
    field_residual, constraint_residual, status)
call assemble_coupled_field_residual_jvp( &
    field_operator, state, source, constraint_operator, constraint_target, &
    field_operator_dot, state_dot, source_dot, constraint_operator_dot, &
    constraint_target_dot, field_residual_dot, constraint_residual_dot, status)
call assemble_coupled_field_residual_vjp( &
    field_operator, state, source, constraint_operator, constraint_target, &
    field_residual_bar, constraint_residual_bar, field_operator_bar, state_bar, &
    source_bar, constraint_operator_bar, constraint_target_bar, status)
```

For a fixed topology, the VJP satisfies the real inner-product identity

\[
  \langle \bar r_f, \dot r_f\rangle+
  \langle \bar r_c, \dot r_c\rangle
  = \langle \bar A,\dot A\rangle+
    \langle \bar C,\dot C\rangle+
    \langle \bar u,\dot u\rangle+
    \langle \bar f,\dot f\rangle+
    \langle \bar g,\dot g\rangle.
\]

`test_coupled_field_residual` checks the independent matrix expression, a
central finite difference of every input, the adjoint identity, and shape
rejection.  The residual is intentionally agnostic about application physics:
equilibrium selection, constitutive laws, gauges, boundary laws, and sparse
global assembly remain caller-owned.
