# Shifted enriched space

`fortfem_shifted_enriched_space` is the matrix-level scalar XFEM/GFEM
composition used to build an enriched basis over a fixed cut topology.  For
base values (N_{i q}) at quadrature or plotting points (x_q), level-set
values (phi_q), and one anchor value (phi_i) per basis function, it
returns

\[
  N^{\mathrm{enr}}_{i q}
  = N_{i q}\bigl(H(\phi_q)-H(\phi_i)\bigr).
\]

The anchor signs provide the shifted partition-of-unity convention.  A zero
level-set or anchor value is rejected as a topology event.  Away from such an
event, the sign mask is locally constant, so the JVP differentiates the base
matrix and the VJP returns its exact reverse contraction while reporting zero
level/anchor derivatives.

The public procedures are:

```fortran
call evaluate_shifted_enriched_space( &
    base_values, level_values, anchor_values, enriched_values, status)
call evaluate_shifted_enriched_space_jvp( &
    base_values, level_values, anchor_values, base_dot, level_dot, &
    anchor_dot, enriched_dot, status)
call evaluate_shifted_enriched_space_vjp( &
    base_values, level_values, anchor_values, enriched_bar, base_bar, &
    level_bar, anchor_bar, status)
```

This is a space constructor, not a global mesh or quadrature builder.  Clients
own cut-cell integration, support activation, connected-component splitting,
Piola transformation, and sparse global numbering.  The same matrix contract
can therefore feed scalar IGA, fitted/unfitted FEM, or a later vector FEEC
constructor without importing application physics.

`test_shifted_enriched_space` checks the independently written sign formula,
central finite differences, the real adjoint identity, and explicit rejection
of a level-set topology event.
