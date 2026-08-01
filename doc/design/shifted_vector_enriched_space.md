# Shifted vector enriched space

`fortfem_shifted_vector_enriched_space` lifts the scalar shifted-Heaviside
space constructor to a complete vector-valued basis matrix.  With component
(c), basis index (i), point (x_q), and fixed anchor (phi_i), it returns

\[
  N^{\mathrm{enr}}_{c i q}
  = N_{c i q}\bigl(H(\phi_q)-H(\phi_i)\bigr).
\]

The shared scalar sign mask keeps all vector components of one basis function
on the same side of the cut.  A zero level-set or anchor value is a topology
event and is rejected.  For a fixed sign pattern, the JVP differentiates the
base vector basis and the VJP contracts each component back to it; level and
anchor cotangents are zero.

```fortran
call evaluate_shifted_vector_enriched_space( &
    base_values, level_values, anchor_values, enriched_values, status)
call evaluate_shifted_vector_enriched_space_jvp( &
    base_values, level_values, anchor_values, base_dot, level_dot, &
    anchor_dot, enriched_dot, status)
call evaluate_shifted_vector_enriched_space_vjp( &
    base_values, level_values, anchor_values, enriched_bar, base_bar, &
    level_bar, anchor_bar, status)
```

The first axis is component, the second is basis function, and the third is
the quadrature or plotting point.  Piola maps, exact-sequence preservation,
cut integration, connected-component activation, and global sparse numbering
remain caller-owned layers.  Thus the same constructor can feed an H(curl),
H(div), vector IGA, or intentionally broken DG composition without selecting
one physical formulation.

`test_shifted_vector_enriched_space` checks an independent component/basis/
point sign formula, central finite differences, the real adjoint identity, and
explicit topology-event rejection.
