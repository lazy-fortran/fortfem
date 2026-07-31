---
title: Harmonic period normalization
---

# Harmonic period normalization

`normalize_harmonic_one_forms` converts a fixed-topology harmonic basis into
caller-selected period or flux units. Let `H` contain edge one-forms and `C`
contain oriented cycle representatives. The cycle-period matrix is

\[
  P=C^T H.
\]

For a target period matrix `T`, the routine solves

\[
  P A=T,
  \qquad H_{\rm normalized}=H A.
\]

Thus `C^T*H_normalized=T` without choosing physical flux units, gauge laws, or
cycle labels. The topology is fixed during differentiation; a cycle or basis
rebuild is a discrete topology event.

```fortran
call normalize_harmonic_one_forms( &
    harmonic_forms, period_cycles, target_periods, normalized_forms, &
    normalization_matrix, status)
call normalize_harmonic_one_forms_jvp( &
    harmonic_forms, period_cycles, target_periods, harmonic_forms_dot, &
    period_cycles_dot, target_periods_dot, normalized_forms_dot, &
    normalization_matrix_dot, status)
call normalize_harmonic_one_forms_vjp( &
    harmonic_forms, period_cycles, target_periods, normalized_forms, &
    normalization_matrix, normalized_forms_bar, normalization_matrix_bar, &
    harmonic_forms_bar, period_cycles_bar, target_periods_bar, status)
```

The JVP differentiates both the cycle-period matrix and the dense solve. The
VJP uses the transpose period solve and returns cotangents for the harmonic
basis, cycle representatives, and target periods. A singular period matrix is
reported rather than silently changing the topological basis.

`test_harmonic_period_normalization` compares the primal map with an explicit
two-by-two inverse oracle, checks the JVP against central differences, verifies
the real VJP dot identity, and rejects incompatible output shapes.
