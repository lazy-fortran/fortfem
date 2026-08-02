# Retained field-split solves

`fortfem_retained_field_split` composes independent field blocks into one
retained block-diagonal solve. Each field owns a square real `csc_t` or complex
`csc_z_t` matrix, a direct factor, and a transpose/adjoint factor for reverse
products. The blocks are concatenated in caller order; no coupled off-diagonal
physics is assumed or silently discarded from a residual.

```fortran
call factor_retained_field_split(field_matrices, split, status)
call apply_retained_field_split(split, rhs, solution, status)
call apply_retained_field_split_jvp( &
    split, field_matrix_dots, solution, rhs_dot, solution_dot, status)
call apply_retained_field_split_vjp( &
    split, solution, solution_bar, rhs_bar, field_matrix_bars, status)
```

The `complex` counterparts use the real-part complex adjoint convention. The
JVP/VJP paths reuse FortSparse's fixed-factor solve products, so repeated
field-split iterations do not refactor. A coupled Schur correction, Krylov
policy, gauge reduction, and field ordering remain caller-owned and can be
composed with `fortfem_block_graph_residual` or the local HDG/wall Schur
primitives.

`test_retained_field_split` uses independent one- and two-by-two block oracles,
checks real and complex solves, and verifies both real and complex adjoint
identities. Call `free_retained_field_split` (or its complex counterpart) when
the retained factors are no longer needed.
