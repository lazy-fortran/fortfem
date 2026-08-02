# Packed N-field block graph residual

`fortfem_block_graph_residual` generalizes the two-field block contract to an
arbitrary fixed field graph without building a monolithic matrix.  Field `i`
has `field_sizes(i)` unknowns.  Edge `e` maps `edge_columns(e)` to
`edge_rows(e)`, and its column-major rectangular block occupies

```text
block_values(block_offsets(e):block_offsets(e+1)-1)
```

The residual is the accumulated edge action

\[
 r_{I(i)} = -f_{I(i)}+
 \sum_{e:\,r(e)=i} A_e x_{I(c(e))}.
\]

Duplicate edges are allowed and are summed.  This is useful when element,
trace, FEM, BEM, DtN, PML, tensor, or Fourier contributions are retained in
their natural local ownership.  The graph is fixed and validated by the
caller; the module allocates no global matrix and does not choose a solver,
Schur complement, gauge, or application closure.

The API supplies the primal residual, product-rule JVP, and real VJP:

```fortran
call assemble_block_graph_residual( &
    field_sizes, edge_rows, edge_columns, block_offsets, block_values, &
    state, rhs, residual, status)
call assemble_block_graph_residual_jvp( &
    field_sizes, edge_rows, edge_columns, block_offsets, block_values, state, rhs, &
    block_values_dot, state_dot, rhs_dot, residual_dot, status)
call assemble_block_graph_residual_vjp( &
    field_sizes, edge_rows, edge_columns, block_offsets, block_values, state, rhs, &
    residual_bar, block_values_bar, state_bar, rhs_bar, status)
```

`test_block_graph_residual` independently materializes a small matrix only as
an oracle, checks a three-field graph with rectangular and duplicate edges,
finite-differences all differentiable inputs, verifies the adjoint identity,
and rejects inconsistent offsets.  Production callers can retain packed local
blocks and use the same derivative contract without the oracle matrix.
