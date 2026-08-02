# Complex packed N-field block graph residual

`fortfem_complex_block_graph_residual` is the frequency-domain counterpart
of the [real packed block graph](block_graph_residual.html).  It evaluates
the same fixed-topology sum of rectangular edge actions for complex fields,
so Helmholtz, curl--curl, FEM/BEM, DtN, PML, wall, and Fourier clients can
retain their natural packed blocks.

The primal and JVP use ordinary complex multiplication.  The VJP follows
FortFEM's real-part complex inner product:

\[
 \langle a,b\rangle_{\mathbb R}=\Re\sum_i\overline{a_i}b_i.
\]

Consequently, a block cotangent uses `residual_bar*conjg(state)` and a state
cotangent uses the conjugate-transpose block action.  This convention is the
same one used by the complex coupled-field and response APIs.

```fortran
call assemble_complex_block_graph_residual( &
    field_sizes, edge_rows, edge_columns, block_offsets, block_values, &
    state, rhs, residual, status)
call assemble_complex_block_graph_residual_jvp( &
    field_sizes, edge_rows, edge_columns, block_offsets, block_values, state, rhs, &
    block_values_dot, state_dot, rhs_dot, residual_dot, status)
call assemble_complex_block_graph_residual_vjp( &
    field_sizes, edge_rows, edge_columns, block_offsets, block_values, state, rhs, &
    residual_bar, block_values_bar, state_bar, rhs_bar, status)
```

`test_complex_block_graph_residual` independently materializes a small matrix
as an oracle, checks value and JVP reassembly, verifies the real-part adjoint
identity, and rejects inconsistent packed offsets.  No global matrix is
formed by the library path.
