# Packed block graph to CSC

`assemble_block_graph_csc` is the explicit sparse boundary between the
matrix-free packed N-field graph and a retained FortSparse factorization. It
accepts the same field sizes, directed rectangular edges, offsets, and packed
column-major blocks as `assemble_block_graph_residual`, but emits either a real
`csc_t` or complex `csc_z_t` matrix through the generic interface.

The adapter allocates triplets proportional to the retained local block graph
and passes them to `csc_from_triplet`. Duplicate `(row,column)` entries are
therefore summed by FortSparse; no dense global matrix is formed. Matrix-free
callers should continue using the packed residual when a factor is not needed.

```fortran
call assemble_block_graph_csc( &
    field_sizes, edge_rows, edge_columns, block_offsets, block_values, matrix, status)
```

The real and complex paths share the same topology contract. The focused test
compares CSC matvecs against independently written dense oracles, checks
duplicate-edge compression through the resulting nonzero count, and rejects
inconsistent offsets. The resulting matrix can be passed directly to FortSparse
retained direct factors and their converged-state derivative actions. The
focused test exercises that complete real and complex graph-to-factor path,
including the JVP/VJP adjoint identity; solver policy, gauge reduction, and
block ordering remain caller-owned.
