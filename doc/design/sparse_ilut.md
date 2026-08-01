---
title: Sparse ILUT contract
---

# Sparse ILUT contract

`build_sparse_ilut` provides a deterministic drop- and fill-controlled
incomplete-LU reference path for nonsymmetric FortSparse CSC matrices. It is
intended for response, transport, and other blocks for which an SPD
Cholesky preconditioner is not valid.

The builder first forms a complete numeric LU factorization in a dense work
array. It then retains the diagonal of \(U\), drops entries below the relative
threshold

\[
  \tau = \text{drop\_tolerance}\,\max(1,\|A\|_\infty),
\]

and keeps at most `max_fill_per_column` largest-magnitude strict entries in
each lower or upper column. Ties are resolved by the ascending row index
because selection scans rows in CSC order. The stored \(L\) and \(U\) factors
are sparse CSC arrays, with unit diagonal in \(L\). A zero or non-finite
pivot, invalid policy, or non-finite input returns an explicit FortSparse
status instead of producing a silently unstable preconditioner.

The dense construction is deliberate for this first neutral contract: it
fixes the reference factor semantics. `build_sparse_ilut_row` is the
memory-scalable companion: it converts CSC input to CSR once, performs the
no-pivot elimination against retained U rows, and stores only row work plus
the retained factor entries. Its drop and fill policy is applied as each row
is finalized, so zero-fill uses the no-fill row diagonal rather than the
complete-LU diagonal. Both builders export the same sparse CSC factor and
the same apply/JVP/VJP contract. The row builder uses O(n + nnz + nnz(L+U))
temporary storage instead of an n-by-n work array.

The apply path itself is sparse and has no factor rebuild.
`apply_sparse_ilut_jvp` and `apply_sparse_ilut_vjp` differentiate only the
active fixed-factor right-hand-side solve; fill selection and pivot decisions
are inactive topology/events.

The focused behavioral tests check a hand-derived four-by-four LU solution,
the dense-reference zero-fill diagonal limit, the row-builder no-fill
diagonal oracle, finite fill-controlled application, the real JVP/VJP
dot-product identity, and rejection of invalid policies and pivots. This
contract complements, rather than replaces, tree--cotree reduction for
gauge-singular curl--curl systems and IC(0) for SPD blocks.
