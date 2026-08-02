---
title: Complex low-rank response contract
---

# Complex low-rank response contract

`complex_low_rank_matrix_t` is a neutral storage and action type for reusable
FEM/BEM/DtN/PML or conducting-wall response blocks. It stores factors in the
explicit convention

\[
  A \simeq U V^{T},
\]

not `V^H`; this matters for complex-symmetric frequency-domain operators.
`apply_complex_low_rank_matrix` evaluates the action without materializing a
dense matrix. `materialize_complex_low_rank_matrix` is retained as a small
dense oracle for reciprocity, passivity, and Schur tests.

`compress_complex_matrix_cross` uses deterministic greedy cross approximation
with a caller-supplied tolerance and rank cap. It records the scale-relative
unresolved residual as `relative_error_bound`; reaching the rank cap is not
silently reported as exact. The factor workspace is bounded before
allocation, preserving a safe small-oracle path for externally supplied
response data.

Factor and input perturbations have analytical JVP and real-complex VJP
actions. Pivot selection is a fixed-topology construction step and is not
differentiated. Production ACA, SVD, H-matrix, or FMM backends can populate
the same factor/action contract without adding application-specific physics or
file formats.
