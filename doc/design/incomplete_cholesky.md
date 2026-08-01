---
title: Incomplete-Cholesky contract
---

# Incomplete-Cholesky contract

`fortfem_incomplete_cholesky` provides a dense IC(0) factor for real,
symmetric positive-definite matrices. It preserves the lower-triangle
nonzero pattern of the input and introduces no fill. The apply operation
solves

\[
    L L^T z=r.
\]

The factor is a preconditioner, not a substitute for a compatible FEEC space
or a tree--cotree gauge. A curl--curl direct solve must remove its gauge
kernel first. For complex or nonsymmetric systems, use an explicitly selected
ILU-family path rather than silently applying IC.

## API

```fortran
type(incomplete_cholesky_factor_t) :: factor
call build_incomplete_cholesky(matrix, factor, status)
call apply_incomplete_cholesky(factor, rhs, solution, status)
```

`fortfem_sparse_incomplete_cholesky` provides the same fixed-pattern contract
directly on a FortSparse CSC matrix. It keeps only the lower CSC pattern,
performs sparse triangular application, and exposes right-hand-side JVP/VJP
actions:

```fortran
type(sparse_incomplete_cholesky_factor_t) :: sparse_factor
call build_sparse_incomplete_cholesky(csc_matrix, sparse_factor, status)
call apply_sparse_incomplete_cholesky(sparse_factor, rhs, solution, status)
call apply_sparse_incomplete_cholesky_jvp(sparse_factor, rhs_dot, solution_dot, status)
call apply_sparse_incomplete_cholesky_vjp(sparse_factor, solution_bar, rhs_bar, status)
```

The sparse builder rejects nonsymmetric or non-positive-pivot CSC inputs and
does not introduce fill. Factor construction is an inactive solver branch;
the fixed-factor actions are the differentiable contract.

`build_sparse_ichol` adds a deterministic drop tolerance and maximum number
of retained strict lower entries per column while reusing the same sparse
fixed-factor apply/JVP/VJP path:

```fortran
call build_sparse_ichol(csc_matrix, drop_tolerance, max_fill_per_column, &
    sparse_factor, status)
```

The reference builder performs its numeric Cholesky phase in a dense work
array, then stores only the selected lower CSC entries. The dense phase fixes
the semantics for an eventual scalable row-oriented implementation; the
preconditioner application itself remains sparse. The diagonal is always
retained, and non-positive or non-finite pivots are reported explicitly.

For nonsymmetric blocks, `fortfem_sparse_incomplete_lu` provides the matching
ILU(0) contract with strict-lower and upper CSC patterns:

```fortran
type(sparse_incomplete_lu_factor_t) :: ilu_factor
call build_sparse_incomplete_lu(csc_matrix, ilu_factor, status)
call apply_sparse_incomplete_lu(ilu_factor, rhs, solution, status)
call apply_sparse_incomplete_lu_jvp(ilu_factor, rhs_dot, solution_dot, status)
call apply_sparse_incomplete_lu_vjp(ilu_factor, solution_bar, rhs_bar, status)
```

The transpose VJP solves with the fixed `U^T L^T` factors. Zero pivots are
reported explicitly; the API does not substitute an SPD preconditioner for a
nonsymmetric response block. Drop tolerances and fill-controlled ILUT remain
higher-level solver options.

The same factor is selectable in the public PCG path with
`solver_options(preconditioner=ichol_preconditioner())` (the aliases `ic` and
`ic0` are accepted). The sparse baseline deliberately densifies a CSC matrix;
the standalone sparse IC(0) API above is available for callers that must keep
the factor sparse. `solver_options(preconditioner=ichol_controlled_preconditioner(),
drop_tolerance=..., fill_level=...)` selects the sparse controlled path in
`solve_sparse`; the exact full-fill one-step PCG fixture is independently
tested. Measured scaling and a row-oriented construction remain benchmark
work. The legacy BiCGSTAB kernel accepts only its ILU
factor contract and therefore reports an explicit unpreconditioned fallback
if ICHOL is selected there.

The builder rejects nonsymmetric, non-finite, non-square, and non-positive
pivot inputs. Factorization branches and sparsity patterns are fixed for
derivative calculations; a breakdown is a reported solver event. The
converged-state PCG JVP/VJP differentiates the exact linear solve, so the
inactive preconditioner iteration does not introduce a hidden derivative
dependency.

## Independent oracle

The focused tests use independent tridiagonal solve oracles for sparse IC(0),
controlled ICHOL, and ILU(0), check their fixed-factor adjoint identities, and
reject nonsymmetry, indefinite diagonals, invalid policies, and zero pivots
independently.
