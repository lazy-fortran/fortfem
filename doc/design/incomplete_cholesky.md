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
array, then stores only the selected lower CSC entries. The public
`build_sparse_ichol_row` companion computes each row from previously retained
lower rows and uses O(n + nnz + nnz(L)) temporary storage; its final factor
has the same CSC apply/JVP/VJP contract. The row path evaluates each diagonal
before dropping strict lower entries, so its zero-fill limit is an independent
Cholesky diagonal oracle. The diagonal is always retained, and non-positive or
non-finite pivots are reported explicitly.

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
higher-level solver options; the sparse ILUT module also exposes a
memory-scalable row-oriented builder.

The same factor is selectable in the public PCG path with
`solver_options(preconditioner=ichol_preconditioner())` (the aliases `ic` and
`ic0` are accepted). `solve_sparse` keeps this path sparse and uses the CSC
IC builder directly; `ichol_controlled_preconditioner()` selects the explicit
drop/fill-controlled spelling. The exact full-fill one-step PCG fixture is
independently tested. The row-oriented ICHOL constructor has an independent
tridiagonal oracle and fixed-factor adjoint check. Measured scaling remains
benchmark work. Sparse BiCGSTAB and GMRES use their matrix-free callback
paths; the dense BiCGSTAB entry point retains its legacy dense-ILU contract.

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
