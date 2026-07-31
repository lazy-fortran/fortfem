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

The builder rejects nonsymmetric, non-finite, non-square, and non-positive
pivot inputs. Factorization branches and sparsity patterns are fixed for
derivative calculations; a breakdown is a reported solver event.

## Independent oracle

The focused test uses a symmetric positive-definite tridiagonal matrix whose
IC(0) factor is exact and checks the residual after triangular application.
Nonsymmetry and an indefinite diagonal are rejected independently.

