---
title: Signed cell-identification contract
---

# Signed cell-identification contract

`fortfem_cell_identification` records the discrete metadata needed when cells
are identified periodically or by a quotient construction. It deliberately
does not modify a mesh, infer a coordinate chart, or impose a physical
periodicity law.

## Contract

For every cell \(i\), the data contain

\[
r_i=\operatorname{representative}(i),
\qquad
s_i=\operatorname{orientation}(i)\in\{-1,+1\}.
\]

Representatives are one-based canonical IDs and must be idempotent:

\[
r_{r_i}=r_i.
\]

The canonical representative itself must have \(s_i=+1\). Noncanonical
members may have either orientation. This keeps periodic edge/face signs
explicit instead of silently dropping them during quotient assembly.

## API

```fortran
type(cell_identification_t) :: identification
integer, allocatable :: classes(:)
integer :: class_count, status

call initialize_cell_identification( &
    identification, representative, orientation, status)
call validate_cell_identification(identification, status)
call cell_identification_classes( &
    identification, classes, class_count, status)
```

`classes` contains compact class IDs in first-appearance order. The API is
safe for zero-cell metadata and rejects mismatched arrays, out-of-range
representatives, zero orientation, representative cycles, and a reversed
canonical representative.

## Scope and next layer

This is a discrete topology contract, not a quotient-mesh builder. A higher
layer must combine it with oriented cell-complex boundary maps, transform
incidences, and prove that the quotient still satisfies boundary-of-boundary
identities. Geometry maps, toroidal mode phases, trace bases, gauges, and
application-owned periodic boundary laws remain outside this primitive.

`test_cell_identification` supplies independent identity, signed-periodic,
cycle-rejection, canonical-orientation, zero-sign, and shape-mismatch oracles.
