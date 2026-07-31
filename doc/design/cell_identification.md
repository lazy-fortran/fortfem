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
call identify_boundary_matrix( &
    lower_identification, upper_identification, boundary, &
    quotient_boundary, status)
```

`classes` contains compact class IDs in first-appearance order. The API is
safe for zero-cell metadata and rejects mismatched arrays, out-of-range
representatives, zero orientation, representative cycles, and a reversed
canonical representative.

`identify_boundary_matrix` pushes an oriented boundary map to quotient
classes. For a lower-cell identification \((r_i,s_i)\), each row contributes
\(s_i B_{ij}\) to its representative row. The canonical column of every upper
class defines the quotient column. All other identified upper columns are then
checked exactly against their declared orientation; inconsistent periodic data
are rejected instead of being averaged.

## Scope and next layer

This remains a discrete topology contract, not a quotient-mesh builder. The
`quotient_cell_complex` layer combines it with oriented boundary maps and
checks the resulting boundary-of-boundary identities. Geometry maps, toroidal
mode phases, trace bases, gauges, and application-owned periodic boundary
laws remain outside both topology primitives.

`test_cell_identification` supplies independent identity, signed-periodic,
interval-to-circle, signed-column, inconsistency-rejection,
cycle-rejection, canonical-orientation, zero-sign, and shape-mismatch oracles.
