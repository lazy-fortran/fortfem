---
title: Oriented cell-complex contract
---

# Oriented cell-complex contract

`fortfem_cell_complex` is the neutral topological foundation for compatible
finite elements, gauges, cuts, and interface graphs. It does not read a mesh
file and it does not select a physical equation.

## Boundary convention

The stored integer matrices represent the chain maps

\[
C_3 \xrightarrow{\partial_3} C_2
    \xrightarrow{\partial_2} C_1
    \xrightarrow{\partial_1} C_0.
\]

Columns identify oriented cells and rows identify their oriented boundary
incidences. A valid complex therefore satisfies

\[
\partial_1\partial_2=0,
\qquad
\partial_2\partial_3=0.
\]

The matrices are integer data. The validation routine checks dimensions and
these identities exactly using a wider integer accumulator, so a reversed
edge or face cannot be hidden by a floating-point tolerance.

## API

```fortran
type(cell_complex_t) :: complex
integer :: status, euler, betti(4)

call initialize_cell_complex( &
    complex, vertex_count, boundary_1, boundary_2, boundary_3, status)
call validate_cell_complex(complex, status)
call cell_complex_euler_characteristic(complex, euler, status)
call cell_complex_betti_numbers(complex, betti, status)
```

`boundary_2` and `boundary_3` are optional. An absent map represents zero
cells in that dimension. Empty dimensions are stored as correctly shaped
zero-column arrays, allowing one-cell, circle, sphere, and torus CW fixtures
to use the same API.

The Betti calculation uses

\[
b_0=n_0-r_1,\quad b_1=n_1-r_1-r_2,\quad
b_2=n_2-r_2-r_3,\quad b_3=n_3-r_3,
\]

where \(r_k\) is the numerical rank of \(\partial_k\). The rank routine is a
small scale-aware elimination diagnostic, not a substitute for a sparse
homology package. It is deliberately sufficient for compact topology
fixtures and nullspace regression tests.

## Independent fixtures

`test_cell_complex` checks an interval, a loop, a one-vertex sphere CW model,
and the standard one-vertex/two-edge/one-face torus CW model against their
known Euler characteristics and Betti numbers. It also supplies a malformed
boundary map and verifies that a nonzero boundary of a boundary is rejected.

Periodic identifications, region adjacency, harmonic basis construction, and
cycle-ledger constraints remain higher-level contracts. They must consume this
orientation convention rather than introduce a second sign convention.
