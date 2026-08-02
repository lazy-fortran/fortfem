---
title: Cell-complex cycle-space contract
---

# Cell-complex cycle-space contract

`cell_complex_cycle_basis` exposes the algebraic one-cycle space of an
oriented cell complex. For the edge coefficient matrix (D_1), its returned
columns (Z) satisfy

\[
    D_1 Z = 0.
\]

The basis is computed with scale-aware real elimination and is independent of
coordinates, metric Hodge factors, finite-element basis functions, and
application fields. The integer boundary matrices remain the source of
orientation truth; the real basis is only a convenient coordinate system for
subsequent cycle, flux, or nullspace operations.

The companion `cell_complex_cocycle_basis` computes the kernel of
(D_2^T). Its columns are edge cochains that annihilate every oriented face
boundary. Together, the two kernels provide the chain and cochain sides of
the first compatible-complex nullspace contract without selecting a metric or
a gauge.

`cell_complex_homology_cycle_basis` removes the span of the oriented face
boundaries from the cycle space and returns representatives of

\[
    H_1 = \ker(D_1) / \operatorname{im}(D_2).
\]

The representatives are selected by a deterministic rank-increase rule. They
retain the original edge ordering and orientation, but do not select a metric,
period normalization, cut, gauge, or physical magnetic flux.

## API

```fortran
call cell_complex_cycle_basis(complex, cycles, cycle_count, status)
call cell_complex_cocycle_basis(complex, cocycles, cocycle_count, status)
call cell_complex_homology_cycle_basis( &
    complex, homology_cycles, homology_count, status)
```

`cycles` and `cocycles` have one row per oriented edge and one column per
independent kernel vector. Zero-edge and acyclic complexes return correctly
shaped zero-width arrays.

## Independent oracle

The focused cell-complex test checks that an interval has no one-cycle, the
one-cell circle has a unit cycle and exact zero boundary, and the two-loop
torus CW complex has two independent cycles whose boundary vanishes. Betti
numbers are checked separately, so the cycle-space and homology dimensions are
not conflated. A filled triangular loop is additionally checked to have zero
first homology, while the torus returns two representatives. The same test
checks the circle and torus cocycle kernels against the independent
transpose-boundary oracle.
