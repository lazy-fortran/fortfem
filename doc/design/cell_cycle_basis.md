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

This routine deliberately returns cycles, not homology representatives or
metric-harmonic forms. Face boundaries, period normalization, cuts, gauges,
and Hodge choices must be supplied by the higher FEEC/topology layer. In
particular, a nonzero cycle count is not by itself a claim about a physical
magnetic flux or a plasma equilibrium degree of freedom.

## API

```fortran
call cell_complex_cycle_basis(complex, cycles, cycle_count, status)
```

`cycles` has one row per oriented edge and one column per independent cycle.
Zero-edge and acyclic complexes return correctly shaped zero-width arrays.

## Independent oracle

The focused cell-complex test checks that an interval has no one-cycle, the
one-cell circle has a unit cycle and exact zero boundary, and the two-loop
torus CW complex has two independent cycles whose boundary vanishes. Betti
numbers are checked separately, so the cycle-space and homology dimensions are
not conflated.

