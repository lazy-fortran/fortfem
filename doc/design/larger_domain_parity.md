# Larger-domain parity contract

`larger_domain_parity_t` is the neutral common-interior diagnostic for an
open-boundary control. A caller evaluates two complex FEM, BEM, DtN, PML, or
other backend fields at the same physical sample points, one with an inner
artificial-boundary distance and one with a farther boundary:

```fortran
call compare_larger_domain_solution(inner_field, outer_field, weights, &
    inner_contract, outer_contract, inner_distance, outer_distance, &
    absolute_tolerance, relative_tolerance, report, status)
```

The routine validates the same equation-space, units, normalization, topology,
and complex-field metadata as `compare_boundary_operator_parity`, while
allowing the backend kind to be the same at both distances. It reports the
positive weighted absolute field difference and a symmetric relative
difference, whose denominator is the larger of the two weighted field norms.
The distance increase, ratio, and relative-difference-per-distance are recorded
as convergence metadata. The latter is a two-domain finite-difference trend,
not an observed convergence order; an order requires a third distance or an
analytical reference.

No mesh, kernel, solver, reader, or application-specific boundary condition is
implemented here. The contract is intended to accompany FEM/BEM/DtN/PML
gallery and benchmark records, including toroidal external surfaces.
