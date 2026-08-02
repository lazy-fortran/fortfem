# Boundary-operator parity contract

`boundary_operator_parity_t` is a neutral validation layer for comparing
exterior-boundary clients on one physical sample set. A manufactured complex
reference field and candidate values from FEM, BEM, DtN, PML, or another
declared backend are reduced with the same positive sample weights. The
report records the common topology/provenance, weighted reference norm, and
absolute/relative error for every backend.

The routine performs no assembly, solve, mesh import, or application-specific
normalization:

```fortran
call compare_boundary_operator_parity(reference, candidates, weights, &
    contracts, absolute_tolerance, relative_tolerance, report, status)
```

`contracts` are the existing `boundary_operator_contract_t` records. The
parity routine rejects invalid metadata, duplicate backends, mixed equation or
trace spaces, mixed units/normalizations, and mixed fixed-topology identifiers.
`candidates(component, sample, backend)` and `reference(component, sample)`
must be finite and use the same physical samples. A backend passes when its
absolute or relative weighted error is below the declared tolerance. This
gives open-boundary examples a reproducible FEM/BEM/DtN/PML parity checkpoint
without coupling FortFEM to an external code or file format.

`test_boundary_operator_parity` uses a hand-computed weighted complex oracle,
checks pass/fail classification, and verifies rejection of a mixed topology.
