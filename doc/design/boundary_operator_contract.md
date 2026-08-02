# Boundary operator contract

`boundary_operator_contract_t` is the typed metadata boundary shared by
interchangeable FEM, BEM, DtN, PML, NESTOR-like, BIEST-like, and virtual-casing
operators.  It records backend kind, equation, trace space, dimensions,
complexity, available matrix-free/assembled actions, derivative and residual
support, units, normalization, fixed-topology identity, and provenance.

The record deliberately contains no procedure pointers and no application
file-format data.  Numerical actions remain caller-owned public APIs, while
this descriptor lets a free-boundary or open-boundary client reject an
incompatible block before composing it.

```fortran
call initialize_boundary_operator_contract( &
    contract, BOUNDARY_OPERATOR_BACKEND_BEM, "helmholtz", "H1-trace", &
    row_count, column_count, complex_valued, matrix_free, assembled, &
    has_jvp, has_vjp, has_residual, units, normalization, provenance, &
    topology_id, status)
valid = validate_boundary_operator_contract(contract, status)
```

Validation requires a known backend kind, positive dimensions, nonempty
equation/space/units/normalization/provenance/topology identifiers, at least
one representation (matrix-free or assembled), and all three derivative/
residual capabilities.  The last requirement is intentional: an operator may
be assembled or applied by a client, but it is not a complete FortFEM boundary
contract if its fixed-topology JVP/VJP or residual contribution is missing.

`test_boundary_operator_contract` checks these invariants independently,
including value-copy semantics and rejection of missing actions or provenance.
