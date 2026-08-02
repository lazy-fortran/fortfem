# Boundary operator contract

`boundary_operator_contract_t` is the typed metadata boundary shared by
interchangeable FEM, BEM, DtN, PML, NESTOR-like, BIEST-like, and virtual-casing
operators.  It records backend kind, equation, trace space, dimensions,
complexity, available matrix-free/assembled actions, derivative and residual
support, units, normalization, fixed-topology identity, and provenance.
It also records a neutral trace channel (`SCALAR`, `NORMAL`, `TANGENTIAL`, or
`MIXED`) and a caller-defined work-pairing label.  Clients can therefore
reject, for example, a normal port where a tangential port is required without
selecting a discretization or physical closure.

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

Existing callers keep the scalar-channel/L2 default.  A client that exposes a
different port opts in through the additive metadata initializer:

```fortran
call initialize_boundary_operator_trace_metadata( &
    contract, BOUNDARY_OPERATOR_TRACE_CHANNEL_TANGENTIAL, &
    "surface-tangential-work", status)
```

The channel constants are `BOUNDARY_OPERATOR_TRACE_CHANNEL_SCALAR`,
`BOUNDARY_OPERATOR_TRACE_CHANNEL_NORMAL`,
`BOUNDARY_OPERATOR_TRACE_CHANNEL_TANGENTIAL`, and
`BOUNDARY_OPERATOR_TRACE_CHANNEL_MIXED`.  `work_pairing` is an opaque,
nonempty label (for example, `surface-tangential-work`) whose interpretation
remains with the caller; clients should compare it when composing work ports.

Validation requires a known backend kind, positive dimensions, nonempty
equation/space/units/normalization/provenance/topology identifiers, at least
one representation (matrix-free or assembled), and all three derivative/
residual capabilities, a supported trace channel, and a nonempty work-pairing
label.  The last requirement is intentional: an operator may
be assembled or applied by a client, but it is not a complete FortFEM boundary
contract if its fixed-topology JVP/VJP or residual contribution is missing.

`test_boundary_operator_contract` checks these invariants independently,
including value-copy semantics and rejection of missing actions or provenance.
`test_boundary_operator_trace_metadata` additionally checks the default,
channel update, work-pairing validation, and rejection without mutation.
