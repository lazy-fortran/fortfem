# Fixed/free-boundary external-field port

`assemble_free_boundary_port_residual` is the neutral physical-space port for
an already evaluated trace.  For quadrature sample `q`, trace component `c`,
positive work weight `w`, caller-owned trace `t`, and supplied exterior/vacuum
target `g`, it evaluates

\[
  r_{qc}=w_q\left(t_{qc}-g_{qc}-k_{qc}\right).
\]

`k` is optional and denotes a supplied sheet-current or jump target expressed
in the same work-conjugate trace coordinates.  FortFEM does not project a
surface current, select an Ampere sign, construct a BEM/DtN/PML map, or choose
an equilibrium/free-boundary law.  An external FEM, BEM, DtN, NESTOR-like,
virtual-casing, wall, or equilibrium adapter owns those decisions and passes
the resulting physical samples to this contract.  Omitting `k` gives the
fixed-boundary/exterior target port; supplying it gives the same residual
shape convention for a fixed-topology free-boundary jump.

The three allocation-free actions are:

```fortran
call assemble_free_boundary_port_residual( &
    physical_trace, external_target, weights, residual, status, sheet_current)
call assemble_free_boundary_port_residual_jvp( &
    physical_trace, external_target, weights, physical_trace_dot, &
    external_target_dot, weights_dot, residual_dot, status, sheet_current, &
    sheet_current_dot)
call assemble_free_boundary_port_residual_vjp( &
    physical_trace, external_target, weights, residual_bar, &
    physical_trace_bar, external_target_bar, weights_bar, status, &
    sheet_current, sheet_current_bar)
```

The JVP includes trace, target, sheet, and weight product rules.  The VJP is
the real transpose under the Euclidean sample/component pairing.  A sheet
target and its derivative/cotangent must be either both supplied or both
omitted.  The base arrays must have a common `(n_quadrature,n_component)`
shape, weights must have `n_quadrature` strictly positive finite entries, and
all data and outputs must be finite and correctly shaped.  Zero-row ports are
rejected explicitly so a missing port cannot silently masquerade as a valid
boundary condition.

This layer deliberately reuses the same fixed-topology trace/work conventions
as `assemble_boundary_trace_residual` and
`assemble_complex_boundary_trace_residual`, but it does not duplicate their
trace-map action: here the physical trace has already been sampled by the
caller.  The independent test uses a direct manufactured component-wise oracle,
central-difference JVP checks, the real dot-product VJP identity, omitted-sheet
coverage, and invalid/zero-row cases.

## Provenance and scope

The contract is the algebraic boundary-port part of the fixed/free-boundary
composition described in [FortFEM's roadmap](../../ROADMAP.md).  It is
compatible with work-conjugate trace coupling used by FEM--BEM and
FEM--DtN formulations, but it intentionally contains no kernel, quadrature,
surface reader, coil model, pressure law, or production equilibrium code.
Those remain external, license-safe application responsibilities.
