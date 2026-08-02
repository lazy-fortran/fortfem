---
title: Eulerian non-nested residual contract
---

# Eulerian non-nested residual contract

`assemble_eulerian_nonnested_residual` is the closure-neutral composition
boundary for fixed-domain clients whose fields are not assumed to remain on
nested surfaces.  It concatenates caller-owned force and divergence residual
vectors and, when supplied, adds a stabilization vector:

\[
  R_E = \begin{bmatrix}R_{\rm force}\\R_{\rm div}\end{bmatrix}+R_{\rm stab}.
\]

The stabilization vector is deliberately precomputed.  A SIESTA-like client
can obtain it from `assemble_pseudo_transient_residual`, while another client
can provide a different fixed-topology continuation or preconditioner block.
This contract never chooses a mass matrix, relaxation closure, pressure law,
field representation, or free-boundary model.

Optional previous/current signed-margin arrays are passed to
`classify_continuation_event`.  The returned code, index, and minimum margin
are diagnostic metadata.  They are discrete and are therefore held fixed by
the JVP/VJP actions; a client must invalidate derivatives when a reported
event changes the topology.

## API

```fortran
call assemble_eulerian_nonnested_residual( &
    force_residual, divergence_residual, residual, status, &
    stabilization_residual, previous_margin, current_margin, event_tolerance, &
    event_code, event_index, minimum_margin)
call assemble_eulerian_nonnested_residual_jvp( &
    force_residual, divergence_residual, force_residual_dot, &
    divergence_residual_dot, residual_dot, status, stabilization_residual, &
    stabilization_residual_dot, previous_margin, current_margin, event_tolerance, &
    event_code, event_index, minimum_margin)
call assemble_eulerian_nonnested_residual_vjp( &
    force_residual, divergence_residual, residual_bar, force_residual_bar, &
    divergence_residual_bar, status, stabilization_residual, &
    stabilization_residual_bar)
```

The stabilization and topology arguments are optional as pairs.  If no
stabilization is supplied the concatenated force/divergence vector is returned.
If no margins are supplied, event code is `CONTINUATION_EVENT_NONE` and the
minimum margin is zero.  In a VJP, a stabilization cotangent is requested by
supplying `stabilization_residual_bar` together with its primal.

`test_eulerian_nonnested_residual` uses an independent concatenation formula,
central differences, and the real Euclidean dot-product identity.  It also
checks sign-crossing metadata, optional arguments, and rejection of malformed
stabilization and tolerance inputs.  No plasma closure or reader is involved.
