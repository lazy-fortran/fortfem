# Real boundary trace residual

`assemble_boundary_trace_residual` is the real counterpart of the complex
boundary-port contract.  For caller-owned normal and tangential trace maps
`N` and `T`, positive work weights, state `u`, and supplied targets it
evaluates

\[
 r_n=w_n(Nu-g_n),\qquad r_t=w_t(Tu-g_t).
\]

The targets may be prescribed scalar traces, jumps, surface-current data, or
total-pressure data. FEM, BEM, DtN, PML, IGA, DG, and wall clients can supply
the maps without selecting a constitutive law or equilibrium convention.

The value, product-rule JVP, and real VJP include both maps, weights, state,
and targets:

```fortran
call assemble_boundary_trace_residual( &
    normal_trace, tangential_trace, normal_weights, tangential_weights, state, &
    normal_target, tangential_target, normal_residual, tangential_residual, status)
call assemble_boundary_trace_residual_jvp( &
    normal_trace, tangential_trace, normal_weights, tangential_weights, state, &
    normal_target, tangential_target, normal_trace_dot, tangential_trace_dot, &
    normal_weights_dot, tangential_weights_dot, state_dot, normal_target_dot, &
    tangential_target_dot, normal_residual_dot, tangential_residual_dot, status)
call assemble_boundary_trace_residual_vjp( &
    normal_trace, tangential_trace, normal_weights, tangential_weights, state, &
    normal_target, tangential_target, normal_residual_bar, tangential_residual_bar, &
    normal_trace_bar, tangential_trace_bar, normal_weights_bar, tangential_weights_bar, &
    state_bar, normal_target_bar, tangential_target_bar, status)
```

The implementation uses no per-call allocatable work arrays.  The independent
test checks direct matrix expressions, finite differences of every input, the
real adjoint identity, and invalid-shape/weight behavior through the shared
status contract.
