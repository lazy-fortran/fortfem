---
title: Pseudo-arclength continuation residual contract
---

# Pseudo-arclength continuation residual contract

`assemble_pseudo_arclength_residual` is the neutral continuation layer for a
caller-owned nonlinear residual. It appends one scalar constraint to
`R(state, parameter)` using a supplied predictor and tangent:

\[
  \widehat R = \begin{bmatrix}R\\
  t_u^T(u-u_0)+t_p^T(p-p_0)-\Delta s\end{bmatrix}.
\]

The routine does not choose a parameter, compute a tangent, perform a Newton
step, or own a free-boundary model. Those decisions remain with the client.
This makes the same contract usable for equilibrium continuation, response
curves, nonlinear elasticity, and structure-preserving wave problems.

The JVP differentiates both the state increment and a moving predictor:

\[
  \dot r_s = \dot t_u^T(u-u_0)+t_u^T(\dot u-\dot u_0)
           +\dot t_p^T(p-p_0)+t_p^T(\dot p-\dot p_0)-\dot{\Delta s}.
\]

The VJP returns the corresponding real Euclidean cotangents. All arrays have
fixed dimensions and finite values. A changed branch, cut topology, or
interface graph is a client-level event and is not differentiated here.

`normalize_pseudo_arclength_tangent` treats the state and continuation
parameter tangent as one Euclidean vector and returns its unit blocks and
norm. Its JVP/VJP products make tangent normalization differentiable while
rejecting a zero or non-finite predictor. Metric-weighted normalization and
tangent construction remain caller policies.

`evaluate_residual_merit` supplies the weighted least-squares scalar
\(\phi(R)=\tfrac12\sum_i w_iR_i^2\) used by line-search and trust-region
clients. Positive caller-owned weights are required; the routine reports the
merit and its analytical JVP/VJP but deliberately does not make a step
acceptance decision.

`assemble_pseudo_transient_residual` supplies the separate continuation
stabilizer
\[
  R_{pt}=R(u)+M(u-u_{old})/\Delta t.
\]
Its matrix, state, previous-state, and positive-step JVP/VJP products are
analytical. This is a solver-neutral pseudo-transient hook and is not a
replacement for the symplectic or dissipative time integrators.

## API

```fortran
call assemble_pseudo_arclength_residual( &
    residual, state, parameter, previous_state, previous_parameter, &
    tangent_state, tangent_parameter, step, augmented_residual, status)
call assemble_pseudo_arclength_residual_jvp( &
    residual, state, parameter, previous_state, previous_parameter, &
    tangent_state, tangent_parameter, step, residual_dot, state_dot, &
    parameter_dot, previous_state_dot, previous_parameter_dot, &
    tangent_state_dot, tangent_parameter_dot, step_dot, augmented_dot, status)
call assemble_pseudo_arclength_residual_vjp( &
    residual, state, parameter, previous_state, previous_parameter, &
    tangent_state, tangent_parameter, step, augmented_bar, residual_bar, &
    state_bar, parameter_bar, previous_state_bar, previous_parameter_bar, &
    tangent_state_bar, tangent_parameter_bar, step_bar, status)

call normalize_pseudo_arclength_tangent( &
    tangent_state, tangent_parameter, normalized_state, normalized_parameter, &
    norm, status)
call normalize_pseudo_arclength_tangent_jvp( &
    tangent_state, tangent_parameter, tangent_state_dot, tangent_parameter_dot, &
    normalized_state_dot, normalized_parameter_dot, norm_dot, status)
call normalize_pseudo_arclength_tangent_vjp( &
    tangent_state, tangent_parameter, normalized_state_bar, normalized_parameter_bar, &
    norm_bar, tangent_state_bar, tangent_parameter_bar, status)

call evaluate_residual_merit(residual, weights, merit, status)
call evaluate_residual_merit_jvp( &
    residual, weights, residual_dot, weights_dot, merit_dot, status)
call evaluate_residual_merit_vjp( &
    residual, weights, merit_bar, residual_bar, weights_bar, status)

call assemble_pseudo_transient_residual( &
    residual, mass, state, previous_state, time_step, augmented_residual, status)
call assemble_pseudo_transient_residual_jvp( &
    residual, mass, state, previous_state, time_step, residual_dot, mass_dot, &
    state_dot, previous_state_dot, time_step_dot, augmented_residual_dot, status)
call assemble_pseudo_transient_residual_vjp( &
    residual, mass, state, previous_state, time_step, augmented_residual_bar, &
    residual_bar, mass_bar, state_bar, previous_state_bar, time_step_bar, status)
```

`test_pseudo_arclength_residual` checks the value against the explicit scalar
formula, the JVP against both the product rule and a central difference, and
the VJP against an independent real dot-product identity. It is intentionally
solver- and application-neutral: free-boundary geometry and force functionals
are supplied by higher layers.
`test_residual_merit` independently checks the weighted scalar formula, a
central-difference JVP, the real dot-product VJP identity, and rejection of a
non-positive weight.
`test_pseudo_transient_residual` checks the matrix/state/time product rule,
central difference, full real adjoint identity, and non-positive-step
rejection against an independent manufactured expression.
