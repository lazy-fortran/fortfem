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
```

`test_pseudo_arclength_residual` checks the value against the explicit scalar
formula, the JVP against both the product rule and a central difference, and
the VJP against an independent real dot-product identity. It is intentionally
solver- and application-neutral: free-boundary geometry and force functionals
are supplied by higher layers.
