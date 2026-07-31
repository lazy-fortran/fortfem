---
title: Dissipative Cayley step
---

# Dissipative Cayley step

`advance_dissipative_cayley` is the structure-aware dissipative counterpart
to FortFEM's Hamiltonian midpoint and symplectic-Euler steps. For a mass
matrix `M`, damping/resistivity matrix `D`, and nonnegative step `h`, it
returns

\[
  x_{n+1}=(M+\tfrac h2D)^{-1}(M-\tfrac h2D)x_n.
\]

When `M` is symmetric positive definite and `D` is symmetric positive
semidefinite, the discrete `M`-energy cannot increase. This is a contract
assumption on caller-supplied blocks, not a hidden projection or a claim that
the map is symplectic. Ideal wave, resistive, viscous, thermal, and absorbing
clients can compose this map in a symmetric split.

```fortran
call advance_dissipative_cayley( &
    mass, damping, time_step, state, state_next, status)
call advance_dissipative_cayley_jvp( &
    mass, damping, time_step, state, mass_dot, damping_dot, time_step_dot, &
    state_dot, state_next_dot, status)
call advance_dissipative_cayley_vjp( &
    mass, damping, time_step, state, state_next, state_next_bar, mass_bar, &
    damping_bar, time_step_bar, state_bar, status)
```

The JVP differentiates both Cayley blocks, the step size, and the incoming
state. The VJP solves the transpose Cayley system and returns mass, damping,
step-size, and state cotangents. Negative steps are rejected so a dissipative
block cannot be mislabeled as a reversible integrator.

`test_dissipative_cayley` checks a direct two-by-two inverse oracle, discrete
energy decay, central-difference JVP, the real adjoint identity, and invalid
state/step rejection.
