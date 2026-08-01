---
title: Symmetric mixed-wave splitting
---

# Symmetric mixed-wave splitting

`advance_mixed_wave_strang` composes two caller-owned ideal mixed-wave
operators with a symmetric Cayley split,

\[
  \mathcal A(\Delta t/2)\,\mathcal B(\Delta t)\,\mathcal A(\Delta t/2).
\]

Each factor is the implicit-midpoint map documented in
[mixed first-order wave time stepping](mixed_wave_time.html).  The split is
therefore reversible for signed time steps and preserves the quadratic
mass-weighted energy of each factor.  The composition is useful for pressure-
velocity acoustics, displacement-momentum elasticity, Maxwell blocks, and
other first-order port-Hamiltonian wave states.  Dissipative, resistive,
viscous, or PML factors remain explicit caller-owned substeps and must not be
silently inserted into this ideal map.

The API also exposes `advance_mixed_wave_strang_jvp` and
`advance_mixed_wave_strang_vjp`.  The JVP propagates state and signed-step
increments through all three factors.  The VJP reverses the factors and
returns the state and time-step cotangents.  The independent test checks
energy preservation, signed-step reversal, a central-difference tangent, and
the real dot-product adjoint identity on nonidentity mass and noncommuting
coupling blocks.

The block matrices and their geometry/constitutive derivatives remain
caller-owned.  This keeps the time integrator composable with FEM, IGA,
mixed, DG, and structure-preserving spatial discretizations without hiding a
plasma model or a particular constitutive closure in the time-step API.
