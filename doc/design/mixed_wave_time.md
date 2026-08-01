---
title: Mixed first-order wave time stepping
---

# Mixed first-order wave time stepping

`advance_mixed_wave_midpoint` advances the compatible semidiscrete state

\[
 M_q\dot q + C^T v = 0, \qquad
 M_v\dot v - Cq = 0.
\]

The state may represent pressure/particle velocity, displacement/momentum,
electric/magnetic variables, or another port-Hamiltonian wave pair. The
implicit-midpoint block is

\[
\begin{bmatrix}
 M_q & \frac{\Delta t}{2}C^T\\
 -\frac{\Delta t}{2}C & M_v
\end{bmatrix}
\begin{bmatrix}q_{n+1}\\v_{n+1}\end{bmatrix}
=
\begin{bmatrix}
 M_q & -\frac{\Delta t}{2}C^T\\
 \frac{\Delta t}{2}C & M_v
\end{bmatrix}
\begin{bmatrix}q_n\\v_n\end{bmatrix}.
\]

For symmetric positive mass blocks, this Cayley map preserves

\[
 E=\tfrac12q^TM_qq+\tfrac12v^TM_vv
\]

and reversing the signed step reverses the state to roundoff. The focused
oscillator test checks the independent Cayley formula, energy preservation,
and reversibility. Dissipative terms are intentionally not hidden in this
routine; they belong in a declared split or metriplectic substep.

The same contract is suitable for tensor-valued pressure, mixed elasticity,
Maxwell, and acoustic-elastic coupling once their compatible mass and
interconnection blocks are assembled.

## Partitioned symplectic Euler

`advance_mixed_wave_symplectic_euler` is the explicit partitioned companion
for the ideal (nondissipative) part.  It first solves

\[
 M_v v_{n+1}=M_vv_n+\Delta t\,Cq_n,
\]

and then solves

\[
 M_q q_{n+1}=M_qq_n-\Delta t\,C^Tv_{n+1}.
\]

This is a first-order symplectic map for the canonical mixed pair (with the
mass blocks defining the corresponding pairing).  Unlike implicit midpoint,
it generally oscillates around the exact quadratic energy rather than
preserving it exactly.  The routine validates all block and state dimensions,
solves both mass systems before committing either part of the state, and
returns a singular status without partially updating the state.

The independent test checks the partitioned oscillator formula and the
canonical two-state symplectic-form identity.  Dissipation, conductivity,
viscosity, and PML terms remain separate declared substeps; they must not be
silently folded into this ideal symplectic update.

The same step exposes analytical derivative products.  Given increments
`(q_dot, v_dot, time_step_dot)`,
`advance_mixed_wave_symplectic_euler_jvp` differentiates the two mass solves
in sequence:

\[
 M_v\dot v^+ = M_v\dot v+hC\dot q+\dot h Cq,
 \qquad
 M_q\dot q^+ = M_q\dot q-hC^T\dot v^+-\dot h C^Tv^+.
\]

`advance_mixed_wave_symplectic_euler_vjp` applies the corresponding real
transpose solves and returns cotangents for `q`, `v`, and `h`.  The focused
test checks the JVP against a central difference of the complete update and
checks the VJP dot-product identity.  Mass and coupling blocks are held fixed
by this contract; a caller differentiating geometry or constitutive data
should differentiate those blocks in its surrounding residual.
