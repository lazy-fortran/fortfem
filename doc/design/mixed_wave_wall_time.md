---
title: Structure-preserving mixed wave--wall time step
---

# Structure-preserving mixed wave--wall time step

`advance_mixed_wave_wall_midpoint` composes a first-order mixed wave port with
a neutral resistive-wall RL block:

\[
 M_q\dot q+C^Tv=0,\qquad
 M_v\dot v-Cq-P^Ti=0,\qquad
 L\dot i+Ri+Pv=0.
\]

The implicit-midpoint block is assembled as one coupled solve. The `P` terms
are skew in the field/wall power pairing, so the independent ledger
`evaluate_mixed_wave_wall_energy_balance` verifies

\[
 H_{n+1}-H_n+h,i_{n+1/2}^TRi_{n+1/2}=0.
\]

The value, full JVP, and real VJP actions cover matrices, port coupling,
initial state, and signed step data. Geometry, wall basis, boundary
orientation, and application normalization remain caller-owned. This is the
composition layer needed by FEM/BEM/DtN field ports and STARWALL-like clients;
the existing ideal mixed-wave symplectic and separate dissipative steps remain
available when a split integrator is preferred.
