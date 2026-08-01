---
title: Complex cycle-period constraints
---

# Complex cycle-period constraints

`assemble_period_constraints` contracts a fixed real cycle basis with a
complex edge field:

\[
  r_p = \sum_e C_{ep} a_e - t_p.
\]

The cycle columns carry orientation and topology; FortFEM does not assign
flux units, choose a gauge, or identify a physical loop.  This makes the
block usable with tree--cotree reductions, harmonic representatives,
Nédélec/IGA edge fields, and external FEM/BEM or linear-response clients.

The JVP differentiates the cycle basis, complex edge values, and target
periods.  The VJP uses the real part of the complex inner product.  In
particular, the real cycle-basis cotangent is

\[
  \bar C_{ep} = \operatorname{Re}(\overline{\bar r_p}a_e),
\]

while `edge_values_bar = C*residual_bar` and
`target_periods_bar = -residual_bar`.

`test_period_constraints` compares the residual with a direct matrix oracle,
checks a finite-difference JVP, verifies the complex adjoint identity, and
rejects incompatible dimensions.  Topology changes and cycle-basis rebuilds
remain discrete events outside this fixed-topology derivative contract.
