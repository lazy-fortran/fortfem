---
title: Nonlinear resistive-MHD block composition
---

# Nonlinear resistive-MHD block composition

`fortfem_nonlinear_resistive_mhd_composition` is a closure-neutral assembly
boundary for nonlinear resistive-MHD clients.  It does not choose a state
vector, constitutive law, normalization, geometry, time integrator, or plasma
closure.  The client supplies one callback that evaluates each named block:

1. Faraday;
2. Ampère;
3. momentum;
4. scalar pressure;
5. tensor pressure/stress;
6. anisotropic transport;
7. conducting-wall;
8. free-boundary.

Each callback receives the same caller-owned state and returns a residual
vector plus three scalar diagnostics: stored energy, signed input power, and
non-negative dissipation.  FortFEM sums those values and reports

\[
  R = \sum_b R_b,\qquad
  E = \sum_b E_b,\qquad
  P = \sum_b P_b,\qquad
  D = \sum_b D_b,\qquad
  \mathcal P = P-D.
\]

The positive-dissipation check is intentionally the only physical policy in
this layer: a callback that reports a negative dissipation is rejected.  The
reported balance is an instantaneous power balance.  A client may compare it
with its own time derivative of `stored_energy`; FortFEM does not silently
declare a static residual to be a time integrator or claim that a dissipative
step is symplectic.

## Derivative contract

The value, JVP, and real VJP interfaces are exact compositions of the
caller-owned callbacks.  The VJP includes cotangents for all four reported
ledger outputs.  Since `balance = input_power - dissipation`, its cotangent is
propagated with the corresponding signs before the callback is called.  This
allows a client to differentiate a residual, an energy objective, a power
constraint, or a passivity diagnostic without differentiating through a
procedural block loop.

The callbacks can internally call compatible H(curl)/H(div) Faraday/Ampère
operators, tensor-valued pressure, field-aligned transport, FEM/BEM/DtN/PML
ports, or wall/free-boundary maps.  Those operators remain independently
testable and externally owned.  This module only supplies their common
nonlinear composition point.

## Scope boundary

Continuation, branch selection, topology events, and nonlinear solve policy
remain with the caller.  The existing pseudo-arclength, merit, pseudo-
transient, and continuation-event primitives can be composed around this
residual.  No JOREK, MHD equilibrium, GEQDSK reader, species model, or
application-specific file format is included.

The independent test
`test_nonlinear_resistive_mhd_composition` uses a nonlinear manufactured law
for all eight blocks and checks the value against a separate elementwise
oracle, the JVP against both an independent derivative and central
differences, and the VJP against the full residual-plus-ledger dot-product
identity.
