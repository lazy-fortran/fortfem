---
title: Conservative FCI power and flux ledger
---

# Conservative FCI power and flux ledger

`evaluate_fci_power_flux_ledger` is the neutral accounting boundary between an
FCI parallel support operator and a caller-owned perpendicular action.  For a
staggered support gradient (g=Qu), coefficient (k_parallel), and positive
support measure (W_s), it reports

\[
  f_\parallel=-k_\parallel g,\qquad
  P_\parallel=-g^T W_s k_\parallel g.
\]

The caller supplies the canonical-cell perpendicular action (a_\perp) and
positive canonical measure (W_c).  Its signed power is

\[
  P_\perp=u^T W_c a_\perp,
  \qquad P=P_\parallel+P_\perp.
\]

Thus a dissipative perpendicular action is negative in the same convention as
the existing support-operator energy identity.  The ledger does not assemble
the FCI map, choose a material trace, impose a boundary condition, or attach a
species/sheath interpretation.  The same contract can therefore consume FEM,
IGA, Fourier, DG, or plane-local actions.

FortSym emits the pointwise value, fixed-topology JVP, and VJP products in
`src/generated/fortfem_fci_power_flux_products.f90`; the wrapper only traverses
caller-owned arrays and accumulates the scalar ledger.  The focused test uses a
separate direct-sum oracle, central differences, and a real transpose identity.

