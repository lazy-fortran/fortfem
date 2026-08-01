---
title: Maxwell FEM--BEM state differentiation
---

# Maxwell FEM--BEM state differentiation

The mixed volume/trace system assembled by FortFEM has the neutral form

\[
A(\mathbf p)x(\mathbf p)=b(\mathbf p),
\]

where `x` may concatenate volume Nédélec coefficients and RWG/BC surface
unknowns. `solve_maxwell_fem_bem_linear_state` is the value map. Its tangent
contract is the implicit relation

\[
A\,\dot x=\dot b-\dot A\,x,
\]

implemented by `solve_maxwell_fem_bem_linear_state_jvp`. The corresponding
real-complex adjoint is exposed by
`solve_maxwell_fem_bem_linear_state_vjp`; it satisfies

\[
\operatorname{Re}(\bar x^H\dot x)=
\operatorname{Re}(\bar A^H\dot A+\bar b^H\dot b).
\]

The layer deliberately does not differentiate mesh topology, quadrature
branching, or a particular FEM/BEM assembly. A caller composes its
FortSym-generated matrix and right-hand-side products here. The toroidal
FEM--BEM value solver uses the same value contract, so a derivative-aware
client cannot accidentally use a different solve convention.

`test_maxwell_fem_bem_state_ad` checks the independent central-difference and
real-complex adjoint identities on a nontrivial complex block. Geometry and
coefficient derivatives of the toroidal assembly remain separate roadmap
work; this is the reusable implicit-solve foundation for them. The companion
`test_maxwell_torus_fem_bem_nonzero_state` drives the actual curved torus
volume/trace matrix with a nonzero manufactured RWG current and verifies both
parts of the recovered state.
