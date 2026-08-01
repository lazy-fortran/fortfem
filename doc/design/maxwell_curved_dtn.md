---
title: Curved Maxwell trace-to-flux maps
---

# Curved Maxwell trace-to-flux maps

`fortfem_maxwell_curved_dtn` builds a finite-dimensional Maxwell
Dirichlet-to-Neumann map from weak boundary operators. The contract is useful
for a curved exterior surface when a pointwise vector spherical or planar
symbol is unavailable.

## Mixed trace convention

Let `j` contain coefficients of an equivalent tangential surface current. The
three input matrices have the following meaning:

\[
 Z_E j = \langle w_E, E_t(j)\rangle,
 \qquad
 Z_H j = \langle w_H, H_t(j)\rangle,
 \qquad
 M_E e = \langle w_E, e\rangle .
\]

`Z_E` and `Z_H` may use different dual test spaces. The trace coefficient
vector `e` is represented in the primal space associated with `M_E`. If
`Z_E` is nonsingular, the weak flux map is

\[
 D = Z_H Z_E^{-1} M_E,
 \qquad
 H_t(D e) = D e.
\]

The result is a weak trace map. Its rows are flux dual coefficients, so the
map must be paired with the corresponding dual basis when it is inserted into
a FEM or IGA weak form. It is not a pointwise interpolation of a tangential
vector field.

The public routines are:

```text
assemble_maxwell_trace_to_flux_map(Z_E, Z_H, M_E, D, status)
apply_maxwell_trace_to_flux(Z_E, Z_H, M_E, e, h, status)
apply_maxwell_trace_to_flux_map(D, e, h, status)
apply_maxwell_trace_to_flux_jvp(...)
apply_maxwell_trace_to_flux_vjp(...)
assemble_maxwell_weak_trace_reconstruction(A, B, R, status)
apply_maxwell_weak_trace_reconstruction(A, B, h, values, status)
apply_maxwell_weak_trace_reconstruction_jvp(...)
apply_maxwell_weak_trace_reconstruction_vjp(...)
```

The matrix-free action solves `Z_E j = M_E e` and then evaluates `Z_H j`.
The JVP solves the tangent equation

\[
 Z_E \dot j = \dot M_E e + M_E \dot e - \dot Z_E j.
\]

The VJP uses the conjugate-transpose solve for the adjoint current. All
complex reverse products use
`real(sum(conjg(output_bar)*output_dot))`.

## Pointwise reconstruction

The weak map deliberately does not pretend that RBC coefficients are point
values. For a caller-owned primal evaluation matrix `B` and weak mass matrix
`A`, where

\[
 A c = h, \qquad u_{\rm points}=B c,
\]

the reconstruction layer returns `R = B A^{-1}` and applies it to a weak
flux vector. `B` may evaluate an RWG, scalar, cut, or IGA trace basis at any
requested points; FortFEM does not choose those points or the basis. The
reconstruction JVP solves the tangent mass system

\[
 A\dot c=\dot h-\dot A c,
 \qquad
 \dot u=\dot Bc+B\dot c,
\]

and the VJP uses the conjugate-transpose mass solve. This makes a physical
field/slice plot or a client FEM coupling explicit without changing the weak
DtN convention.

## Exact-curved torus composition

`assemble_maxwell_torus_curved_dtn_rwg_3d` assembles:

* the exact-parametric torus RWG EFIE as `Z_E`;
* the one-sided extrapolated torus MFIE in the RBC dual space as `Z_H`;
* the torus RWG mass matrix as `M_E`.

The wrapper returns the corresponding RWG-to-RBC weak map. Coincident and
adjacent panel quadrature, the surface Piola geometry, and the RBC dual trace
remain owned by the existing torus BEM operators. The new layer only composes
their compatible matrices and exposes the same trace convention to FEM, IGA,
and larger-domain clients.

This construction follows the Calderon trace framework of Scroggs et al.
([arXiv:1703.10900](https://arxiv.org/abs/1703.10900)) and the dual-complex
construction of Buffa and Christiansen
([doi:10.1090/S0025-5718-07-01965-5](https://doi.org/10.1090/S0025-5718-07-01965-5)).
The torus MFIE limit and EFIE use the conventions documented in the exact
curved torus Maxwell example.

## Verification boundary

`test_maxwell_curved_dtn` checks the matrix-free action against an assembled
map, the complete JVP against a central difference, and the VJP against the
complex adjoint identity on independent nonsymmetric complex matrices. It
also assembles and applies the exact-curved torus map. The algebraic test is
independent of the torus quadrature, while the torus call verifies that the
existing RWG, RBC, EFIE, MFIE, and mass operators have compatible dimensions.

The map requires a fixed surface topology and a nonsingular electric form.
Resonance handling, regularization, and topology changes remain explicit
client responsibilities. The pointwise reconstruction layer also requires a
fixed basis topology and nonsingular weak mass; it does not regularize a
resonant BEM form or differentiate topology changes.
