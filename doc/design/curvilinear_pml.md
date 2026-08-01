---
title: Curvilinear PML coefficient contract
---

# Curvilinear PML coefficient contract

`fortfem_curvilinear_helmholtz_pml` is the local coefficient boundary for a
complex, nonsingular three-dimensional coordinate stretch. It is neutral: the
caller supplies the stretch produced by a mesh, an IGA map, or a curved-layer
construction; FortFEM does not choose a PML geometry or an application wave
number.

For a complex stretch matrix \(S\), the scalar Helmholtz weak coefficients are

\[
 G = \det(S)S^{-1}S^{-T},\qquad M=\det(S).
\]

The curl--curl Maxwell coefficients use the covariant pullback

\[
 C=\det(S)^{-1}S^TS,\qquad
 M_c=\det(S)S^{-1}S^{-T}.
\]

The transpose is deliberately not a Hermitian transpose. A coordinate stretch
is differentiated holomorphically; complex conjugation enters only the
real-vector reverse product. The public primal, JVP, and VJP routines are:

```fortran
call curvilinear_scalar_helmholtz_pml_coefficients( &
    stretch, gradient_coefficient, mass_coefficient, status)
call curvilinear_scalar_helmholtz_pml_coefficients_jvp( &
    stretch, stretch_dot, gradient_dot, mass_dot, status)
call curvilinear_scalar_helmholtz_pml_coefficients_vjp( &
    stretch, gradient_bar, mass_bar, stretch_bar, status)

call curvilinear_curl_curl_pml_coefficients( &
    stretch, curl_coefficient, mass_tensor, status)
call curvilinear_curl_curl_pml_coefficients_jvp( &
    stretch, stretch_dot, curl_dot, mass_tensor_dot, status)
call curvilinear_curl_curl_pml_coefficients_vjp( &
    stretch, curl_bar, mass_tensor_bar, stretch_bar, status)
```

The reverse routines use

\[
 \operatorname{Re}\sum_i\overline{y_{\!bar,i}}\,\dot y_i
 =
 \operatorname{Re}\sum_i\overline{S_{\!bar,i}}\,\dot S_i,
\]

so they are suitable for real and imaginary geometry-design parameters. A
nonfinite or numerically singular stretch is rejected with a nonzero status;
the output arrays are initialized to zero on failure. The singularity guard is
scale-aware and uses the determinant relative to the largest matrix entry.

The diagonal case is exactly the existing Cartesian tensor path, while the
off-diagonal entries retain shear and curvilinear metric coupling. This keeps
IGA and geometry-generated layers from being forced into a diagonal
approximation. `test_curvilinear_helmholtz_pml` independently checks the
closed-form values, two-sided directional differences, the real complex
adjoint identities, diagonal reduction, and singular-input rejection.

The same tensors are now consumed by the tetrahedral first-kind Nédélec local
form. The public
`assemble_tetra_nedelec_curvilinear_pml_element`, `_jvp`, and `_vjp` routines
differentiate the physical covariant Piola map, determinant quadrature,
wave-number term, and every stretch entry together. A diagonal stretch is
checked against the established Cartesian element assembly, while a full
shear stretch is checked by reassembly finite differences and the real complex
adjoint identity in
`test_tetra_nedelec_curvilinear_pml_element_ad`.

For a mesh with one full stretch per tetrahedron, the corresponding CSC
assembly routines are `assemble_tetra_nedelec_curvilinear_pml_csc`, `_jvp`,
and `_vjp`. They preserve the existing merged sparsity pattern and accumulate
shared-vertex and per-cell stretch cotangents. The global behavioral oracle
`test_tetra_nedelec_curvilinear_pml_csc_ad` checks pattern preservation,
complete reassembly differences, and the shared-geometry complex adjoint.

The solved-state contract is exposed as
`solve_tetra_nedelec_curvilinear_pml`, `_jvp`, and `_vjp`. It applies the same
fixed-topology constrained direct solve and optional boundary-form coupling as
the Cartesian PML state path, but retains the full per-cell stretch tensor.
`test_tetra_nedelec_curvilinear_pml_state_ad` checks the manufactured solve,
independent re-solves for the state JVP, and the real-complex state adjoint.

The module is a coefficient contract, not a complete curved-object PML
assembler. Layer geometry, active-cell topology, quadrature, and reflection
diagnostics remain caller-owned. The solver-neutral weighted reflection and
error metrics used to compare those samples are documented in
`wave_reflection_diagnostics.md`; geometry and layer-topology chains remain
caller-owned and are tracked in the open-boundary roadmap.
