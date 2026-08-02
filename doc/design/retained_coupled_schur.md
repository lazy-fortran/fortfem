---
title: Retained coupled Schur contract
---

# Retained coupled Schur contract

The assemble_retained_coupled_schur operation eliminates a retained field
block without materializing a monolithic matrix. For

\[
\begin{bmatrix}E&C\\D&F\end{bmatrix}
\begin{bmatrix}y\\x\end{bmatrix}
=
\begin{bmatrix}f\\g\end{bmatrix},
\]

where the lower-right block is represented by retained_field_split_t or
retained_complex_field_split_t, the returned system is

\[
S y = r,\qquad S=E-C D^{-1}F,\qquad r=f-C D^{-1}g.
\]

The exterior blocks E, C, F, and their right-hand sides are caller-owned.
The retained split owns only the square per-field CSC factors. This keeps
FEM/BEM, DtN, PML, wall, HDG, Fourier, and MHD composition policies outside
the numerical primitive.

The value, JVP, and VJP paths are public through fortfem_api. The JVP
propagates derivatives of all exterior blocks, right-hand sides, and retained
CSC values. The VJP returns the corresponding block cotangents and one CSC
cotangent per retained field. Complex VJPs use the real-part inner product,
consistent with the frequency-domain contracts.

The independent test test_retained_coupled_schur compares value and
central-difference JVP results with a diagonal dense elimination oracle and
checks real and complex adjoint identities. The operation is a solver
composition contract, not an equilibrium or MHD model.
