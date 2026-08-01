---
title: Manufactured toroidal Maxwell FEM--BEM solution
---

# Manufactured toroidal Maxwell FEM--BEM solution

`maxwell_torus_fem_bem_solution` exercises the neutral volume/surface block
on a small exact solid torus. For a constant physical field (E_0), the
lowest-order covariant Nédélec coefficients are the oriented edge integrals

\[
e_i = E_0\mathbin{\cdot}(x_{i,2}-x_{i,1}).
\]

The coupled system is assembled as

\[
\begin{bmatrix} A & C^T \\ C & -Z_E \end{bmatrix}
\begin{bmatrix} e \\ j \end{bmatrix}
= b,
\]

and the manufactured right-hand side is generated as (b=[A\,C^T;C\,-Z_E]
[e_0;0]). The solve therefore has an independent solved-state oracle:
the volume field must recover (e_0), the RWG surface current must vanish,
and the assembled matrix must remain complex symmetric.

The fixture also reconstructs the solved physical field with the covariant
Piola map. Its first plot is a cross-sectional magnitude, followed by a
vector-arrow plot, an absolute-error map, a one-dimensional cross-section,
and a curved three-dimensional torus view. It is deliberately neutral: no
plasma closure, coil model, or external-equilibrium reader is involved.
