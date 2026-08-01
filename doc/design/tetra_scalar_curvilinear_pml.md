---
title: Tetrahedral scalar curvilinear PML assembly
---

# Tetrahedral scalar curvilinear PML assembly

`fortfem_assembly_tetra_lagrange_curvilinear_pml_3d` provides the P1
scalar-Helmholtz volume path for a full complex curvilinear stretch.  For each
tetrahedron it consumes the scalar tensors

\[
G=\det(S)S^{-1}S^{-T},\qquad M=\det(S),
\]

and assembles

\[
K_{rc}=|T|\left(\nabla\phi_r^T G\nabla\phi_c
                  - k^2 M\,m_{rc}\right).
\]

The local and global public contracts are:

```fortran
call assemble_tetra_lagrange_curvilinear_pml_element( &
    vertices, degree, quadrature_degree, stretch, wave_number, matrix, status)
call assemble_tetra_lagrange_curvilinear_pml_element_jvp( &
    vertices, degree, quadrature_degree, stretch, wave_number, vertices_dot, &
    stretch_dot, wave_number_dot, matrix_dot, status)
call assemble_tetra_lagrange_curvilinear_pml_element_vjp( &
    vertices, degree, quadrature_degree, stretch, wave_number, matrix_bar, &
    vertices_bar, stretch_bar, wave_number_bar, status)

call assemble_tetra_lagrange_curvilinear_pml_csc( &
    mesh_vertices, tetrahedra, degree, stretch, wave_number, matrix, status)
call assemble_tetra_lagrange_curvilinear_pml_csc_jvp( &
    mesh_vertices, tetrahedra, degree, stretch, wave_number, &
    mesh_vertices_dot, stretch_dot, wave_number_dot, matrix_dot, status)
call assemble_tetra_lagrange_curvilinear_pml_csc_vjp( &
    mesh_vertices, tetrahedra, degree, stretch, wave_number, matrix_values_bar, &
    mesh_vertices_bar, stretch_bar, wave_number_bar, status)
```

The current contract intentionally follows the existing executable scalar
PML path and supports P1 tetrahedra.  It preserves the merged CSC pattern and
scatters shared-vertex and per-cell stretch cotangents.  The local derivatives
compose the FortSym-backed curvilinear coefficient products with analytical
tetrahedral determinant and inverse maps.  The independent
`test_tetra_lagrange_curvilinear_pml_ad` test checks the closed-form P1 oracle,
central differences, local and global complex adjoint identities, and CSC
pattern preservation.

The constrained solved-state wrapper is
`solve_tetra_lagrange_curvilinear_pml`, with `_jvp` and `_vjp` variants.  It
uses the same fixed CSC pattern and the complex implicit direct-solve contract
as the vector PML path.  The state test independently re-solves perturbed
systems for the JVP and checks the real-part complex adjoint identity for
volume load, constrained values, geometry, stretch, and wave number.
