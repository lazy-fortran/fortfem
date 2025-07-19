---
title: Quick Start Guide
---

# Quick Start Guide

This guide will help you get started with FortFEM quickly.

## Installation

### Prerequisites

- Modern Fortran compiler (gfortran 9+ or ifort)
- LAPACK and BLAS libraries
- [Fortran Package Manager (fpm)](https://github.com/fortran-lang/fpm)
- Optional: SuiteSparse for advanced sparse solvers

### Building FortFEM

```bash
# Clone the repository
git clone https://github.com/itpplasma/fortfem.git
cd fortfem

# Build the library
fpm build

# Run tests
fpm test

# Build and run examples
fpm run --example
```

## Your First FortFEM Program

Here's a simple example solving the Poisson equation on a unit square:

```fortran
program poisson_example
    use fortfem
    implicit none
    
    type(mesh_2d_t) :: mesh
    type(p1_space_t) :: V
    type(sparse_matrix_t) :: A
    real(dp), allocatable :: b(:), u(:)
    
    ! Create a mesh
    call create_rectangular_mesh(mesh, nx=10, ny=10, &
                                 x_min=0.0_dp, x_max=1.0_dp, &
                                 y_min=0.0_dp, y_max=1.0_dp)
    
    ! Create function space
    call V%init(mesh)
    
    ! Assemble system (weak form: ∫∇u·∇v dx = ∫fv dx)
    call assemble_poisson_2d(mesh, V, A, b)
    
    ! Apply boundary conditions (u = 0 on boundary)
    call apply_dirichlet_bc(A, b, V%boundary_dofs)
    
    ! Solve
    allocate(u(V%ndof))
    call solve_sparse(A, b, u)
    
    ! Output solution
    call write_solution_vtk('solution.vtk', mesh, u)
    
end program
```

## Key Concepts

### 1. Meshes
FortFEM supports triangular and quadrilateral meshes in 1D and 2D:
- `mesh_1d_t`: 1D interval meshes
- `mesh_2d_t`: 2D triangular meshes

### 2. Function Spaces
Available finite element spaces:
- `p1_space_t`: Linear Lagrange elements
- `p2_space_t`: Quadratic Lagrange elements  
- `rt0_space_t`: Raviart-Thomas elements
- `nedelec_space_t`: Nédélec edge elements

### 3. Assembly
The library provides routines for assembling common operators:
- Mass matrix: `∫uv dx`
- Stiffness matrix: `∫∇u·∇v dx`
- Load vector: `∫fv dx`

### 4. Solvers
Multiple solver options:
- Direct: LAPACK for dense, UMFPACK for sparse
- Iterative: GMRES with preconditioning

## Next Steps

- Explore the [examples](examples/index.html) for more complex problems
- Read the [design documentation](design/index.html) to understand the architecture
- Check the [API reference](../modules.html) for detailed function documentation