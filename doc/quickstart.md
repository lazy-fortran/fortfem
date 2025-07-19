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

### Current Simple API
```fortran
program poisson_simple
    use fortfem
    implicit none
    
    type(mesh_2d_t) :: mesh
    type(sparse_matrix_t) :: A
    real(dp), allocatable :: u(:), f(:)
    integer :: info
    
    ! Create mesh
    call create_unit_square_mesh(mesh, n=20)
    
    ! Assemble system -∆u = f
    call assemble_poisson_2d(mesh, A, f)
    
    ! Apply boundary conditions u = 0 on boundary
    call apply_zero_bc(mesh, A, f)
    
    ! Solve
    allocate(u(size(f)))
    u = f
    call solve_lapack_dense(A, u, info)
    
    ! Save solution
    call write_vtk("solution.vtk", mesh, u)
    
end program
```

### Future FEniCS-style API (in development)
```fortran
program poisson_forms
    use fortfem
    implicit none
    
    type(mesh_t) :: mesh
    type(function_space_t) :: Vh
    type(trial_function_t) :: u
    type(test_function_t) :: v
    type(function_t) :: uh, f
    type(dirichlet_bc_t) :: bc
    
    ! Define problem
    mesh = unit_square_mesh(32, 32)
    Vh = function_space(mesh, "Lagrange", 1)
    
    u = trial_function(Vh)
    v = test_function(Vh)
    f = constant(1.0_dp)
    
    ! Weak form: ∫∇u·∇v dx = ∫fv dx
    a = inner(grad(u), grad(v))*dx
    L = f*v*dx
    
    ! Solve with BC: u = 0 on boundary
    bc = dirichlet_bc(Vh, 0.0_dp, boundary)
    uh = function(Vh)
    call solve(a == L, uh, bc)
    
end program
```

## Key Concepts

### 1. Meshes
Create meshes with simple factory functions:
```fortran
mesh = unit_square_mesh(n=20)        ! n×n uniform grid on [0,1]²
mesh = rectangle_mesh([0,0], [2,1], nx=40, ny=20)  ! Custom rectangle
```

### 2. Function Spaces  
Currently available through simplified API:
- P1 Lagrange elements (linear, continuous)
- P2 Lagrange elements (quadratic, 6 DOFs per triangle)
- Nédélec edge elements (H(curl) conforming)
- Raviart-Thomas elements (H(div) conforming)

Future clean API:
```fortran
Vh = function_space(mesh, "Lagrange", 1)  ! P1 elements
Wh = function_space(mesh, "Lagrange", 2)  ! P2 elements
Eh = function_space(mesh, "Nedelec", 1)   ! Edge elements
```

### 3. Assembly
High-level assembly hides matrix details:
```fortran
call assemble_poisson_2d(mesh, A, f)     ! Current: -∆u = f
call apply_zero_bc(mesh, A, f)           ! u = 0 on boundary
```

Future forms-based assembly:
```fortran
a = inner(grad(u), grad(v))*dx    ! Bilinear form
L = f*v*dx                        ! Linear form
```

### 4. Solvers
Simple solver interface:
```fortran
call solve_lapack_dense(A, u, info)      ! Current: dense solver
call solve(a == L, uh, bc)               ! Future: automatic solver
```

## Next Steps

- Explore the [examples](examples/index.html) for more complex problems
- Read the [design documentation](design/index.html) to understand the architecture
- Check the [API reference](../modules.html) for detailed function documentation