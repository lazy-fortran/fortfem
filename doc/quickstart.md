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

### Current FortFEM API

The library provides a clean, FEniCS-inspired API for finite element problems:

```fortran
program poisson_forms
    use fortfem_api
    implicit none
    
    type(mesh_t) :: mesh
    type(function_space_t) :: Vh
    type(trial_function_t) :: u
    type(test_function_t) :: v
    type(function_t) :: uh, f
    type(dirichlet_bc_t) :: bc
    type(form_expr_t) :: a, L
    
    ! Define problem with minimal code
    mesh = unit_square_mesh(20)
    Vh = function_space(mesh, "Lagrange", 1)
    
    u = trial_function(Vh)
    v = test_function(Vh)
    f = constant(1.0_dp)
    
    ! Weak form: ∫∇u·∇v dx = ∫fv dx
    a = inner(grad(u), grad(v))*dx
    L = f*v*dx
    
    ! Solve and plot
    uh = function(Vh)
    bc = dirichlet_bc(Vh, 0.0_dp)
    call solve(a == L, uh, bc)
    call plot(uh, title="Poisson Solution")
    
end program
```

## Key Concepts

### 1. Meshes
Create meshes with simple factory functions:
```fortran
mesh = unit_square_mesh(20)          ! 20×20 uniform grid on [0,1]²
mesh = rectangle_mesh(20, 30, [0.0, 2.0, 0.0, 3.0])  ! Rectangle mesh
mesh = unit_disk_mesh(5)             ! Unit disk with resolution 5
```

### 2. Function Spaces  
FortFEM supports multiple element types with clean API:
- P1 Lagrange elements (linear, continuous)
- Nédélec edge elements (H(curl) conforming for electromagnetic problems)

```fortran
Vh = function_space(mesh, "Lagrange", 1)        ! P1 scalar elements
Vh = vector_function_space(mesh, "Nedelec", 1)  ! Edge elements for vectors
```

### 3. Forms and Assembly
Natural mathematical notation for weak forms:
```fortran
! Scalar problems
a = inner(grad(u), grad(v))*dx    ! Laplacian bilinear form
L = f*v*dx                        ! Source linear form

! Vector problems  
a = inner(curl(E), curl(F))*dx + inner(E, F)*dx  ! Curl-curl + mass
L = inner(J, F)*dx                                 ! Current source
```

### 4. Solvers and Visualization
Automatic solver selection and one-line plotting:
```fortran
call solve(a == L, uh, bc)                        ! Auto solver selection
call plot(uh, title="Solution", colormap="plasma") ! Instant visualization
call plot(Eh, plot_type="streamplot")             ! Vector field plots
```

## Next Steps

- Explore the [examples](examples/index.html) for more complex problems
- Read the [design documentation](design/index.html) to understand the architecture
- Check the [API reference](../modules.html) for detailed function documentation