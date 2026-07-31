---
title: Design Documentation
---

# FortFEM Design Documentation

## Overview

FortFEM is designed with the following principles:
- **Natural mathematical notation**: Express weak forms as you would write them mathematically
- **Modularity**: Clear separation between mesh, function spaces, assembly, and solvers
- **Performance**: Efficient sparse matrix operations and assembly routines
- **Extensibility**: Easy to add new element types and problem formulations

The [differentiation strategy](differentiation.html) defines the common
primal/JVP/VJP API, analytical FortSym path, optional Enzyme tournament, and
implicit sparse-solve adjoints.

The [FCI parallel support-operator contract](fci_parallel_operator.html)
defines the PARALLAX-aligned mapped gradient and its conservative
volume-weighted adjoint without coupling FortFEM to a particular field-line
tracer.

The [mixed first-order wave step](mixed_wave_time.html) documents the common
pressure/velocity, displacement/momentum, and electromagnetic midpoint/Cayley
contract with independent energy and reversibility checks.

The [CGL pressure tensor](cgl_pressure_tensor.html) documents the generated
gyrotropic tensor constitutive block and its JVP/VJP contract.

## Architecture

### Core Components

1. **Mesh Module** (`fortfem_mesh_*`)
   - Handles mesh data structures and connectivity
   - Supports triangular and quadrilateral elements
   - Provides mesh I/O and generation utilities

2. **Function Space Module** (`fortfem_function_space`)
   - Defines finite element spaces on meshes
   - Manages degrees of freedom (DOF) numbering
   - Handles boundary condition identification

3. **Basis Functions** (`fortfem_basis_*`)
   - Implements shape functions and their derivatives
   - Supports P1/P2/Q1 and verified lowest-order Nédélec/RT0 kernels
   - Handles reference-to-physical element transformations

4. **Assembly Module** (`fortfem_assembly_*`)
   - Assembles global matrices and vectors from element contributions
   - Uses efficient quadrature rules
   - Supports various bilinear and linear forms

5. **Solver Module** (`fortfem_solver`)
   - Interfaces to LAPACK and includes custom GMRES implementation
   - Provides both direct and iterative methods
   - Includes GMRES implementation

## Design Patterns

### Object-Oriented Approach
FortFEM uses Fortran's object-oriented features:
```fortran
type :: mesh_2d_t
    real(dp), allocatable :: vertices(:,:)
    integer, allocatable :: triangles(:,:)
    integer :: nvert, ntri
contains
    procedure :: init => mesh_2d_init
    procedure :: destroy => mesh_2d_destroy
end type
```

### Generic Programming
The library uses interfaces for flexibility:
```fortran
interface assemble
    module procedure assemble_mass_p1
    module procedure assemble_stiffness_p1
    module procedure assemble_mass_p2
end interface
```

## Element Types

### Lagrange Elements
- **P1**: Linear polynomials on triangles (implemented)
- **P2**: Quadratic polynomials with 6 DOFs per triangle (implemented)

### Vector Elements

- **Nédélec**: First- and second-kind triangular families through order four,
  with covariant Piola maps and oriented sparse curl-mass forms.
- **Tetrahedral Nédélec**: First-kind reference bases through order five have
  verified edge, face, and cell moments and curls. The affine covariant map
  and local curl-mass assembly accept all implemented orders. Canonical edge
  signs and generated/runtime face transforms provide orientation-aware
  global topology and direct `fortsparse` CSC assembly through order five.
- **Raviart--Thomas and BDM**: Triangular families through order four, with
  contravariant Piola maps and oriented sparse divergence-mass forms.
- Weighted vector forms compile to local matrices and `fortsparse` CSC
  matrices. The order-one Nédélec solve supports cellwise scalar curl and mass
  coefficients, constant or cellwise vector sources, constant physical
  tangential data, owned nonconstant edge moments, and a sparse direct
  interior solve. Its three-material Fourier-mode magnetic convergence has an
  analytical oracle, and a smooth manufactured curl-mass problem attains
  first-order field convergence. General pointwise sources, tensor
  coefficients, and higher-order high-level solves remain experimental.

## Weak Form Framework

The design allows natural expression of weak forms:
```fortran
! Poisson equation: -∆u = f
! Weak form: ∫∇u·∇v dx = ∫fv dx
call assemble_stiffness(mesh, V, A)
call assemble_load(mesh, V, f, b)
```

## Memory Management

- Uses allocatable arrays for dynamic memory
- Automatic cleanup through finalizers
- Efficient sparse storage formats (CSR/CSC)

## Future Directions

- 3D element support
- Adaptive mesh refinement
- Parallel assembly and solvers
- More element types (DG, HDG)
