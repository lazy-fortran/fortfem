---
title: Module Overview
---

# FortFEM Module Overview

FortFEM is organized into several key modules, each responsible for specific functionality. This page provides an overview of the main modules and their purposes.

## Core Modules

### Mesh Modules
- **fortfem_mesh_2d**: 2D triangular mesh data structures and operations
- **fortfem_mesh_1d**: 1D interval mesh handling
- **fortfem_mesh_utils**: Utility functions for mesh manipulation

### Function Space Modules  
- **fortfem_function_space**: Base functionality for finite element spaces
- **fortfem_p1_space**: P1 (linear) Lagrange elements
- **fortfem_p2_space**: P2 (quadratic) Lagrange elements
- **fortfem_rt0_space**: Raviart-Thomas H(div) elements
- **fortfem_nedelec_space**: Nédélec H(curl) edge elements

### Basis Function Modules
- **fortfem_basis_p1_1d**: 1D linear basis functions
- **fortfem_basis_p1_2d**: 2D linear basis functions on triangles
- **fortfem_basis_p2_2d**: 2D quadratic basis functions
- **fortfem_basis_rt0**: RT0 vector basis functions
- **fortfem_basis_nedelec**: Nédélec edge basis functions

### Assembly Modules
- **fortfem_assembly_1d**: 1D finite element assembly
- **fortfem_assembly_2d**: 2D finite element assembly
- **fortfem_boundary_assembly_2d**: Boundary integral assembly
- **fortfem_vector_assembly**: Assembly for vector-valued problems

### Solver Modules
- **fortfem_solver**: Main solver interface
- **fortfem_sparse_matrix**: Sparse matrix storage and operations
- **fortfem_gmres**: GMRES iterative solver implementation

### Integration Modules
- **fortfem_gauss_quadrature**: Gaussian quadrature rules
- **fortfem_quadrature_1d**: 1D quadrature specialization
- **fortfem_quadrature_2d**: 2D quadrature on triangles

### Utility Modules
- **fortfem_kinds**: Precision definitions and constants
- **fortfem_utils**: General utility functions
- **fortfem_io**: Input/output routines

## Module Dependencies

The modules follow a layered architecture:

1. **Foundation Layer**: kinds, utils
2. **Mesh Layer**: mesh_1d, mesh_2d, mesh_utils  
3. **Basis Layer**: basis functions for each element type
4. **Space Layer**: function spaces built on meshes and basis functions
5. **Assembly Layer**: matrix and vector assembly using spaces
6. **Solver Layer**: solution of assembled systems

## Using the Modules

Most users will interact with FortFEM through the main module:

```fortran
use fortfem
```

This provides access to all public types and procedures. For more control, individual modules can be imported:

```fortran
use fortfem_mesh_2d
use fortfem_p1_space
use fortfem_assembly_2d
```

## Module Documentation

Detailed documentation for each module, including:
- Public types and their components
- Public procedures and their interfaces
- Usage examples
- Implementation notes

Can be found in the [module reference](../modules.html) section of this documentation.