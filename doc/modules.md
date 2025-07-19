---
title: Module Overview
---

# FortFEM API Overview

FortFEM provides a clean, high-level API that hides implementation details while maintaining the power and flexibility of finite element methods.

## Current API Structure

### High-Level API Module
- **fortfem_api**: Clean FEniCS-style types and constructors
- **fortfem_simple_solvers**: Simplified solver functions

### User-Facing Types (with `_t` suffix)
- **mesh_t**: Mesh wrapper with clean constructors
- **function_space_t**: Finite element function spaces  
- **function_t**: Functions that hold values
- **trial_function_t**: Trial functions for weak forms
- **test_function_t**: Test functions for weak forms
- **form_t**: Base type for weak forms
- **bilinear_form_t**: Forms like `∫∇u·∇v dx`
- **linear_form_t**: Forms like `∫fv dx`
- **dirichlet_bc_t**: Dirichlet boundary conditions

### Factory Functions (lowercase_underscore)
- **unit_square_mesh()**: Create uniform unit square mesh
- **rectangle_mesh()**: Create rectangular mesh with custom bounds
- **function_space()**: Create function space on mesh
- **trial_function()**: Create trial function from space
- **test_function()**: Create test function from space
- **dirichlet_bc()**: Create boundary condition

### High-Level Solvers
- **create_unit_square_mesh()**: Simplified mesh creation
- **assemble_poisson_2d()**: Assemble Poisson system
- **apply_zero_bc()**: Apply homogeneous Dirichlet BC
- **solve_lapack_dense()**: Simple dense solver
- **write_vtk()**: Save solution in VTK format

## API Design Principles

1. **Clean Types**: All user-facing types use `_t` suffix following Fortran conventions
2. **Factory Functions**: Constructor functions use `lowercase_underscore` naming
3. **Hidden Implementation**: Users never see low-level assembly or DOF management
4. **Mathematical Syntax**: Forms use natural mathematical notation where possible
5. **Automatic Memory Management**: Types handle their own allocation/deallocation

## Using the API

Users interact with FortFEM through the single main module:

```fortran
use fortfem
```

This provides access to the clean high-level API. The old complex module imports are no longer needed.

### Current Simple API Example
```fortran
type(mesh_2d_t) :: mesh
type(sparse_matrix_t) :: A
real(dp), allocatable :: u(:), f(:)

call create_unit_square_mesh(mesh, n=20)
call assemble_poisson_2d(mesh, A, f)
call apply_zero_bc(mesh, A, f)
call solve_lapack_dense(A, u, info)
```

### Future Form-based API
```fortran
type(mesh_t) :: mesh
type(function_space_t) :: Vh
type(trial_function_t) :: u
type(test_function_t) :: v

mesh = unit_square_mesh(32, 32)
Vh = function_space(mesh, "Lagrange", 1)
u = trial_function(Vh)
v = test_function(Vh)

a = inner(grad(u), grad(v))*dx
L = f*v*dx
call solve(a == L, uh, bc)
```

## Implementation Status

- ✅ **Current**: Simplified solver functions working with existing codebase
- 🚧 **In Progress**: Form algebra implementation (`inner()`, `grad()`, operators)
- 📋 **Planned**: Complete boundary condition system, automatic solver selection