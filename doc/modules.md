---
title: Module Overview
---

# FortFEM API Overview

FortFEM provides a clean, high-level API that hides implementation details while maintaining the power and flexibility of finite element methods.

## Current API Structure

### Main Module
- **fortfem_api**: Complete FEniCS-style API with types and functions

### Core Types
- **mesh_t**: Mesh container with topology and geometry
- **function_space_t**: Scalar finite element function spaces  
- **vector_function_space_t**: Vector finite element function spaces
- **function_t**: Scalar functions with nodal values
- **vector_function_t**: Vector functions with edge values
- **trial_function_t** / **test_function_t**: Symbolic functions for forms
- **vector_trial_function_t** / **vector_test_function_t**: Vector versions
- **form_expr_t**: Expression trees for weak forms
- **form_equation_t**: Equation a == L for solving
- **dirichlet_bc_t**: Scalar Dirichlet boundary conditions
- **vector_bc_t**: Vector boundary conditions

### Factory Functions
**Mesh Creation:**
- **unit_square_mesh(n)**: Create n×n uniform unit square mesh
- **rectangle_mesh(nx, ny, domain)**: Create rectangular mesh with bounds
- **unit_disk_mesh(resolution)**: Create unit disk mesh
- **mesh_from_boundary(boundary, resolution)**: Create mesh from boundary

**Function Spaces:**
- **function_space(mesh, family, degree)**: Create scalar FE space
- **vector_function_space(mesh, family, degree)**: Create vector FE space
- **function(space)**: Create function in space
- **trial_function(space)** / **test_function(space)**: Create symbolic functions

**Forms and Operators:**
- **inner(a, b)**: Inner product of expressions
- **grad(u)**: Gradient operator
- **curl(u)**: Curl operator for vector fields
- **dx**: Integration measure (global variable)

**Boundary Conditions:**
- **dirichlet_bc(space, value)**: Create scalar Dirichlet BC
- **vector_bc(space, values, bc_type)**: Create vector BC

### Solving and Visualization
- **solve(equation, uh, bc)**: Solve weak form equation
- **plot(uh, filename, title, colormap)**: Plot scalar function
- **plot(Eh, filename, title, plot_type)**: Plot vector field
- **plot(mesh, filename, title)**: Plot mesh triangulation

## API Design Principles

1. **Clean Types**: All user-facing types use `_t` suffix following Fortran conventions
2. **Factory Functions**: Constructor functions use `lowercase_underscore` naming
3. **Hidden Implementation**: Users never see low-level assembly or DOF management
4. **Mathematical Syntax**: Forms use natural mathematical notation where possible
5. **Automatic Memory Management**: Types handle their own allocation/deallocation

## Using the API

Users interact with FortFEM through the main module:

```fortran
use fortfem_api
```

This provides access to the complete FEniCS-style API.

### Complete Working Example
```fortran
program poisson_example
    use fortfem_api
    implicit none
    
    type(mesh_t) :: mesh
    type(function_space_t) :: Vh
    type(trial_function_t) :: u
    type(test_function_t) :: v
    type(function_t) :: uh, f
    type(dirichlet_bc_t) :: bc
    type(form_expr_t) :: a, L
    
    ! Create mesh and function space
    mesh = unit_square_mesh(20)
    Vh = function_space(mesh, "Lagrange", 1)
    
    ! Define trial and test functions
    u = trial_function(Vh)
    v = test_function(Vh)
    f = constant(1.0_dp)
    
    ! Define weak form
    a = inner(grad(u), grad(v))*dx
    L = f*v*dx
    
    ! Create solution function and boundary conditions
    uh = function(Vh)
    bc = dirichlet_bc(Vh, 0.0_dp)
    
    ! Solve and visualize
    call solve(a == L, uh, bc)
    call plot(uh, "solution.png", "Poisson Solution", "viridis")
    call plot(mesh, "mesh.png", "FEM Mesh")
    
end program
```

## Implementation Status

- ✅ **Implemented**: Complete FEniCS-style API with form algebra
- ✅ **Working**: Scalar Poisson problems with P1 Lagrange elements
- ✅ **Working**: Vector curl-curl problems with Nédélec edge elements
- ✅ **Working**: Mesh generation and visualization with fortplotlib
- 📋 **Planned**: More element types, advanced boundary conditions, parallel solvers