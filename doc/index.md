---
title: FortFEM Documentation
---

# FortFEM Documentation

Welcome to the FortFEM documentation! FortFEM is a modern Fortran finite element library designed for ease of use, inspired by FreeFEM and FEniCS.

## Getting Started

- [Quick Start Guide](quickstart.html) - Get up and running with FortFEM
- [Examples](examples/index.html) - Learn from example programs
- [Module Overview](modules.html) - Understand the library structure

## User Guide

- [Design Documentation](design/index.html) - Architecture and design decisions
- [API Reference](../modules.html) - Detailed module documentation
- [Source Files](../lists/files.html) - Browse source code

## Features

- **Natural mathematical notation** for expressing weak forms
- **Multiple element types**: P1, P2, RT0, Nédélec
- **Efficient sparse solvers**: Direct and iterative methods
- **Modern Fortran**: Object-oriented design with clear interfaces
- **Extensible architecture**: Easy to add new elements and problems

## Example Code

### Current Simple API
```fortran
program poisson_simple
    use fortfem
    implicit none
    
    type(mesh_2d_t) :: mesh
    type(sparse_matrix_t) :: A
    real(dp), allocatable :: u(:), f(:)
    
    ! Create mesh
    call create_unit_square_mesh(mesh, n=20)
    
    ! Assemble and solve -∆u = f
    call assemble_poisson_2d(mesh, A, f)
    call apply_zero_bc(mesh, A, f)
    call solve_lapack_dense(A, f, info)
    
    call write_vtk("solution.vtk", mesh, f)
end program
```

### Target FEniCS-style API
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
    type(form_t) :: a, L
    
    ! Create mesh and function space
    mesh = unit_square_mesh(32, 32)
    Vh = function_space(mesh, "Lagrange", 1)
    
    ! Define variational problem
    u = trial_function(Vh)
    v = test_function(Vh)
    f = constant(1.0_dp)
    
    a = inner(grad(u), grad(v))*dx
    L = f*v*dx
    
    ! Solve with boundary conditions
    bc = dirichlet_bc(Vh, 0.0_dp, boundary)
    call solve(a == L, uh, bc)
end program
```

## Contributing

FortFEM is open source and welcomes contributions! Visit our [GitHub repository](https://github.com/itpplasma/fortfem) to:
- Report issues
- Submit pull requests
- Join discussions

## License

FortFEM is released under the MIT License.