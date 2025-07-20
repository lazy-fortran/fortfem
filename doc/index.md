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
- [API Reference](../lists/modules.html) - Detailed module documentation
- [Source Files](../lists/files.html) - Browse source code

## Features

- **Natural mathematical notation** for expressing weak forms
- **Multiple element types**: P1 Lagrange and Nédélec edge elements
- **Built-in visualization**: Plotting with fortplotlib integration
- **Modern Fortran**: Object-oriented design with clear interfaces
- **Clean API**: FEniCS-inspired interface for ease of use

## Example Code

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
    
    ! Define variational problem
    u = trial_function(Vh)
    v = test_function(Vh)
    f = constant(1.0_dp)
    
    a = inner(grad(u), grad(v))*dx
    L = f*v*dx
    
    ! Solve with boundary conditions
    uh = function(Vh)
    bc = dirichlet_bc(Vh, 0.0_dp)
    call solve(a == L, uh, bc)
    
    ! Visualize results
    call plot(uh, "solution.png", "Poisson Solution", "viridis")
    call plot(mesh, "mesh.png", "FEM Mesh")
end program
```

## Contributing

FortFEM is open source and welcomes contributions! Visit our [GitHub repository](https://github.com/itpplasma/fortfem) to:
- Report issues
- Submit pull requests
- Join discussions

## License

FortFEM is released under the MIT License.