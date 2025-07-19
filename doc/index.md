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

```fortran
program poisson
    use fortfem
    implicit none
    
    type(mesh_2d_t) :: mesh
    type(p1_space_t) :: V
    type(sparse_matrix_t) :: A
    real(dp), allocatable :: b(:), u(:)
    
    ! Create mesh and function space
    call create_rectangular_mesh(mesh, nx=20, ny=20)
    call V%init(mesh)
    
    ! Assemble and solve -∆u = f
    call assemble_poisson_2d(mesh, V, A, b)
    call apply_dirichlet_bc(A, b, V%boundary_dofs)
    call solve_sparse(A, b, u)
end program
```

## Contributing

FortFEM is open source and welcomes contributions! Visit our [GitHub repository](https://github.com/itpplasma/fortfem) to:
- Report issues
- Submit pull requests
- Join discussions

## License

FortFEM is released under the MIT License.