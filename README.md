[![CI](https://github.com/itpplasma/fortfem/actions/workflows/ci.yml/badge.svg)](https://github.com/itpplasma/fortfem/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/itpplasma/fortfem/branch/main/graph/badge.svg?token=CODECOV_TOKEN)](https://codecov.io/gh/itpplasma/fortfem)

A modern Fortran finite element library designed for ease of use, inspired by FreeFEM and FEniCS.

## Features

- **FEniCS-style API**: Natural mathematical notation inspired by FEniCS/FreeFEM
- **Simple Forms Syntax**: Define weak forms using `inner(grad(u), grad(v))*dx`
- Support for various element types:
  - **P1 Lagrange**: Linear elements with optimal convergence
  - **P2 Lagrange**: Quadratic elements with O(h³) L2 convergence
  - **Nédélec edge elements**: H(curl) conforming for electromagnetic problems
- Triangular and rectangular meshes with full connectivity
- High-order Gaussian quadrature rules (up to order 7)
- Built-in visualization with fortplotlib
- Test-driven development with comprehensive test suite

## Quick Start

```bash
# Build the library
fpm build

# Run tests
fpm test

# Run examples
fpm run --example plot_basis
```

## Usage Example

FortFEM provides a clean, FEniCS-inspired API for defining finite element problems:

```fortran
program poisson_example
    use fortfem
    
    ! Create mesh and function space
    mesh = unit_square_mesh(32)
    Vh = function_space(mesh, "Lagrange", 1)
    
    ! Define trial and test functions
    u = trial_function(Vh)
    v = test_function(Vh)
    
    ! Define weak form: inner(grad(u), grad(v))*dx
    a = inner(grad(u), grad(v))*dx
    L = inner(f, v)*dx
    
    ! Solve the system
    uh = function(Vh)
    bc = dirichlet_bc(Vh, 0.0_dp)
    solve(a == L, uh, bc)
end program
```

## Examples

Explore the [examples/](https://github.com/itpplasma/fortfem/tree/main/example) directory for complete working examples:

- [Poisson equation solver](https://github.com/itpplasma/fortfem/blob/main/example/poisson_2d.f90)
- [Nédélec element visualization](https://github.com/itpplasma/fortfem/blob/main/example/plot_basis.f90)
- [Mesh generation and visualization](https://github.com/itpplasma/fortfem/blob/main/example/mesh_2d_demo.f90)

## Project Structure

- [`src/`](https://github.com/itpplasma/fortfem/tree/main/src) - Core library modules
- [`test/`](https://github.com/itpplasma/fortfem/tree/main/test) - Comprehensive test suite
- [`example/`](https://github.com/itpplasma/fortfem/tree/main/example) - Example programs
- [`doc/`](https://github.com/itpplasma/fortfem/tree/main/doc) - Documentation
- [`app/`](https://github.com/itpplasma/fortfem/tree/main/app) - Main applications

## Navigation

- [API Documentation](https://itpplasma.github.io/fortfem/modules.html) - Detailed module and procedure documentation
- [Source Files](https://itpplasma.github.io/fortfem/sourcefile/index.html) - Browse the source code
- [Program Structure](https://itpplasma.github.io/fortfem/program/index.html) - Main programs and examples
- [Design Documentation](https://itpplasma.github.io/fortfem/page/design/index.html) - Architecture and design decisions

## Quick Links

- [Getting Started](https://itpplasma.github.io/fortfem/page/quickstart.html) - Build and run your first FE simulation
- [Examples](https://itpplasma.github.io/fortfem/page/examples/index.html) - Complete working examples
- [Module Reference](https://itpplasma.github.io/fortfem/modules.html) - Complete API reference

## Contributing

1. Check the [TODO.md](https://github.com/itpplasma/fortfem/blob/main/TODO.md) for current development priorities
2. Follow strict TDD: write tests first
3. See [CLAUDE.md](https://github.com/itpplasma/fortfem/blob/main/CLAUDE.md) for development guidelines

## License

See [LICENSE](https://github.com/itpplasma/fortfem/blob/main/LICENSE) for details.