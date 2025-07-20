[![CI](https://github.com/itpplasma/fortfem/actions/workflows/ci.yml/badge.svg)](https://github.com/itpplasma/fortfem/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/itpplasma/fortfem/branch/main/graph/badge.svg?token=CODECOV_TOKEN)](https://codecov.io/gh/itpplasma/fortfem)

A modern Fortran finite element library designed for ease of use, inspired by FreeFEM and FEniCS.

## Features

- **FEniCS-style API**: Natural mathematical notation inspired by FEniCS/FreeFEM
- **Simple Forms Syntax**: Define weak forms using `inner(grad(u), grad(v))*dx`
- **One-line plotting**: `call plot(uh, title="Solution", colormap="viridis")`
- Support for various element types:
  - **P1 Lagrange**: Linear elements with optimal convergence
  - **Nédélec edge elements**: H(curl) conforming for electromagnetic problems
- **Vector problems**: Curl-curl equations with GMRES iterative solver
- **Built-in visualization**: Automatic plotting with fortplotlib integration
- **Minimal code**: Solve PDEs in ~10 lines of meaningful code
- Test-driven development with comprehensive test suite

## Quick Start

```bash
# Build the library
fpm build

# Run tests
fpm test

# Run examples
fpm run --example simple_poisson
```

## Usage Example

FortFEM provides a clean, FEniCS-inspired API for defining finite element problems:

```fortran
program poisson_example
    use fortfem_api
    
    ! Create mesh and function space
    mesh = unit_square_mesh(20)
    Vh = function_space(mesh, "Lagrange", 1)
    
    ! Define trial and test functions
    u = trial_function(Vh)
    v = test_function(Vh)
    f = constant(1.0_dp)
    
    ! Define weak form using natural mathematical notation
    a = inner(grad(u), grad(v))*dx  ! Bilinear form: ∫ ∇u·∇v dx
    L = f*v*dx                      ! Linear form:   ∫ f v dx
    
    ! Solve and plot in one line
    uh = function(Vh)
    bc = dirichlet_bc(Vh, 0.0_dp)
    call solve(a == L, uh, bc)
    call plot(uh, title="Poisson Solution", colormap="viridis")
end program
```

## Examples

Explore the [examples/](https://github.com/itpplasma/fortfem/tree/main/example) directory for complete working examples:

- [**Simple Poisson solver**](https://github.com/itpplasma/fortfem/blob/main/example/simple_poisson.f90) - FEniCS-style API demonstration with plotting
- [**Curl-curl electromagnetic**](https://github.com/itpplasma/fortfem/blob/main/example/curl_curl_example.f90) - Vector problems with Nédélec elements and GMRES solver  
- [**Plotting demonstration**](https://github.com/itpplasma/fortfem/blob/main/example/plotting_demo.f90) - Comprehensive plotting API showcase

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