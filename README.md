[![CI](https://github.com/itpplasma/fortfem/actions/workflows/ci.yml/badge.svg)](https://github.com/itpplasma/fortfem/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/itpplasma/fortfem/branch/main/graph/badge.svg?token=CODECOV_TOKEN)](https://codecov.io/gh/itpplasma/fortfem)

A modern Fortran finite element library designed for ease of use, inspired by FreeFEM and FEniCS.

## Features

- Natural mathematical notation for weak form definition
- Support for various element types:
  - **P1 Lagrange**: Linear elements with optimal convergence
  - **P2 Lagrange**: Quadratic elements with O(h³) L2 convergence
  - **Nédélec edge elements**: H(curl) conforming for electromagnetic problems
- Triangular meshes with edge connectivity
- High-order Gaussian quadrature rules (up to order 7)
- Built-in visualization with fortplotlib
- Test-driven development with >90% coverage

## Quick Start

```bash
# Build the library
fpm build

# Run tests
fpm test

# Run examples
fpm run --example plot_basis
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

## Documentation

See the [doc/](https://github.com/itpplasma/fortfem/tree/main/doc) folder for detailed documentation and [DESIGN.md](https://github.com/itpplasma/fortfem/blob/main/DESIGN.md) for architecture overview.

## Contributing

1. Check the [TODO.md](https://github.com/itpplasma/fortfem/blob/main/TODO.md) for current development priorities
2. Follow strict TDD: write tests first
3. See [CLAUDE.md](https://github.com/itpplasma/fortfem/blob/main/CLAUDE.md) for development guidelines

## License

See [LICENSE](https://github.com/itpplasma/fortfem/blob/main/LICENSE) for details.
