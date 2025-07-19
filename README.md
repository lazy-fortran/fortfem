# FortFEM

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

## Documentation

See the [doc/](doc/) folder for detailed documentation and [DESIGN.md](DESIGN.md) for architecture overview.
