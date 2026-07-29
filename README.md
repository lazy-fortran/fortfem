[![CI](https://github.com/itpplasma/fortfem/actions/workflows/ci.yml/badge.svg)](https://github.com/itpplasma/fortfem/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/itpplasma/fortfem/branch/main/graph/badge.svg?token=CODECOV_TOKEN)](https://codecov.io/gh/itpplasma/fortfem)
[![Documentation](https://github.com/itpplasma/fortfem/actions/workflows/docs.yml/badge.svg)](https://itpplasma.github.io/fortfem/)

> **Note**: This project is experimental and subject to major changes. APIs may change without notice.

A modern Fortran finite element library designed for ease of use, inspired by FreeFEM and FEniCS.

## Implementation Status

Verified numerical paths:

- P1 and P2 triangular Lagrange elements, Q1 quadrilaterals, and scalar
  Poisson solvers.
- Arbitrary-order triangular Lagrange, first- and second-kind Nédélec,
  Raviart--Thomas, BDM, and discontinuous scalar families through order four.
  Their affine Piola maps, global moment orientations, commuting projections,
  and sparse differential or mass forms have analytical tests.
- Iterative dense solvers and real or complex sparse direct solves through
  `fortsparse`.
- A C/C++ API for oriented triangle meshes, retained complex sparse factors,
  and RT0 coefficient transfer.
- Planar and circular outgoing Helmholtz DtN operators, dense two-dimensional
  Laplace and Helmholtz boundary operators, a Helmholtz CFIE solve, and
  symmetric Laplace FEM/BEM transmission coupling.
- Executable symbolic P1 scalar mass, stiffness, and constant-load forms.
  Weighted H(curl) and H(div) forms compile to element matrices or
  `fortsparse` CSC matrices for all four triangular vector families.
- A sparse RT0-DG0 mixed Poisson solve with homogeneous Dirichlet pressure
  data preserves the source balance independently on every cell.
- Mesh generation, refinement, and plotting through `fortplot`.

Experimental interfaces:

- Symbolic mixed forms, nonconstant scalar loads, boundary measures, and
  general pointwise vector sources are not compiled yet.
- The high-level order-one Nédélec path compiles cellwise scalar curl and mass
  coefficients, constant or cellwise physical vector sources, and constant
  physical or explicitly supplied tangential edge moments. Its direct solve
  uses a `fortsparse` interior block and converges to the analytical
  three-material, Fourier-mode magnetic solution. A smooth manufactured
  curl-mass problem attains first-order field convergence.

Tensor and complex coefficients, higher-order high-level solves,
three-dimensional elements, electromagnetic BEM coupling, and validation
against the unavailable proprietary MEPHIT cases and the not-yet-ported
acoustic-paper fixtures remain. See the
[FEEC, MEPHIT, and open-boundary roadmap](doc/roadmap_mephit_feec_bem.md).

## Quick Start

```bash
# Build the library
fpm build

# Run tests
fpm test

# Run examples
fpm run --example simple_poisson
```

For C and C++ consumers, CMake builds the focused shared API and fetches the
pinned `fortsparse` dependency when needed:

```bash
cmake -S . -B build/cmake -DCMAKE_BUILD_TYPE=Release
cmake --build build/cmake
cmake --install build/cmake --prefix /desired/prefix
```

Downstream CMake projects can then use `find_package(fortfem CONFIG REQUIRED)`
and link `fortfem::capi`. The installed
[`fortfem.h`](include/fortfem.h) documents the zero-based array layouts,
process-local handles, and status arguments.

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

- [**Simple Poisson solver**](https://github.com/itpplasma/fortfem/blob/main/example/simple_poisson/simple_poisson.f90) - FEniCS-style API demonstration with plotting
- [**Curl-curl electromagnetic prototype**](https://github.com/itpplasma/fortfem/blob/main/example/curl_curl/curl_curl.f90) - Experimental vector API and plotting path
- [**Plotting demonstration**](https://github.com/itpplasma/fortfem/blob/main/example/plotting/plotting.f90) - Plotting API examples
- [**Mesh plotting**](https://github.com/itpplasma/fortfem/blob/main/example/plot_mesh/plot_mesh.f90) - Mesh visualization examples

## Project Structure

- [`src/`](https://github.com/itpplasma/fortfem/tree/main/src) - Core library modules
- [`test/`](https://github.com/itpplasma/fortfem/tree/main/test) - Comprehensive test suite
- [`example/`](https://github.com/itpplasma/fortfem/tree/main/example) - Example programs
- [`doc/`](https://github.com/itpplasma/fortfem/tree/main/doc) - Documentation
- [`app/`](https://github.com/itpplasma/fortfem/tree/main/app) - Main applications

## Navigation

- [API Documentation](https://itpplasma.github.io/fortfem/lists/modules.html) - Detailed module and procedure documentation
- [Source Files](https://itpplasma.github.io/fortfem/lists/files.html) - Browse the source code
- [Program Structure](https://itpplasma.github.io/fortfem/lists/procedures.html) - Main programs and procedures
- [Design Documentation](https://itpplasma.github.io/fortfem/page/design/index.html) - Architecture and design decisions

## Quick Links

- [Getting Started](https://itpplasma.github.io/fortfem/page/quickstart.html) - Build and run your first FE simulation
- [Examples](https://itpplasma.github.io/fortfem/page/examples/index.html) - Complete working examples
- [Module Reference](https://itpplasma.github.io/fortfem/lists/modules.html) - Complete API reference

## Contributing

1. Check the [GitHub Issues](https://github.com/itpplasma/fortfem/issues) for current development priorities
2. Follow strict TDD: write tests first
3. See [CLAUDE.md](https://github.com/itpplasma/fortfem/blob/main/CLAUDE.md) for development guidelines

## License

See [LICENSE](https://github.com/itpplasma/fortfem/blob/main/LICENSE) for details.
