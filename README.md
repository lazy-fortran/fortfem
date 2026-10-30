# FortFEM

A Fortran finite-element and boundary-element library with a high-level API
inspired by FEniCS and FreeFEM. It provides the structure-preserving (FEEC)
discretisations, open-boundary operators and sparse solves needed by plasma
physics and electromagnetics codes.

Status: experimental research code; APIs may change without notice.

## Features

- Finite-element exterior calculus: arbitrary-order triangular and tetrahedral
  H1, H(curl) (Nédélec), H(div) (Raviart-Thomas, BDM) and discontinuous
  families, with Piola maps, global orientations and commuting projections.
- Weak forms written as expressions (`inner(grad(u), grad(v))*dx`), assembled to
  dense element matrices or sparse CSC matrices.
- Open boundaries: planar, circular and spherical Dirichlet-to-Neumann maps, PML,
  and Laplace, Helmholtz and Maxwell boundary-integral operators; FEM/BEM
  coupling.
- Sparse direct solves through [fortsparse](https://github.com/lazy-fortran/fortsparse)
  (native CPU backend, real and complex), quadrature and special functions
  through [fortnum](https://github.com/lazy-fortran/fortnum).
- Mesh generation, refinement and plotting through
  [fortplot](https://github.com/lazy-fortran/fortplot).
- JVP/VJP derivative actions for many operators, and a C/C++ API
  ([`include/fortfem.h`](include/fortfem.h)).

A more detailed list of what is verified and what is still experimental is in
[doc/status.md](doc/status.md).

## Building

Requirements: a Fortran 2018 compiler, BLAS/LAPACK, OpenMP and
[fpm](https://github.com/fortran-lang/fpm). fpm fetches fortnum, fortsparse and
fortplot.

```bash
fpm build
fpm test
fpm run --example simple_poisson
```

The CMake build produces the C/C++ shared library `fortfem::capi` (plus a few
standalone component libraries) and fetches its pinned dependencies:

```bash
cmake -S . -B build/cmake -DCMAKE_BUILD_TYPE=Release
cmake --build build/cmake
cmake --install build/cmake --prefix /desired/prefix
```

Downstream CMake projects use `find_package(fortfem CONFIG REQUIRED)` and link
`fortfem::capi`. CMake options and the sparse solver configuration are described
in [doc/building.md](doc/building.md).

## Usage

From [`example/simple_poisson`](example/simple_poisson/simple_poisson.f90),
solving -Δu = 1 on the unit square with homogeneous Dirichlet data:

```fortran
program simple_poisson
    use fortfem_core, only: dirichlet_bc_t, function_space_t, mesh_t, &
        unit_square_mesh
    use fortfem_feec, only: constant, dirichlet_bc, dx, form_expr_t, function, &
        function_space, function_t, grad, inner, operator(*), &
        operator(==), solve, test_function, test_function_t, trial_function, &
        trial_function_t
    use fortfem_kinds, only: dp
    use fortfem_plot, only: plot
    implicit none

    type(mesh_t) :: mesh
    type(function_space_t) :: Vh
    type(trial_function_t) :: u
    type(test_function_t) :: v
    type(function_t) :: f, uh
    type(dirichlet_bc_t) :: bc
    type(form_expr_t) :: a, L

    mesh = unit_square_mesh(20)
    Vh = function_space(mesh, "Lagrange", 1)
    u = trial_function(Vh)
    v = test_function(Vh)
    f = constant(1.0_dp)

    a = inner(grad(u), grad(v))*dx
    L = f*v*dx
    bc = dirichlet_bc(Vh, 0.0_dp)
    uh = function(Vh)

    call solve(a == L, uh, bc)
    call plot(uh, filename="poisson_solution.png", colormap="viridis")
end program simple_poisson
```

More programs, including H(curl) problems, DtN and BEM examples, are under
[`example/`](example); the
[example gallery](https://lazy-fortran.github.io/fortfem/page/examples/index.html)
shows their source and plots.

## Documentation

- Documentation site: <https://lazy-fortran.github.io/fortfem/> (generated with
  FORD from `doc/` and the sources)
- [Quick start](doc/quickstart.md), [module overview](doc/modules.md),
  [design notes](doc/design/index.md)
- [Building and sparse solver configuration](doc/building.md)
- [Implementation status](doc/status.md)

## Related projects

- [fortsparse](https://github.com/lazy-fortran/fortsparse): sparse direct solvers.
- [fortnum](https://github.com/lazy-fortran/fortnum): quadrature, special
  functions and other numerics.
- [fortsym](https://github.com/lazy-fortran/fortsym): symbolic generation of the
  committed reference-element kernels under `src/generated`.
- [fortplot](https://github.com/lazy-fortran/fortplot): plotting.

## License

MIT. See [LICENSE](LICENSE).
