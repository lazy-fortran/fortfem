---
title: Implementation status
---

# Implementation status

FortFEM is experimental. This page summarises which numerical paths have
analytical or independent tests on the default branch and which interfaces are
still incomplete. The design notes under [design pages](https://lazy-fortran.github.io/fortfem/page/design/index.html) give the
details of individual components.

## Tested

### Elements and spaces

- P1 and P2 triangular Lagrange elements, Q1 quadrilaterals and scalar Poisson
  solvers.
- Triangles: arbitrary-order Lagrange, first- and second-kind Nédélec,
  Raviart-Thomas, BDM and discontinuous scalar families, tested through order
  four. Affine Piola maps, global moment orientations, commuting projections and
  the sparse differential and mass forms have analytical tests.
- Tetrahedra: first-kind Nédélec reference bases have no fixed order limit; the
  solve, topology, orientation and sparse assembly paths are tested through
  order six. Raviart-Thomas H(div) and discontinuous L2 spaces are tested
  through degree six, including an RT6-DG6 mixed Poisson solve with
  conservation tests.
- Reference-element tables are generated with fortsym through order four. Order
  five uses monomial moments computed at run time, and higher orders use
  fortnum's orthogonal simplex (Dubiner/Koornwinder) modes.

### Forms and solves

- Symbolic P1 scalar mass, stiffness and constant-load forms. Weighted H(curl)
  and H(div) forms compile to element matrices or fortsparse CSC matrices for
  the four triangular vector families.
- Sparse RT0-DG0 and matching-degree RT-DG mixed Poisson solves: the source
  balance holds on every cell, RT0-DG0 converges at first order and RT1-DG1 at
  second order in flux and pressure.
- Tensor-weighted tetrahedral assembly with physical vector loads reproduces the
  reference value of the magnetic-field benchmark from
  [CPC 277 (2022) 108401](https://doi.org/10.1016/j.cpc.2022.108401) at its
  reference resolution, and the reconstructed field converges toward an
  independent series solution.
- Iterative dense solvers and real or complex sparse direct solves through
  fortsparse.

### Open boundaries

- Planar and circular outgoing Helmholtz DtN operators, dense two-dimensional
  Laplace and Helmholtz boundary operators, a Helmholtz CFIE solve and symmetric
  Laplace FEM/BEM transmission coupling.

### Other

- A C/C++ API for oriented triangle meshes, retained complex sparse factors and
  RT0 coefficient transfer.
- Mesh generation, refinement and plotting through fortplot.

## Experimental or missing

- Symbolic mixed forms, nonconstant scalar loads, boundary measures and general
  pointwise vector sources are not compiled yet.
- The high-level order-one Nédélec path supports cellwise scalar curl and mass
  coefficients, constant or cellwise vector sources and prescribed tangential
  edge moments. It converges to an analytical three-material Fourier-mode
  magnetic solution and at first order on a manufactured curl-mass problem.
- Not yet available: two-dimensional tensor and complex coefficients,
  higher-order high-level solves, other three-dimensional FEEC families and
  electromagnetic BEM coupling.

## Generated kernels

The normal library build needs no computer algebra system. Generators live in
`tools/codegen`, pin their tested fortsym revision in
`tools/codegen/fortsym.lock` and commit their output under `src/generated`.
`tools/codegen/check_generated.sh` regenerates into a temporary directory and
compares the result byte for byte with the committed kernels.
