---
title: FortFEM Documentation
---

# FortFEM Documentation

FortFEM is a modern Fortran library for finite-element and boundary-element
methods. It includes scalar and vector-valued discretizations, open-boundary
operators, sparse solvers, and analytical verification examples.

## Getting Started

- [Quick Start Guide](quickstart.html) - Get up and running with FortFEM
- [Examples](examples/index.html) - Learn from example programs
- [Module Overview](modules.html) - Understand the library structure

## User Guide

- [Design Documentation](design/index.html) - Architecture and design decisions
- [Differentiation](design/differentiation.html) - JVP/VJP contracts,
  FortSym generation, Enzyme validation, and sparse adjoints
- [FCI parallel operator](design/fci_parallel_operator.html) - PARALLAX-aligned
  mapped gradients and conservative support divergences
- [FCI boundary-patch mortar](design/fci_boundary_patch_mortar.html) -
  conservative transfer between FCI and fitted or cut boundary traces
- [Mixed wave time step](design/mixed_wave_time.html) - structure-preserving
  first-order pressure/velocity and port-Hamiltonian updates
- [Quadratic average-vector-field step](design/quadratic_average_vector_field.html)
  - energy-preserving ideal update with analytical JVP/VJP actions
- [CGL pressure tensor](design/cgl_pressure_tensor.html) - generated
  gyrotropic tensor pressure with JVP/VJP actions
- [Material-surface flux](design/nonlinear_surface_flux.html) - application
  callbacks with tagged residual, ledger, JVP, and VJP assembly
- [Roadmap](https://github.com/lazy-fortran/fortfem/blob/main/ROADMAP.md) -
  Structure-preserving FEM, BEM, IGA, Fourier-FEM, and plasma-physics ingredients
- [API Reference](../lists/modules.html) - Detailed module documentation
- [Source Files](../lists/files.html) - Browse source code

## Features

- **Finite-element exterior calculus**: arbitrary-order triangular and
  tetrahedral H1, H(curl), H(div), and discontinuous families
- **Open boundaries**: planar, circular, and spherical DtN maps, PML, and
  scalar and Maxwell boundary-integral operators
- **FEM/BEM coupling**: scalar transmission formulations and compatible
  Nédélec/RWG Maxwell trace coupling
- **Shared numerics**: `fortnum` quadrature and special functions plus
  `fortsparse` CSC assembly, factorization, and solves
- **Reproducible examples**: analytical comparisons and fortplot-generated
  gallery artifacts
- **Independent FEM oracles**: the
  [interoperability benchmark protocol](interoperability_benchmarks.html)
  defines license-safe comparisons with FEniCSx, FreeFEM, and MFEM

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

FortFEM is open source and welcomes contributions! Visit our
[GitHub repository](https://github.com/lazy-fortran/fortfem) to:
- Report issues
- Submit pull requests
- Join discussions

## License

FortFEM is released under the MIT License.
