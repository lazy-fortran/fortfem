---
title: Design Documentation
---

# FortFEM Design Documentation

## Overview

FortFEM is designed with the following principles:
- **Natural mathematical notation**: Express weak forms as you would write them mathematically
- **Modularity**: Clear separation between mesh, function spaces, assembly, and solvers
- **Performance**: Efficient sparse matrix operations and assembly routines
- **Extensibility**: Easy to add new element types and problem formulations

The [differentiation strategy](differentiation.html) defines the common
primal/JVP/VJP API, analytical FortSym path, optional Enzyme tournament, and
implicit sparse-solve adjoints.

The [spherical-harmonic contract](spherical_harmonics.html) records the
FortNum-backed normalization, angular domains, and pole convention used by
Fourier-FEM and boundary operators.

The [toroidal-harmonic contract](toroidal_harmonics.html) records the
FortNum-backed Hobson-normalized half-integer P/Q branches used by exact
toroidal analytical solutions and DtN fixtures.

The [FCI parallel support-operator contract](fci_parallel_operator.html)
defines the PARALLAX-aligned mapped gradient and its conservative
volume-weighted adjoint without coupling FortFEM to a particular field-line
tracer.

The [FCI interpolation-map contract](fci_interpolation_map.html) defines
linear, bilinear, quadratic, and FortSym-generated cubic and quartic
fixed-stencil maps,
including their topology-aware JVP/VJP actions and independent polynomial
oracles.

The [FCI boundary-patch mortar contract](fci_boundary_patch_mortar.html)
defines conservative cross-mass transfer between an FCI background trace and a
fitted or cut boundary patch, including constant reproduction, ownership, and
the weighted adjoint identity.

The [oriented cell-complex contract](cell_complex.html) defines integer chain
boundary maps, exact boundary-of-boundary validation, Euler characteristic, and
small homology diagnostics for later FEEC, gauge, cut, and interface graphs.

The [cell-complex cycle-space contract](cell_cycle_basis.html) exposes the
real cycle and cocycle kernels while keeping face quotients, metric harmonic
representatives, and gauges in higher layers.

The [harmonic period normalization contract](harmonic_period_normalization.html)
maps a fixed metric-harmonic basis to caller-selected cycle periods or flux
units, with dense-solve JVP/VJP actions and no application-specific labels.

The [tree-cotree gauge contract](tree_cotree_gauge.html) provides a
fixed-topology direct-solve reduction for curl--curl nullspaces, including the
control-mesh rule needed by high-order FEEC and IGA/mortar spaces.

The [FEEC exact-sequence diagnostic](feec_exact_sequence.html) exposes the
metric-independent `curl(grad)` and `div(curl)` composition defects with
value/JVP/VJP actions for simplicial, IGA, multipatch, and periodic incidence
maps.

The [FEEC commuting-projection diagnostic](feec_commuting_projection.html)
compares discrete differentials with projected continuous differentials and
provides complete value/JVP/VJP actions for fitted, cut, enriched, IGA, and
periodic maps.

The [Fourier mode registry](fourier_mode_registry.html) defines fixed-topology
field-period phase, normalization, conjugate packing, triad lookup, radial
regularity, and complex-coordinate derivative conventions for Fourier-FEM and
IGA clients.

The [neutral equilibrium interchange contract](equilibrium_interchange.html)
defines the sample/provenance boundary for external equilibrium adapters
without implementing GEQDSK, COCOS, or plasma-specific readers and closures.

The [incomplete-Cholesky contract](incomplete_cholesky.html) provides an
explicit SPD IC(0) factor/apply primitive for compatible iterative paths and a
fixed-pattern sparse IC(0) path; its companion sparse ILU(0) API is the
explicit nonsymmetric counterpart.

The [region/interface graph contract](region_interface_graph.html) adds
oriented plus/minus region incidence, periodic self-identifications, and
connectivity labels without importing an application-specific interface law.

The [internal-manifold graph contract](internal_manifold_graph.html) adds
explicit open/closed manifold endpoints and junction incidence for later
surface-current divergence and interface-balance operators.

The [signed cell-identification contract](cell_identification.html) records
canonical representatives and orientation signs for quotient or periodic
metadata without constructing a mesh or interpreting application coordinates.

The [quotient cell-complex contract](quotient_cell_complex.html) composes
those signed maps with integer boundary matrices and validates the resulting
periodic or multipatch chain complex before metric assembly.

The [surface-current trace contract](surface_current.html) provides generic
Ampere jump algebra, an integrated current ledger, and fixed-topology JVP/VJP
actions without embedding a material or plasma boundary law.

The same [surface-current trace contract](surface_current.html) documents the
geometry-to-edge-flux contraction and its analytic JVP/VJP, which composes
caller-owned fitted, cut, DG, or IGA edge quadrature with the topology-only
surface-edge balance.

The surface-current trace contract also exposes a neutral residual for an
independent tangential-current test/trial basis, with analytic JVP/VJP actions
for fitted, cut, DG/HDG, and IGA trace ownership.

The [mixed first-order wave step](mixed_wave_time.html) documents the common
pressure/velocity, displacement/momentum, and electromagnetic midpoint/Cayley
contract with independent energy and reversibility checks.

The [dissipative Cayley step](dissipative_cayley.html) provides the separate
energy-contracting mass/damping update and its analytic JVP/VJP actions for
resistive, viscous, thermal, and absorbing clients.

The [CGL pressure tensor](cgl_pressure_tensor.html) documents the generated
gyrotropic tensor constitutive block and its JVP/VJP contract.

The [tensor volume work contract](tensor_volume_work.html) assembles the
caller-owned `stress:grad(test)` residual and its geometry/tensor/weight
JVP/VJP actions for CGL, Maxwell, anisotropic, and elastic clients.

The [compatible flux elimination contract](compatible_flux_elimination.html)
provides the differentiable local Schur reduction and recovery map for
RT/BDM/Nédélec, compatible IGA, and HDG flux blocks without owning a global
skeleton numbering.

The [interface normal-traction balance](interface_traction_balance.html)
documents the neutral traction-jump residual that composes tensor pressure,
elastic stress, or Maxwell-stress blocks without selecting an application law.

The [shifted Heaviside enrichment](heaviside_enrichment.html) documents the
first XFEM/GFEM partition-of-unity enrichment and its explicit topology-event
differentiation contract, including scalar and componentwise vector bases for
H(curl)/H(div) composition and the corrected-XFEM blending ramp.

The [enrichment-support activation contract](enrichment_support_activation.html)
defines CSR support sign activation, unique extrema and signed margins, with
fixed-topology JVP/VJP actions and explicit topology-event rejection.

The [vector enrichment-support conditioning contract](enrichment_vector_support_diagnostics.html)
defines the physical vector/metric Gram contraction for Piola-mapped FEEC and
IGA enrichments, with geometry-compatible JVP/VJP actions and SPD rejection.

The [Piola-enriched vector composition contract](piola_enriched_vector.html)
maps 2D/3D covariant or contravariant reference vectors before applying a
pointwise enrichment, including geometry/reference/activation JVP and VJP
actions.

The [nonlinear material-surface flux contract](nonlinear_surface_flux.html)
keeps application wall and sheath laws separate from orientation-preserving
trace assembly and its generated-compatible value/JVP/VJP bookkeeping.

## Architecture

### Core Components

1. **Mesh Module** (`fortfem_mesh_*`)
   - Handles mesh data structures and connectivity
   - Supports triangular and quadrilateral elements
   - Provides mesh I/O and generation utilities

2. **Function Space Module** (`fortfem_function_space`)
   - Defines finite element spaces on meshes
   - Manages degrees of freedom (DOF) numbering
   - Handles boundary condition identification

3. **Basis Functions** (`fortfem_basis_*`)
   - Implements shape functions and their derivatives
   - Supports P1/P2/Q1 and verified lowest-order Nédélec/RT0 kernels
   - Handles reference-to-physical element transformations

4. **Assembly Module** (`fortfem_assembly_*`)
   - Assembles global matrices and vectors from element contributions
   - Uses efficient quadrature rules
   - Supports various bilinear and linear forms

5. **Solver Module** (`fortfem_solver`)
   - Interfaces to LAPACK and includes custom GMRES implementation
   - Provides both direct and iterative methods
   - Includes GMRES implementation

## Design Patterns

### Object-Oriented Approach
FortFEM uses Fortran's object-oriented features:
```fortran
type :: mesh_2d_t
    real(dp), allocatable :: vertices(:,:)
    integer, allocatable :: triangles(:,:)
    integer :: nvert, ntri
contains
    procedure :: init => mesh_2d_init
    procedure :: destroy => mesh_2d_destroy
end type
```

### Generic Programming
The library uses interfaces for flexibility:
```fortran
interface assemble
    module procedure assemble_mass_p1
    module procedure assemble_stiffness_p1
    module procedure assemble_mass_p2
end interface
```

## Element Types

### Lagrange Elements
- **P1**: Linear polynomials on triangles (implemented)
- **P2**: Quadratic polynomials with 6 DOFs per triangle (implemented)

### Vector Elements

- **Nédélec**: First- and second-kind triangular families through order four,
  with covariant Piola maps and oriented sparse curl-mass forms.
- **Tetrahedral Nédélec**: First-kind reference bases through order five have
  verified edge, face, and cell moments and curls. The affine covariant map
  and local curl-mass assembly accept all implemented orders. Canonical edge
  signs and generated/runtime face transforms provide orientation-aware
  global topology and direct `fortsparse` CSC assembly through order five.
- **Raviart--Thomas and BDM**: Triangular families through order four, with
  contravariant Piola maps and oriented sparse divergence-mass forms.
- Weighted vector forms compile to local matrices and `fortsparse` CSC
  matrices. The order-one Nédélec solve supports cellwise scalar curl and mass
  coefficients, constant or cellwise vector sources, constant physical
  tangential data, owned nonconstant edge moments, and a sparse direct
  interior solve. Its three-material Fourier-mode magnetic convergence has an
  analytical oracle, and a smooth manufactured curl-mass problem attains
  first-order field convergence. General pointwise sources, tensor
  coefficients, and higher-order high-level solves remain experimental.

## Weak Form Framework

The design allows natural expression of weak forms:
```fortran
! Poisson equation: -∆u = f
! Weak form: ∫∇u·∇v dx = ∫fv dx
call assemble_stiffness(mesh, V, A)
call assemble_load(mesh, V, f, b)
```

## Memory Management

- Uses allocatable arrays for dynamic memory
- Automatic cleanup through finalizers
- Efficient sparse storage formats (CSR/CSC)

## Future Directions

- 3D element support
- Adaptive mesh refinement
- Parallel assembly and solvers
- More element types (DG, HDG)
