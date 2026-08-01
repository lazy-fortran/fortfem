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

The [curved Maxwell trace-to-flux contract](maxwell_curved_dtn.html) composes
EFIE, MFIE/RBC, and a primal trace mass matrix into a weak curved-surface DtN
map with matrix-free, JVP, and VJP actions. Its exact-curved torus wrapper
keeps the RWG and RBC trace conventions explicit; the same contract includes
caller-owned weak-to-point trace reconstruction for physical field evaluation.

The [manufactured toroidal Maxwell FEM--BEM solution](maxwell_torus_fem_bem_solution.html)
composes the same curved RWG trace with a volume Nédélec solve and checks the
coupled field/current state against a constant analytical (H(\mathrm{curl}))
field.

The [Maxwell FEM--BEM state differentiation contract](maxwell_fem_bem_state_ad.html)
provides the reusable complex implicit value/JVP/VJP map for concatenated
volume and surface unknowns. Assembly-specific geometry products compose into
this layer rather than differentiating through a factorization.

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

The [partition layout contract](partition_layout.html) records fixed
local-to-global IDs, owner/ghost metadata, and deterministic communicator-free
reductions as the serial boundary for future distributed assembly.

The [deterministic property-testing contract](property_testing.html) provides
per-case seeds, bounded samples, and a reproducible integer shrinker through
the ordinary `fo test` helper path.

The [cell-complex cycle-space contract](cell_cycle_basis.html) exposes the
real cycle and cocycle kernels while keeping face quotients, metric harmonic
representatives, and gauges in higher layers.

The [harmonic period normalization contract](harmonic_period_normalization.html)
maps a fixed metric-harmonic basis to caller-selected cycle periods or flux
units, with dense-solve JVP/VJP actions and no application-specific labels.

The [complex cycle-period constraint](period_constraints.html) contracts a
fixed oriented cycle basis with a complex edge field and provides the matching
JVP/VJP block for tree--cotree, harmonic, FEM/BEM, and IGA compositions.

The [tree-cotree gauge contract](tree_cotree_gauge.html) provides a
fixed-topology direct-solve reduction for curl--curl nullspaces, including the
control-mesh rule needed by high-order FEEC and IGA/mortar spaces.

The [broken and skeleton space layout contract](broken_skeleton_spaces.html)
provides deterministic ownership maps for broken H¹/H(curl)/H(div)/L² spaces
and shared scalar, tangential, or normal facet traces. Basis evaluation,
metric maps, and physical interface laws remain caller-owned.

The [broken FEEC sequence contract](broken_feec_sequence.html) embeds
cell-local gradient, curl, and divergence maps without introducing
inter-cell continuity, and exposes their exact-sequence compositions together
with value/JVP/VJP actions for DG, HDG, cut, and IGA composition.

The [fitted trace constraint contract](fitted_trace_constraint.html) assembles
the signed multiplier block for independently discretized plus/minus traces,
with complete value/JVP/VJP actions and no embedded constitutive jump law.

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

The [linear-response interchange contract](linear_response_interchange.html)
defines the modal metadata, complex block composition, provenance, and
real-part complex JVP/VJP convention for external ideal/resistive, vacuum, and
wall response adapters without importing their application physics.

The [generalized complex eigen-residual](generalized_eigen_residual.html)
provides the neutral `K u - lambda M u` modal residual and analytic JVP/VJP
actions for Fourier-FEM, curl--curl, GLISS-like, and other linear-response
clients without owning an eigensolver or application-specific interchange
format.

The [coupled field residual](coupled_field_residual.html) composes a
caller-owned rectangular field operator and constraint operator into one
residual with exact real JVP/VJP actions. It is the neutral assembly boundary
for multi-field FEM, BEM, DtN, PML, tensor, and interface clients.

The [incomplete-Cholesky contract](incomplete_cholesky.html) provides an
explicit SPD IC(0) factor/apply primitive for compatible iterative paths and a
fixed-pattern sparse IC(0) path; its companion sparse ILU(0) API is the
explicit nonsymmetric counterpart.

The [sparse ILUT contract](sparse_ilut.html) adds deterministic drop and
per-column fill control for nonsymmetric response blocks, with sparse fixed-
factor apply and value/JVP/VJP actions. Its row-oriented companion provides
the same CSC factor contract with O(n + nnz) construction work storage.

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

The [symmetric mixed-wave split](mixed_wave_strang.html) composes ideal
midpoint/Cayley factors as A(Δt/2)-B(Δt)-A(Δt/2), with analytical JVP/VJP
actions and independent energy, reversibility, finite-difference, and adjoint
checks.

The [shifted enriched space](shifted_enriched_space.html) builds a complete
scalar XFEM/GFEM basis-value matrix from base values and fixed level-set anchor
signs, with topology-event rejection and exact JVP/VJP actions.

The [shifted vector enriched space](shifted_vector_enriched_space.html) lifts
that constructor to component-by-basis-by-point vector values while leaving
Piola choice and the de Rham sequence policy explicit for H(curl), H(div), IGA,
or DG clients.

The [dissipative Cayley step](dissipative_cayley.html) provides the separate
energy-contracting mass/damping update and its analytic JVP/VJP actions for
resistive, viscous, thermal, and absorbing clients.

The [CGL pressure tensor](cgl_pressure_tensor.html) documents the generated
gyrotropic tensor constitutive block and its JVP/VJP contract.

The [tensor volume work contract](tensor_volume_work.html) assembles the
caller-owned `stress:grad(test)` residual and its geometry/tensor/weight
JVP/VJP actions for CGL, Maxwell, anisotropic, and elastic clients.

The [tensor diffusion matrix contract](tensor_diffusion_matrix.html) provides
the analogous tensor-weighted gradient pairing for anisotropic scalar,
elastic, and compatible field blocks, with analytic JVP/VJP actions.

The [mixed elasticity residual](mixed_elasticity_residual.html) provides the
neutral first-order `C sigma - E u` constitutive block and `D sigma - f`
equilibrium block, including complete analytic JVP/VJP actions for compatible,
TDNNS, Hellinger--Reissner, and anisotropic elasticity clients.

The [elasticity weak-symmetry constraint](elasticity_symmetry_constraint.html)
adds the independent `W sigma - g` multiplier/trace block with complete
fixed-topology JVP/VJP actions.

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
