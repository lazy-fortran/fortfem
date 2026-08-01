# FortFEM roadmap

Status: living architecture and verification plan, 2026-08-01

FortFEM is a Fortran library for finite-element, boundary-element, and
compatible discretizations. The long-term goal is to provide the reusable
mathematical ingredients from which three-dimensional equilibrium, linear
stability, and nonlinear resistive-MHD applications can be built. This file
is a roadmap for those ingredients. It is not a promise to reimplement
VMEC, GVEC, SPEC, GPEC, MARS, GLISS, JOREK, CHEASE, FreeGS, or any other
production code in this repository.

The target architecture is:

\[
\boxed{\text{FEEC}+\text{IGA}+\text{Fourier-FEM}+\text{DG/HDG}
      +\text{XFEM/XIGA}+\text{explicit interfaces}
      +\text{distributional sources}+\text{FortSym}
      +\text{differentiable residuals}+\text{structure preservation}}
\]

The public result should be a small, composable library with many toy
problems, analytical solutions, and independent oracles. A future equilibrium
or MHD application can then select the appropriate spaces, geometry, weak
forms, solvers, and time integrator without changing the numerical contracts.

## FortFEM and application boundaries

FortFEM owns numerical foundations. It does not implement plasma applications.
In particular, the repository will not contain GEQDSK or COCOS readers and
writers, equilibrium profile models, coil or vessel physics, Braginskii or
sheath closures, species and reaction networks, neutral or impurity models,
or production CHEASE, FreeGS, VMEC, GVEC, SPEC, GPEC, MARS, GLISS, JOREK,
GRILLIX, BOUT++, or PARALLAX algorithms.

External applications may provide those data and laws through documented
callbacks or adapters. FortFEM supplies the typed geometry, field, trace,
tensor, Fourier, residual, derivative, solver, and structure-preserving
contracts that consume them. License-safe examples use manufactured fields,
small neutral arrays, or data produced by an external adapter. A reference
code name in this roadmap identifies a parity target or oracle, not a planned
FortFEM module.

## How to read this document

The labels below describe the state of a deliverable.

| Label | Meaning |
| --- | --- |
| **done** | Present on `main` and covered by an independent behavioral oracle |
| **active** | Partly present or currently being closed with tests first |
| **planned** | Required architecture or example that has not been implemented |
| **reference** | External code or literature target, not a FortFEM implementation |

The detailed current audit is in [the FEEC, MEPHIT, and open-boundary
roadmap](doc/roadmap_mephit_feec_bem.md). Derivative contracts are specified in
[the differentiation design](doc/design/differentiation.md), and the existing
JOREK-related coverage is tracked in
[the JOREK operator document](doc/jorek_iga_operator_coverage.md).

## 1. Development principles

These principles apply to every phase. They are part of the product design,
not only contributor etiquette.

### 1.1 Mathematical and scientific principles

1. **Structure preservation is a first-class requirement.** Exact sequences,
   orientations, conservation laws, gauge structure, symplectic or Poisson
   structure, and correct dissipation must be represented explicitly. A result
   that looks plausible but violates the relevant invariant is not accepted as
   a reference implementation.
2. **Residuals come first.** Every formulation is an explicit residual
   \(R(u,p)=0\), with declared fields, spaces, traces, constraints, and
   boundary terms. The assembled matrix is a product of this residual, not a
   separate hand-written equation.
3. **FortSym is the source of truth wherever algebra is repetitive.** Weak-form
   reduction, tensor contractions, element kernels, enrichment derivatives,
   manufactured forcing, and analytical comparison data should be derived and
   generated. Generated files are revision-pinned, checked byte for byte, and
   never edited by hand.
4. **Independent oracles are mandatory.** A test must compare against an
   analytical solution, a separately written reference expression, a discrete
   identity, a conservation law, a cross-formulation result, or an external
   code. Checking that a file matches its own generator is a gate, not a
   behavioral test.
5. **Differentiation is part of the API.** New operators expose primal, JVP,
   and VJP actions with documented real and complex inner-product conventions.
   Shape, coefficient, geometry, and interface parameters are differentiated
   when the topology is fixed. Topology changes report piecewise-smooth or
   invalid derivative regions instead of silently returning a misleading
   gradient.
6. **Compatible spaces before stabilizing patches.** Use de Rham complexes,
   commuting projections, Piola maps, and orientation-aware traces before
   adding divergence cleaning or penalty terms. Stabilization must have a
   declared consistency and conservation effect.
7. **The code remains Fortran-first.** Do not move performance-critical work
   to C. Existing interoperability ABIs may remain where they are already
   required, but new numerical kernels are Fortran and generated Fortran.
8. **No Lean dependency is required.** FortSym identities, independent
   numerical oracles, exact-sequence checks, compiler diagnostics, and Enzyme
   checks are the normal proof and verification stack. A formal proof tool can
   be introduced later for a narrowly defined theorem if it provides a clear
   benefit, but it is not a prerequisite for the roadmap.

### 1.2 Engineering principles

- Make the smallest complete change. Add the test and oracle before the
  implementation.
- Keep focused tests short enough for roughly a ten-second feedback loop.
  `fo test` therefore keeps a ten-second per-target wall-clock limit. Tests
  that intentionally exercise large tetrahedral systems, high-order
  convergence, or curved torus operators carry a `_slow` suffix and are run
  explicitly with `fo test --all`; they are not silently allowed to consume the
  fast-suite budget. Run full `fo` suites and expensive examples asynchronously
  in CI or a controlled background job. The build-test workflow uses the same
  `fo test` runner for its per-test timeout and coverage path, rather than an
  unbounded serial `fpm test`. Do not wait for a complete gallery build before
  making an independent, low-risk progress increment.
- Work directly on FortFEM `main` because it is the upstream library. Work on
  PR branches for downstream repositories such as MEPHIT and FortPlot. Never
  duplicate a fix that has already landed upstream. Rebase downstream work on
  the current upstream branch before opening or updating a PR.
- Commit coherent increments, push them promptly, and record the exact source
  revisions, compiler, flags, hardware, and external-code versions for every
  benchmark.
- Keep generated plots and large benchmark output out of the source tree.
  GitHub Actions regenerates the gallery. Large or license-sensitive oracle
  data belongs in the separate benchmark-data repository, with a small
  manifest and immutable revision recorded here.
- Prefer pure, elemental, allocatable-safe Fortran. Avoid hidden global state,
  pointer ownership, aliasing, and allocation inside element kernels.
- Profile before optimizing. Use Fortran and existing `fortnum` or
  `fortsparse` facilities first. Do not add a new dependency merely to make a
  benchmark table look better.
- Treat performance as an acceptance criterion. Generated kernels should be
  allocation-free in their hot path, vectorizable, batchable, and measurable
  with a roofline or operation-count explanation. FortFEM targets competitive
  or better performance for each named workload. A claim such as "faster than
  the competition" is valid only for a declared problem, compiler, hardware,
  thread count, and reference implementation.
- Keep the implementation serial-first and distributed-ready. Global IDs,
  ownership, ghost data, deterministic reductions, local and global index
  maps, and communicator operations are part of the data model before a full
  MPI backend is added. Residual and element kernels remain independent of the
  communicator. This preserves a short single-node feedback loop while
  avoiding a later rewrite of mesh, trace, and sparse interfaces.
- Add a deterministic property-testing path to `fo` and the test helpers.
  Every randomized case records its seed, generated topology, and tolerance.
  Seeded random tests provide fast coverage of orientations, geometry
  degeneracies, mode sets, constitutive tensors, and solver dimensions. They
  complement analytical and manufactured oracles and never replace them.

## 2. Current baseline

The following capabilities are already on FortFEM `main` or in the verified
documentation baseline. The list is intentionally conservative.

| Area | Current state | Next gate |
| --- | --- | --- |
| Scalar FEM | P1/P2/Q1 and arbitrary-order triangular scalar paths, Poisson and diffusion forms, boundary conditions, plotting | General field and coefficient callbacks in the symbolic form compiler |
| FEEC | Oriented triangular and tetrahedral H1, H(curl), H(div), and DG families, Piola maps, commuting tests, sparse assembly, mixed RT-DG Poisson | General multi-field block composition and arbitrary multipatch assembly |
| IGA | Nonuniform B-splines, rational maps, two- and three-dimensional de Rham incidence complexes, cylindrical and toroidal Fourier blocks, initial JOREK magnetic-flux residual/JVP | General patch graphs, trimming, enrichment, and the remaining coupled JOREK variables |
| Special functions | FortNum quadrature, ordinary and associated Legendre P/Q, orthonormal complex spherical harmonics with angular derivatives, Hobson-normalized toroidal P/Q branches, Bessel/Hankel paths, and a FortFEM Fourier mode registry with phase/radial derivative contracts; the pinned spherical and toroidal APIs are re-exported through `fortfem_api` and checked by integration oracles | Stable high-order and near-cut half-integer continuation, spherical-harmonic products, and cross-geometry special-function oracles |
| Sparse algebra | FortSparse CSC assembly, retained factors, real and complex solves, sparse products, tree--cotree CSC direct reductions with fixed-map JVP/VJP, and CG, PCG, GMRES, and BiCGSTAB converged-state derivative contracts; dense and standalone sparse IC(0)/ILU(0) factor/apply paths plus deterministic sparse fixed-factor ILUT, memory-scalable row-oriented ILUT, and controlled ICHOL paths are public | Scalable row-oriented ICHOL construction, measured scaling, flexible Krylov products, and block solver derivatives |
| Open boundaries | Planar, circular, and spherical scalar Helmholtz DtN paths, scalar BEM, Maxwell trace and PML components | General curved Maxwell DtN, toroidal exterior maps, and robust FEM/BEM/DtN comparison fixtures |
| PML | Scalar and curl-curl Cartesian complex-stretching tensors with slab, triangular, and tetrahedral examples | Automated curved-object layers, reflection/error metrics, and derivative coverage for all geometry parameters |
| Differentiation | Analytical FortSym paths, selected Enzyme checks, sparse matrix products, converged CG/PCG/GMRES/BiCGSTAB solves, toroidal coordinate and DtN products | Complete operator inventory, JVP/VJP parity for all public operators, and shape derivatives |
| Parallel readiness | Serial local kernels and deterministic focused tests | Owned/ghost mesh and field data, partition-independent numbering, halo exchange, distributed assembly, checkpointing, and MPI-enabled solver backends |
| Examples | Generated documentation pages for Poisson, Maxwell, Helmholtz, BEM, IGA, torus, PML, and solver examples | Ordered gallery beginning with simple Poisson and adding 1D, 2D, and 3D application-oriented toy models with manufactured or external oracle data |

The aggregate full test suite can expose resource-sensitive failures when many
independent numerical targets run in parallel. Focused repetitions are the
first diagnostic and remain the merge gate for a local increment. CI remains
the final cross-platform gate. The current fast suite has 430 targets; 25
deliberately long regression targets are classified as `_slow` and retain a
separate explicit run path. A result formatter truncation is never treated as
a pass: batch or named runs are used when a report exceeds the display limit.

## 3. Stable contracts before physics applications

### 3.1 Field, space, and residual contracts

Every public PDE object must declare:

- domain dimension, coordinate chart, orientation, and physical units;
- trial, test, trace, and multiplier spaces;
- scalar, real-vector, complex-vector, or tensor value type;
- primal residual and active constraints;
- boundary and internal-interface terms;
- nullspace or gauge treatment;
- assembled and matrix-free actions;
- JVP and VJP actions, including parameter and geometry derivatives;
- independent manufactured forcing and an analytical or cross-code oracle.

The standard interface is conceptually:

```text
residual(state, parameter) -> R
jvp(state, parameter, state_dot, parameter_dot) -> R_dot
vjp(state, parameter, R_bar) -> state_bar, parameter_bar
solve(residual, constraints) -> converged_state
implicit_jvp(converged_state, ...) -> state_dot
implicit_vjp(converged_state, ...) -> state_bar, parameter_bar
```

The solve derivatives assume a converged state and hold iteration and stopping
branches inactive. This contract must be visible in each solver's documentation.
For complex fields, the VJP convention and conjugation must be stated next to
the API and checked by a dot-product identity.

### 3.2 Topology and geometry contracts

Meshes, spline patches, Fourier charts, and internal manifolds share a common
orientation model. Geometry objects provide values, Jacobians, inverse maps,
surface normals, measures, and derivatives with respect to control parameters.
The same physical problem can be represented by a fitted interface or by an
unfitted level set without changing the interface-law API.

### 3.3 Foundation gap register

The current algebraic slices are reusable building blocks. The following
contracts still have to be composed before FortFEM can serve as the basis for
an arbitrary-topology three-dimensional MHD or edge application.

| Foundation | Contract still missing | Required independent oracle |
| --- | --- | --- |
| Topological complex | A region and cell-complex graph with periodic identifications, orientations, homology, cohomology, harmonic representatives, cuts, and gauges; metric-harmonic one-forms now have a fixed-topology cycle-period normalization map with JVP/VJP | Chain-complex identities, Euler characteristic, cycle and flux integrals, normalized periods, and nullspace dimension on slab, cylinder, sphere, and torus cells |
| Sheet-current interface | Neutral open/closed internal-manifold graph, integrated-current junction ledger, fixed-topology loop-current constraints, differentiable geometry-to-edge-flux contraction, topology-only edge-flux balance, normal-traction jump, and an independent test/trial tangential surface-current trace residual are public; constitutive pressure laws and flux/helicity constraints remain | Ampere jump, surface-current conservation, loop current, pressure jump, and regularized-layer limits |
| Cut FEEC spaces | The scalar shifted-Heaviside activation, a 3D vector-enrichment curl/divergence product-rule diagnostic, the physical vector/metric support Gram contraction, batched 2D/3D covariant/contravariant Piola-enrichment composition, and a rectangular commuting-projection audit with value/JVP/VJP actions are public; Piola-aware vector-compatible XFEM/XIGA and DG spaces that preserve or explicitly report the de Rham sequence across cuts remain | Actual enriched-space construction, curl-gradient and divergence-curl identities on every generated space, and fitted versus unfitted convergence |
| Coupled field residuals | Generic composable blocks for vector fields, tensor constitutive laws, interfaces, constraints, and boundary operators; the neutral tensor volume-work contraction is public. Plasma state assembly remains in an external client | FortSym manufactured residuals, block-Jacobian products, energy or power balance, and cross-formulation parity |
| Equilibrium interchange | The neutral external-adapter schema for mapped coordinates, physical samples, named scalar/vector/tensor coefficients and scalar profiles, segmented boundaries, provenance, units, and normalization is public and validated. GEQDSK and COCOS parsing remain outside FortFEM | Analytic manufactured data now covers toroidal mapped samples, named fields, segmented boundaries, vector/tensor component offsets, copy semantics, and rejection; license-safe CHEASE and FreeGS outputs sampled on a common physical grid remain |
| Fourier and toroidal modes | The fixed-topology mode registry now provides field-period phase, normalization, conjugate packing, retained-triad lookup, radial regularity, complex coordinate derivatives, and a generic three-factor vector/tensor Fourier product with JVP/VJP; FortNum now exposes independently tested ordinary Legendre Q values/derivatives, orthonormal complex spherical harmonics with angular derivatives, and Hobson-normalized half-integer toroidal P/Q values/derivatives; model-specific mode operators and branch data remain | Stable high-order recurrences/continuation, symmetry and de-aliasing checks, and independent mode-by-mode energy |
| Edge and SOL equations | Equation-as-data fields, generic coefficient and boundary callbacks, conservative sources, FCI events, and target ledgers. Species and closures remain client-owned | Manufactured source terms, mass and energy balances, terminal flux tallies, and a reproducible FCI map |
| Mixed waves and elasticity | A common compatible port-Hamiltonian state for pressure, velocity, displacement, momentum, and tensor stress, including boundary power ports, plus a separate dissipative Cayley block | Discrete energy, symplectic-form or passivity tests, dispersion, reversibility, and mixed versus second-order parity |
| Open boundaries | Curved vector FEM/BEM/DtN/PML coupling on toroidal external surfaces with larger-domain controls | Reciprocity, passivity, far-field, reflection, and interior-field agreement across all four paths |
| Verification and delivery | Seeded random tests, external-code adapters, provenance manifests, mesh-completeness checks, and Pages health checks | Repeated seeds, license and revision records, independent samplers, HTTP link checks, and FortPlot image regression |

Closing a row requires a public API, a focused test, an independent oracle, a
FortSym or provenance entry, and a gallery fixture when the row affects a
visible numerical result. A current primitive does not imply that its coupled
application residual is complete.

## 4. Structure-preserving numerics

Structure preservation is measured, tested, and shown in the gallery.

### 4.1 Spatial structure

- Maintain the discrete de Rham chain
  \(H^1\xrightarrow{\nabla}H(\mathrm{curl})
  \xrightarrow{\nabla\times}H(\mathrm{div})
  \xrightarrow{\nabla\cdot}L^2\) whenever the formulation requires it.
- Generate and test incidence matrices independently of metric matrices.
- Check \(\nabla\times\nabla=0\) and
  \(\nabla\cdot\nabla\times=0\) algebraically, including multipatch and
  periodic identifications.
- Preserve normal magnetic flux in H(div), tangential electric or magnetic
  traces in H(curl), and the declared jump when a surface current is present.
- Use exact or compatible divergence control. A cleaning term is never a
  substitute for a space that has the required topology.
- Keep skew, symmetric, positive, and constraint blocks identifiable in block
  systems. These properties drive solver choice and preconditioning.

### 4.2 Time structure

The time integrator is selected from the continuous structure.

- Ideal or nondissipative Hamiltonian components use variational, symplectic,
  Poisson, or multisymplectic updates. Implicit midpoint, Cayley, splitting,
  and composition methods are appropriate building blocks.
- Resistive, viscous, and thermal terms are dissipative. Their discrete energy
  law must be monotone or have the declared balance. They are not described as
  symplectic merely because they are solved implicitly.
- Coupled ideal plus dissipative systems use a structure-aware split, such as
  a symmetric composition of reversible and dissipative substeps, with order
  and invariant tests.
- Time-dependent examples report energy, magnetic helicity, cross-helicity,
  mass, divergence, charge, and constraint defects as appropriate.
- Generic Runge--Kutta stepping is allowed for a deliberately non-geometric
  transport example, but it is not the default for a Hamiltonian or
  structure-preserving path.
- Spacetime FEM is reserved for a later phase. The present roadmap requires a
  sound space-time contract without implementing spacetime FEM now.

### 4.3 Structure tests

Each time or operator increment adds at least one of:

1. an exact algebraic identity;
2. a conservation or dissipation balance;
3. a reversibility or symplectic-form test;
4. a manufactured solution and convergence rate;
5. an independent energy or adjoint identity.

### 4.4 Mixed first-order waves and mechanics

Wave equations are not forced into a second-order pressure-only form when a
first-order mixed form exposes more structure. The reusable state contract is
of the form

\[
\partial_t q + C^T v = 0,\qquad
M_v\partial_t v - Cq = f,
\]

where `q` may be pressure, displacement, electric flux, or another potential,
`v` is its flux or velocity companion, and `C` is a compatible discrete
gradient, divergence, curl, or elasticity-complex map. The semidiscrete
Hamiltonian or port-Hamiltonian form must expose its skew interconnection and
positive mass or constitutive blocks.

The roadmap therefore includes:

- mixed first-order acoustics with pressure and particle velocity;
- mixed elastodynamics with displacement, velocity, and stress or momentum;
- Maxwell and elastic wave systems with compatible electric, magnetic,
  traction, and displacement traces;
- acoustic-elastic and electromagnetic-elastic interface coupling;
- symplectic Euler, implicit midpoint, Cayley, variational, and partitioned
  symplectic updates for the nondissipative part;
- dissipative viscosity, damping, conductivity, and absorbing layers as
  metriplectic or energy-decaying substeps;
- exact semidiscrete energy conservation or a declared discrete power balance,
  plus dispersion and phase-error diagnostics.

This contract applies to acoustics and to the general wave family. A
second-order reduction remains available as a compatibility path, but it must
be demonstrably equivalent to the mixed first-order residual and must not hide
the conserved pairing needed by a geometric time integrator.

### 4.5 Tensor-valued pressure, stress, and constitutive laws

Pressure is a scalar only for an isotropic closure. FortFEM must support a
tensor-valued pressure or stress contribution \(\boldsymbol\Pi(x,u,p)\), with

\[
\boldsymbol\sigma=-\boldsymbol\Pi+\boldsymbol\tau,
\qquad
\mathbf f_{\mathrm{int}}=\nabla\!\cdot\boldsymbol\sigma,
\]

where the symmetric, nonsymmetric, deviatoric, gyrotropic, or field-aligned
parts are declared by the constitutive model. The form compiler must preserve
tensor symmetry and produce stress, traction, JVP, VJP, and shape derivatives.
Required capabilities are:

- anisotropic elasticity and nearly incompressible mixed elasticity;
- Hellinger--Reissner stress-displacement and weak-symmetry formulations;
- tensor constitutive coefficients with field-aligned, gyrotropic, or
  caller-defined anisotropy;
- symmetric and nonsymmetric stress with explicit angular-momentum balance;
- traction and normal-flux traces on fitted, cut, and BEM interfaces;
- constitutive tensors with spatial, parameter, and magnetic-field dependence.

For a field-aligned constitutive oracle, the CGL-shaped tensor is a first test
case,

\[
\mathbf P=p_\perp I+(p_\parallel-p_\perp)\,\mathbf b\mathbf b^T,
\qquad \mathbf b=\mathbf B/|\mathbf B|,
\]

with optional caller-supplied gyrotropic terms. The parallel and perpendicular
coefficients, their force divergence, and their work pairing are separate
residual terms. FortFEM does not provide a plasma closure or a Braginskii
model.

The first constitutive slice is now public on `main`: FortSym generates the
six independent symmetric CGL components and their JVP/VJP, while the
`fortfem_cgl_pressure_tensor` wrapper validates the unit magnetic direction,
packs a full symmetric tensor, and combines full-matrix off-diagonal
cotangents. The generated tensor now feeds a compositional traction
`t=P n` with value/JVP/VJP, including the `P^T` normal cotangent; independent
tests cover the closed-form traction oracle, central differences, adjoint
identities, and invalid directions. Generic force-volume assembly, caller-
supplied correction tensors, and field-aligned assembly remain separate active
work. These blocks are not a claim that an anisotropic MHD application is
implemented.

The elasticity complex is treated as a structure-preserving extension of the
de Rham complex. The mixed weak-symmetry construction of
[Arnold, Falk, and Winther](https://arxiv.org/abs/math/0701506) is the initial
literature oracle. The tangential-displacement, normal-normal-stress (TDNNS)
work from the [JKU Linz group](https://www.numa.uni-linz.ac.at/Talks/abstract/154/)
and [Pechstein and Schöberl](https://doi.org/10.1002/nme.3319) is the
anisotropic and nearly incompressible implementation reference. The first
FortFEM examples are a manufactured anisotropic elasticity patch, a nearly
incompressible block, and a tensor-pressure wave that compares
displacement-only and mixed stress formulations.

## 5. Geometry, topology, and interfaces

### 5.1 Embedded geometry

Implement first-class geometry for:

- level sets \(\phi(x,p)=0\), with gradients, Hessians, Boolean operations,
  periodicity, and parameter derivatives;
- cut triangles, tetrahedra, spline cells, and Fourier cells;
- adaptive, moment-fitting, recursive, and exact polynomial cut-cell
  quadrature;
- closed or open curves and surfaces, intersections, junctions, nested
  interfaces, separatrices, X-points, and rational magnetic surfaces;
- fitted refinement, anisotropic refinement, h, p, and hp refinement,
  hierarchical and truncated hierarchical B-splines;
- shape derivatives with explicit diagnostics when a cut crosses a node or a
  topology changes.

An interface is an object rather than a boundary marker:

```text
interface
    geometry
    plus_domain, minus_domain, orientation
    volume_spaces, trace_spaces, skeleton_spaces
    jump_laws, constraints, surface_unknowns
    quadrature, shape_derivatives
```

### 5.2 Broken and trace spaces

Provide broken versions of H1, H(curl), H(div), and L2. Their trace API must
support scalar, normal, tangential, rotated-tangential, plus, minus, average,
and jump operators with one orientation convention. Skeleton spaces support
scalar traces, normal-flux traces, tangential-field traces, surface current,
pressure-jump multipliers, and interface displacement.

### 5.3 Interface laws and distributional sources

The weak-form compiler must represent both volume and surface terms:

\[
\int_\Omega f v\,dV + \int_\Gamma g\,\operatorname{tr}(v)\,dS.
\]

For electromagnetics, the surface-current law is

\[
\mathbf K=\mathbf n\times[\![\mathbf H]\!],
\qquad
\mathbf J_{\mathrm{sheet}}=\mathbf K\,\delta_\Gamma.
\]

The implementation should support duplicated degrees of freedom, explicit
surface-current unknowns, Lagrange multipliers, mortar coupling, Nitsche
coupling, hybridization, cohesive-like jump laws, and regularized-layer
comparisons. A delta source is not approximated numerically when an explicit
surface integral is available.

### 5.4 Arbitrary topology and sheet-current composition

The interface contract becomes an application-ready foundation only after the
geometry and spaces are connected through a topological region graph. The
graph must represent volume regions, oriented facets, boundary components,
periodic identifications, open edges, closed loops, and junctions. It assigns
stable IDs to cells, facets, traces, cycles, and cohomology representatives so
that a mesh rebuild does not silently change a gauge or a Fourier mode.

The graph owns the following structure-preserving data:

- boundary, coboundary, and metric matrices with separate incidence and Hodge
  factors;
- harmonic one- and two-form bases for multiply connected domains and their
  flux normalization;
- gauge and nullspace constraints for vector potentials, scalar potentials,
  and pressure-like multipliers;
- plus/minus trace ownership for every interface facet, including a single
  owner for cut or overlapping FCI terminal pieces;
- topology-event records when an interface intersects a vertex, an edge, a
  periodic seam, or another interface.

The first quotient-complex slice is now public on `main`: signed vertex,
edge, face, and volume representative maps are composed with integer chain
boundaries, and the output is checked with exact boundary-of-boundary
identities. The interval-to-circle test supplies the independent Betti-number
oracle. This is still metadata only; quotient coordinates, metric Hodge
operators, gauges, and application-owned periodic laws remain later layers.

The cell-complex cycle-space slice is also public: scale-aware real kernels of
the oriented edge boundary and its face-boundary transpose are available for
subsequent cuts, flux normalization, and gauge construction. It intentionally
does not label those cycles or cocycles as harmonic forms or physical fluxes
until metric Hodge data and application constraints are supplied.

The metric-harmonic one-form slice is public as well. Given a positive edge
metric it returns the closed and co-closed kernel. The fixed-topology
`normalize_harmonic_one_forms` slice now maps that basis to caller-selected
cycle periods or flux units with dense-solve JVP/VJP actions, while leaving
cycle labels, physical units, and gauge ownership to the caller. The
fixed-topology tree--cotree gauge slice supplies a spanning-forest split,
cotree restriction/prolongation, direct dense-system reduction, and fixed-map
JVP/VJP actions. Its selector is frozen for differentiation; graph rebuilds
are topology events.
For high-order FEEC and IGA, the tree is built on the lowest-order
mesh/control graph and higher-order moments are retained by an external DOF
map.

The direct-gauge design follows [Manges and Cendes
(1995)](https://doi.org/10.1109/20.376275), [Bíró, Preis, and Richter
(1996)](https://ieeexplore.ieee.org/document/558631), and the high-order
tree-gauge analysis in [the tree-gauge review](https://doi.org/10.3390/j5010004).
The mortared IGA extension is tracked against [Kapidani et al.
(2022)](https://doi.org/10.1016/j.cma.2022.114654) and its
[preprint](https://arxiv.org/abs/2110.15860). These are literature and oracle
references, not imported source code.

An electromagnetic sheet is a coupled unknown, not only a post-processed
delta source. The implementation must support equivalent formulations and
their conversion:

1. a jump in tangential \(\mathbf H\);
2. an explicit tangential \(\mathbf K\) trace space;
3. a distributional volume current \(\mathbf K\delta_\Gamma\);
4. a resolved resistive layer approaching the sheet limit.

Each formulation carries the same orientation, units, edge balance, and
current-loop normalization. On a closed interface the surface divergence and
global current constraints are explicit. On an open interface the edge flux is
owned by a boundary law. For MHD interfaces, the sheet law is coupled to normal
magnetic-flux continuity, tangential field jumps, total-pressure balance,
field-line or loop-flux constraints, and optional helicity constraints.

The cut or enriched H(curl) and H(div) spaces must expose whether an
enrichment preserves the discrete de Rham sequence. A scalar Heaviside basis
cannot be promoted to a compatible vector space by naming convention alone.
The compiler and tests therefore record the commuting diagram, the required
trace regularity, and the exact identity that is intentionally relaxed when a
physical jump is present.

## 6. Approximation families and discontinuities

FortFEM supports two complementary representations of every important
discontinuity.

| Mechanism | Use case | Required ingredients |
| --- | --- | --- |
| Fitted interface | Known rational surface, material boundary, shielding sheet, multi-region equilibrium | Duplicated traces, interface elements, mortar or Nitsche, multipliers, surface current, independent constraints |
| Unfitted XFEM/GFEM | Moving or emerging interface, island separatrix, geometry optimization | Level sets, cut quadrature, partition of unity, shifted enrichment, blending correction, conditioning diagnostics |
| XIGA, immersed, trimmed, CutIGA | Spline geometry with interfaces not aligned to knot lines | NURBS/B-spline enrichments, connected-component activation, trimmed-support stabilization, shape derivatives |
| DG | Shocks, steep layers, nonmatching cells, discontinuous material or pressure | SIPG, NIPG, IIPG, upwind, entropy-stable fluxes, broken FEEC spaces |
| HDG | Local elimination and interface unknowns | Skeleton traces, static condensation, compatible H(curl)/H(div) fluxes |

### 6.1 XFEM and XIGA library

The initial enrichment library includes Heaviside jumps, signed-distance and
absolute-value kinks, logarithmic and algebraic singularities, crack-tip and
corner modes, user-defined analytical enrichments, vector and tensor
enrichments, helical Fourier factors, and resonant radial forms near
\(m-nq=0\). Shifted enrichments
\(N_i(x)(\Psi(x)-\Psi(x_i))\) preserve interpolation where possible.

The corrected blending-element construction from
[Fries (2008)](https://doi.org/10.1002/nme.2259) is required. Enrichment
activation must include support-size thresholds, rank-revealing checks,
orthogonalization or pruning, diagonal scaling, local condensation, and
condition-number diagnostics.

XIGA must support H1, H(curl), H(div), L2, multipatch, and trimmed supports.
Connected components of a spline support receive independent branches when
the physical support is disconnected. Trimmed B-splines must also support the
small-support stabilization described in
[stable trimmed IGA work](https://www.sciencedirect.com/science/article/abs/pii/S0045782516308222).
Physical jump enrichment and trimmed-space stabilization remain separate
composable mechanisms.

### 6.2 DG and hybridization

The planned DG layers are scalar elliptic interior penalty, hyperbolic central
and upwind fluxes, Rusanov/HLL-style fluxes where appropriate, broken vector
forms with tangential and normal penalties, and HDG with global skeleton
unknowns. Mixed CG-DG and XFEM-enriched DG are valid combinations. Every flux
declares whether it preserves conservation, entropy, energy, or only
consistency.

## 7. IGA, Fourier-FEM, and toroidal geometry

### 7.1 IGA and patch graphs

Move from the current two-patch checks to an arbitrary patch graph with
periodic identifications, toroidal closure, orientation reversals,
extraordinary junctions, nonmatching knot vectors, mortar coupling, and local
refinement. Continuity is field-specific: smooth regions may use
\(C^{p-1}\), interfaces may use \(C^0\), repeated knots may encode derivative
jumps, and traces may be fully discontinuous.

Trimmed and cut IGA are needed for coils, conductors, vacuum vessels, islands,
moving rational surfaces, and non-parametric interfaces.

### 7.2 Fourier-FEM contract

Fourier-FEM combines a Fourier or helical expansion in periodic angles with
finite elements or splines in the remaining coordinates. The canonical
representation is documented as

\[
u(s,\theta,\zeta)=
\sum_{m,n}u_{mn}(s)\,e^{i(m\theta-nN_{\mathrm{fp}}\zeta)}.
\]

The implementation must define phase, field-period, real/conjugate packing,
mode ordering, normalization, and de-aliasing. It must support:

- radial FE or B-spline spaces with regular-axis and edge constraints;
- 2D poloidal FE plus toroidal Fourier modes;
- mode-diagonal linear operators;
- nonlinear convolution with explicit retained-triad policy;
- complex block systems and their real equivalent;
- matrix-free mode actions and adjoints;
- parameter derivatives of mode numbers, geometry, and Fourier coefficients;
- Fourier-mode plots and energy per mode.

This is the common numerical ingredient for CHEASE-like axisymmetric models,
GLISS and CAS3D-like linear calculations, GPEC/MARS-like perturbations, and
JOREK-like reduced-MHD operators. It is not a claim that all those codes use
the same formulation.

### 7.3 Special-function and geometry scope

FortNum now provides the real ordinary Legendre second-kind contract
`legendre_q`/`legendre_q_derivative` for x > 1, in addition to Ferrers
associated `P_l^m`, and its independently tested Hobson-normalized
half-integer toroidal `P/Q` API. It also provides standard orthonormal complex
`spherical_harmonic` values and analytical theta/phi derivatives on the closed
polar interval (with derivative evaluation intentionally undefined at its
poles). The [pinned FortNum API
revision](https://github.com/lazy-fortran/fortnum/blob/b5afac7/docs/api.md)
records the domains, normalization, derivative convention, DLMF provenance,
and the current moderate-degree limitation. The remaining foundation is
stable high-order and near-cut continuation, spherical harmonics and products,
and cross-geometry special-function oracles. Analytical solutions and DtN
maps are tested on slabs, cylinders, spheres, and exact toroidal surfaces. The
FortFEM now re-exports the pinned toroidal P/Q values and derivatives through
`fortfem_toroidal_harmonics`, with an independent branch/derivative/domain
oracle. The torus API must state whether a function is a toroidal harmonic, a
Fourier mode in toroidal coordinates, or a numerical continuation of a special
function. These are not interchangeable names.

## 8. PDE and boundary-operator ingredients

The same residual and derivative contracts cover scalar and vector-valued
equations.

### 8.1 Core equations

- Poisson, diffusion-reaction, and generic axisymmetric elliptic/Fourier
  operators that an external Grad-Shafranov client can instantiate;
- scalar Helmholtz with FEM, BEM, DtN, and PML boundaries;
- Ampere, curl-curl, eddy-current, Maxwell, and anisotropic H(curl) forms;
- H(div) flux, mixed Poisson, magnetic induction, and divergence constraints;
- generic linearized field, constitutive, interface, and wall-response blocks
  that an external ideal or resistive MHD client can compose;
- generic skew brackets, energy functionals, and compatible time operators for
  client-supplied state variables.

Pressure-like and stress-like quantities are represented as scalar, vector, or
tensor fields supplied by the caller. A tensor-valued coefficient is not
silently projected to its trace. Its symmetry, divergence, traction, and work
pairing are tested with manufactured data. A CGL-shaped tensor is one generated
constitutive oracle, not a plasma-closure implementation.

### 8.2 Strongly anisotropic and field-aligned operators

Plasma transport and wave models often contain coefficients whose parallel
and perpendicular scales differ by many orders of magnitude. FortFEM must
represent this explicitly rather than hiding it in an isotropic scalar
diffusion coefficient. For a unit magnetic direction
\(\mathbf b=\mathbf B/|\mathbf B|\), the basic decomposition is

\[
\mathbf K
 = k_\parallel\,\mathbf b\mathbf b^T
 + k_\perp\,(I-\mathbf b\mathbf b^T)
 + \mathbf K_{\mathrm{Hall}},
\]

with optional nonsymmetric Hall or gyroviscous terms. The same idea applies to
viscosity, thermal conduction, resistivity, pressure, elasticity, and wave
impedance. Required work includes:

- covariant and contravariant tensor pullbacks on curvilinear meshes;
- field-line and flux-coordinate-independent (FCI) parallel derivatives;
- perpendicular plane operators and parallel line operators with separate
  quadrature and tolerances;
- robust assembly for \(k_\parallel/k_\perp\) ratios spanning many decades;
- anisotropy-aware preconditioners, field-split Krylov methods, and fGMRES;
- exact energy or dissipation tests for symmetric and skew tensor parts;
- B-field, tensor, and geometry JVP/VJP actions;
- limiting tests as \(k_\parallel/k_\perp\to 1\), \(k_\perp\to0\), and
  \(k_\parallel\to\infty\).

The parallel operator may be represented by a compatible H(curl)/H(div)
complex, a Fourier derivative, or an FCI field-line map. These choices share a
residual and oracle interface but are not assumed algebraically identical.
The generic pointwise field-aligned flux

\[
  \mathbf F=k_\perp\mathbf g+(k_\parallel-k_\perp)\mathbf b
  (\mathbf b\cdot\mathbf g)
\]

is now public with FortSym-generated value/JVP/VJP products and independent
unit-direction, finite-difference, and dot-product tests. It is the common
constitutive block for future anisotropic diffusion, conduction, resistivity,
and wave assemblies.
The first PARALLAX-aligned algebraic slice is now on `main`: a dependency-light
RK4 field-line tracer provides geometry endpoints; mapped upper and lower plane
interpolation matrices assemble a sparse staggered gradient; the support
divergence is constructed as its negative volume-weighted adjoint; and a
matrix-free (P K_\parallel Q) diffusion action is public. The matrix-free
gradient action exposes FortSym-generated value, JVP, and VJP kernels.
Independent tests cover a constant/exponential tracing oracle, a linear-field
map oracle, constants, an explicit flux-balance vector, a central-difference
JVP, the weighted adjoint identity, and the weighted negative-energy identity.
The `fci_parallel_diffusion` gallery example now runs the same matrix-free
action on a manufactured open-line cosine profile and publishes 1D FortPlot
profiles plus CSV values for the mass-rate and dissipation oracles.
The field-only VJP of this diffusion action is also public and is checked by
an independent dot-product identity. The complete fixed-topology diffusion
JVP/VJP now covers interpolation maps, line lengths, parallel coefficients,
canonical and staggered volumes, and the field through pinned FortSym local
contribution kernels; a central-difference and full real dot-product oracle
guard this contract. The public 1D linear interpolation-map builder now checks
partition of unity, affine reproduction, fixed-topology JVP/VJP dot products,
and Cartesian bilinear affine reproduction. A generated quadratic Lagrange
map now accepts explicit three-node stencils and reproduces quadratic fields
on nonuniform slices, with generated fixed-stencil JVP/VJP dot-product and
finite-difference oracles. Higher-order interpolation Jacobians, higher-order
curved support-volume measures, and anisotropy-aware preconditioning remain
active work. A generated fixed-topology quadrilateral area map now supplies positive
unstructured plane-cell measures with value/JVP/VJP actions, independent
shoelace/finite-difference/adjoint oracles, and a FortPlot gallery fixture. A
generic fixed-topology polygon map now accepts boundary-ordered cells with
any number of vertices at least three, with a generated per-edge contribution,
shared-vertex JVP/VJP accumulation, independent pentagon oracle, and gallery
fixture. Its quadratic Bezier-edge extension now supplies arbitrary curved
polygon value/JVP/VJP measures with independent Gauss--Green and gallery
oracles. A batched Cartesian bilinear adapter now
turns traced forward/backward endpoint arrays into the per-segment FCI map
tensors used by the support operator, with source-grid accumulation in its
fixed-topology JVP/VJP and independent finite-difference/dot-product oracles.
An unstructured fixed-cell triangle adapter now supplies barycentric maps and
geometry/target JVP/VJP products with affine and dot-product oracles; higher
order triangle and moving-cell connectivity remain active. Its batched
forward/backward endpoint adapter now emits the per-segment tensors consumed
by the support operator and accumulates shared-vertex cotangents.
The positive diagonal of `-W_c^{-1}Q^TW_sK_\parallel Q` is now public as a
FortSym-generated per-stencil Jacobi baseline, with an explicit Q-squared
oracle and a validated matrix-free diagonal apply; plane multigrid and
field-split preconditioners remain active.
The split action now also accepts one validated square CSC perpendicular block
per canonical plane and combines it with the matrix-free FCI parallel term,
with a negative-energy anisotropic oracle. Its combined positive diagonal is
also public: the generated support diagonal is combined with the negative
diagonal of each oriented plane block and checked against an independent
oracle, giving a reproducible scalar preconditioner baseline while plane
multigrid and stronger field splitting remain active work. A convenience
matrix-free anisotropic Jacobi apply now uses the same combined diagonal;
iterative callers may cache the diagonal for repeated solves. A two-level plane
V(1,1) cycle now supplies CSC restriction/prolongation and a replaceable direct
coarse solve; deeper hierarchies and retained coarse factors remain active.

### 8.3 Equilibrium, perturbation, and edge equation data

The named application families share interchangeable data and residual
contracts. FortFEM supplies those contracts and small fixtures. It does not
ship their production solvers.

#### 2D equilibrium data

An external equilibrium adapter may provide a poloidal mesh or spline,
coefficient fields, axis or boundary metadata, coil or wall traces, units, and
normalization. FortFEM consumes this neutral data through a documented schema.
It does not parse or write GEQDSK, interpret COCOS conventions, or implement
equilibrium profiles. A generated manufactured field is the FortFEM baseline.
CHEASE and FreeGS can provide optional numerical oracle data through an
external, license-safe adapter sampled on the same physical points.

#### Linear ideal and resistive response data

The linear schema must distinguish an eigenproblem from a forced response and
must carry complex frequency, Fourier mode set, equilibrium coefficients,
vacuum and wall regions, interface traces, singular-layer matching data,
normalization, and response-matrix conventions. The residual exposes the
energy-principle, inertia, resistive, wall, and vacuum blocks separately. A
response matrix has reciprocity, passivity, and normalization tests. GPEC,
MARS-F, GLISS, and STARWALL remain versioned oracle or interchange targets.

#### Nonlinear MHD composition

The reusable nonlinear state schema must admit generic scalar, vector, and
tensor fields, constitutive coefficients, interfaces, and constraint
multipliers. An external application can map density, velocity, magnetic field,
pressure, current, or other variables into that schema. FortFEM provides the
residual composition, energy or power ledger, surface-current coupling, and
derivative actions. Plasma closures, profile laws, wall physics, and state
selection remain outside the library. Every optional numerical block has an
independent manufactured case before an application combines it.

#### Edge and SOL equation data

An edge application can map its fields and coefficients into generic equation,
source, conservative-flux, boundary-event, and material-ledger callbacks. The
core owns FCI geometry, anisotropic operators, trace spaces, and residual/JVP/
VJP composition. Species, Braginskii coefficients, sheath or Bohm laws,
recycling, radiation, neutral or impurity sources, and target material data
remain client-owned. The callback ABI records units and signs so that a
terminal ledger can be compared with the volume balance. This is the
foundation needed by GRILLIX, BOUT++, and PARALLAX-style clients without
copying or implementing their models.

### 8.4 FEM/BEM, DtN, and PML

The open-boundary layer must support scalar and vector equations in 2D and 3D,
including toroidal external surfaces:

- Laplace and Helmholtz single, double, hypersingular, Calderon, CFIE, and
  symmetric transmission blocks;
- Maxwell EFIE, MFIE, CFIE, RWG, BC/RBC, and higher-order trace pairings;
- planar, circular, spherical, and exact-curved-torus DtN and NtD maps;
- FEM/BEM and FEM/DtN coupling with consistent work-conjugate traces;
- Cartesian and curvilinear PML for scalar Helmholtz and curl-curl Maxwell;
- analytical Mie, cylindrical, toroidal, and manufactured field oracles;
- far-field, flux, reciprocity, passivity, and reflection diagnostics.

Every boundary path has an explicit larger-domain comparison: move the
artificial boundary away, compare FEM/BEM, DtN, and PML fields on a common
interior region, and report the expected nonzero but convergent difference.

## 9. Solver, constraints, and differentiation roadmap

### 9.1 Sparse and nonlinear solvers

Complete the common solver product layer for direct solves, CG, PCG, GMRES,
BiCGSTAB, flexible Krylov, block systems, Schur complements, static
condensation, matrix-free actions, field splits, and retained factorizations.
Each converged-state solver has an implicit JVP/VJP. Iteration counts and
stopping decisions are inactive in that derivative contract.

For curl--curl and other gauge-singular systems, direct solves must use an
explicit tree--cotree or equivalent compatible nullspace reduction before
factorization. The selector is topology metadata and is not differentiated
through. Iterative paths must expose compatible nullspace projection and
measured ILU/ILUT, incomplete Cholesky (IC(0)), or incomplete Cholesky with
level/fill control (ICHOL) options where the matrix symmetry and definiteness
permit them. Complex or nonsymmetric blocks use an explicitly declared
replacement (for example ILU rather than ICHOL), with factor breakdown and
loss of positive definiteness reported rather than hidden.

The dense IC(0) factor/apply primitive and the PCG `ichol`/`ic0` option are now
public and independently tested. Standalone sparse IC(0) and ILU(0) paths now
consume FortSparse CSC directly, preserve the input lower/upper patterns
without fill, and expose fixed-factor right-hand-side JVP/VJP actions,
including the transpose solve for ILU. PCG integration, ILUT/fill-controlled
ICHOL, and scaling benchmarks remain active solver work. The deterministic
`build_sparse_ilut`/`apply_sparse_ilut` path now supplies drop tolerance,
per-column fill selection, fixed-factor JVP/VJP, and explicit pivot status;
`build_sparse_ichol` supplies the matching SPD drop/fill path on the existing
sparse IC apply/JVP/VJP contract; `solve_sparse` now accepts the
`ichol_controlled` preconditioner with an exact full-fill one-step PCG oracle.
The public `build_sparse_ilut_row` constructor now performs no-pivot ILUT
against retained row factors with O(n + nnz) temporary work storage and an
independent no-fill diagonal oracle; it exports the same sparse CSC apply and
fixed-factor derivative contract as the dense reference builder. Controlled
row-oriented ICHOL and measured scaling remain. The
converged-state PCG JVP/VJP differentiates the exact solve independently of
the inactive preconditioner iteration path; factor rebuilds, breakdowns, and
graph changes are reported events rather than silently differentiated.

Nonlinear infrastructure includes Newton, damped Newton, Newton-Krylov,
pseudo-transient continuation, trust regions, deflation, parameter
continuation, bifurcation diagnostics, and explicit failure reasons.

Constraints include flux, helicity, mass, current, cross-helicity, angular
momentum, gauge, mode normalization, interface volume, and wall-current
conditions. Nullspace support covers constants, electromagnetic gauges, rigid
modes, harmonic representatives, and cohomology bases on multiply connected
domains.

### 9.2 FortSym and Enzyme pipeline

Every new formulation follows this sequence:

1. Define fields, spaces, geometry, and strong equations in FortSym.
2. Integrate by parts and expose volume, boundary, interface, and
   distributional terms.
3. Derive residual, Jacobian, JVP, VJP, and parameter or shape derivatives.
4. Simplify tensor symmetries and generate Fortran kernels.
5. Generate manufactured forcing, an independently simplified oracle, and
   equation documentation.
6. Check generated output byte for byte against the pinned source revision.
7. Use Enzyme as a second derivative route where compiler support is available.
8. Compare analytical derivatives with finite differences away from topology
   events, using step-size sweeps rather than one arbitrary step.

Generated kernels cover basis values, gradients, curls, divergences, Hessians,
Piola maps, traces, jumps, singular limits, cut-cell moments, Fourier triads,
special-function products, and geometry derivatives. Transient generated
files, plots, and compiler output are gitignored.

## 10. External application parity matrix

The following table defines the generic numerical foundation that an external
application could exercise. It does not authorize copying source code, adding
plasma-specific readers or closures, or shipping a replacement for the named
code. A row is complete when the generic contract and an independent oracle
work. Input conversion and application physics remain outside FortFEM.

| Reference target | Physics class | FortFEM ingredient parity target |
| --- | --- | --- |
| [CHEASE](https://crppwww.epfl.ch/~sauter/chease/), [paper](https://doi.org/10.1016/0010-4655(96)00046-X) | 2D fixed-boundary axisymmetric toroidal equilibrium | Generic axisymmetric elliptic/Fourier forms, spline/FEM geometry, axis and boundary trace contracts, and a common sampler. No COCOS or GEQDSK implementation |
| [FreeGS](https://freegs.readthedocs.io/en/stable/creating_equilibria.html) | 2D free-boundary axisymmetric equilibrium | Generic nonlinear residual, external boundary and coil-trace callbacks, X/O-point metadata fields, and manufactured profiles. Coil physics and GEQDSK conversion remain external |
| [VMEC/PARVMEC](https://github.com/ORNL-Fusion/PARVMEC), [VMEC++ numerics](https://arxiv.org/abs/2502.04374) | 3D nested-surface variational ideal equilibrium | Fourier angles, radial FE/IGA, generic energy and constraint blocks, shape JVP/VJP, and an external-data sampler |
| [GVEC](https://gvec.readthedocs.io/develop/index.html), [DESC](https://github.com/PlasmaControl/DESC) | Flexible 3D variational equilibrium and optimization | General coordinate maps, radial B-splines, Fourier modes, multiple interfaces, and exact residual derivatives. Input and profile models remain external |
| [SPEC](https://princetonuniversity.github.io/SPEC/) | Multi-region relaxed MHD with ideal interfaces | Region graph, independent fields, generic curl-eigenproblem and constraint blocks, total-pressure trace law, and interface shape derivatives. Beltrami and profile selection remain client code |
| [GPEC](https://princetonuniversity.github.io/GPEC/), [references](https://princetonuniversity.github.io/GPEC/references.html) | Linear ideal, kinetic, and resistive perturbed response | Fourier coupling, singular outer/inner layer contracts, vacuum and wall response, response matrices, normalization, and reciprocity. Equilibrium import remains external |
| [MARS-F response work](https://doi.org/10.1016/j.cpc.2006.09.003) | Linear toroidal ideal/resistive MHD and wall response | Linear block interfaces, complex frequency, generic plasma-vacuum-wall trace coupling, Fourier-FEM assembly, and resistive-layer matching. MARS physics remains external |
| [GLISS](https://github.com/itpplasma/GLISS) | Global linear ideal-MHD stability in 3D toroidal equilibria | Compatible radial spline FE, Fourier mode topology, generic eigenvalue and inertia contracts, and derivative boundaries. GVEC/VMEC input adapters remain external |
| [JOREK](https://www.jorek.eu/), [overview paper](https://arxiv.org/abs/2011.09120) | Nonlinear extended and resistive MHD | 2D FE plus toroidal Fourier blocks, coupled residuals, anisotropic transport, implicit structure-aware stepping, wall and free-boundary traces, operator-level parity tests |
| MEPHIT and STARWALL | Electromagnetic response and resistive-wall coupling | H(curl)/H(div) FEEC, surface traces, FEM/BEM/DtN wall blocks, retained complex factors, interface currents, reciprocity and energy tests |
| [BOUT++](https://bout-dev.readthedocs.io/en/stable/user_docs/introduction.html) | General 3D curvilinear fluid framework with model-specific clients | Equation-as-data residuals, curvilinear metric and boundary contracts, field-aligned operators, mixed conservative fluxes, and model-level JVP/VJP. Fluid models remain external |
| [GRILLIX](https://physik.uni-greifswald.de/ag-manz/forschung/codes/grillix/), [FCI paper](https://doi.org/10.1088/1361-6587/aaa373) | 3D edge and scrape-off-layer use of flux-coordinate-independent operators | FCI field-line tracing and interpolation, parallel/perpendicular operator split, immersed boundaries, anisotropic transport, generic material traces, and manufactured MMS fixtures. Braginskii and sheath laws remain external |
| [PARALLAX](https://gitlab.mpcdf.mpg.de/phoenix-public/parallax), [elliptic solver paper](https://arxiv.org/abs/2509.11831) | FCI mesh, magnetic-field handling, 2D elliptic solves, matrix-free 3D actions, and multigrid for GRILLIX and GENE-X | A Fortran-compatible geometry and operator adapter, matrix-free sparse products, anisotropy-aware multigrid contracts, and independent Poisson/Ampere timings. PARALLAX is LGPL-3.0, so no source is copied into FortFEM. |
| Linear elasticity and wave FEM literature | Mixed stress-displacement elasticity and symplectic mixed acoustics | Elasticity-complex spaces, tensor pressure/stress, first-order wave state, port-Hamiltonian pairing, symplectic time step, and cross-physics manufactured tests |

The local `/home/ert/code/GLISS` checkout was inspected for this roadmap. Its
README identifies version 0.0.2 as MIT licensed, with fixed-boundary FEEC
spectra, Mercier diagnostics, symmetric GVEC/VMEC input, radial B-splines,
Fourier angular modes, and Enzyme-generated derivative actions. Its
`PROVENANCE.md` records source-level traceability and a compatibility policy.
The planned FortFEM integration is therefore an oracle and interchange layer,
not a source import.

The name “Brandon Shennahan” did not resolve to a unique plasma-code project
in the local checkouts or the public search results. The relevant public
lineage is nevertheless clear enough to track: BOUT++ provides a general
curvilinear fluid framework, GRILLIX provides flux-coordinate-independent
edge/SOL turbulence, and PARALLAX provides the shared FCI geometry and
elliptic-solver layer. If a specific Shennahan project is intended, add its
repository or paper to this table before claiming parity.

## 11. Ordered example gallery

Examples are generated from one command, run in increasing complexity, and
publish plots plus machine-readable CSV or JSON. Images are build artifacts,
not committed source files. Every example includes at least a field plot, a
mesh or patch plot, an error or convergence plot, and a derivative or
conservation diagnostic when applicable.

### Foundations

1. **1D Poisson and Sturm--Liouville.** P1, high-order, exact polynomial
   reproduction, JVP/VJP, and a simple line plot.
2. **1D Legendre and spherical functions.** Ordinary, associated, and
   half-integer functions with recurrence and differential-equation oracles.
3. **2D simple Poisson.** First gallery entry, with mesh, solution, error,
   and timing plots.
4. **2D anisotropic diffusion.** Tensor coefficients, flux conservation, and
   manufactured solution.
5. **2D FEM/BEM Laplace.** Circle and polygon transmission with a larger-domain
   comparison.

### Interfaces, vectors, and open boundaries

6. **Curved Poisson interface.** Fitted CG, XFEM, XIGA, Nitsche, SIPG, and a
   regularized interface compared against a known jump solution.
7. **Surface delta source.** Explicit interface integral, duplicated traces,
   and narrow regularized source.
8. **H(curl) Ampere jump.** Analytical tangential magnetic jump and explicit
   surface current \(\mathbf K\), with Nitsche and fitted alternatives.
9. **H(div) magnetic flux.** Exact normal-flux continuity and divergence
   diagnostics on a fitted and cut interface.
10. **Slab Helmholtz.** FEM, BEM, DtN, and PML against an analytical plane wave.
11. **Cylinder and sphere Helmholtz.** Bessel/Hankel and Mie oracles, DtN, PML,
    larger-domain comparison, and 1D radial plus 2D/3D plots.
12. **Torus scalar exterior problem.** Exact curved toroidal surface, torus
    harmonics, FEM/BEM, DtN, PML, and far-field error.
13. **Curl-curl torus scattering.** Nédélec FEM, RWG/BC BEM traces, vector DtN,
    PML, reciprocity, and Ampere performance data.

### Application-oriented numerical ingredients

14. **Cylindrical axisymmetric elliptic fixture.** Manufactured coefficients,
    axis regularity, Fourier-FEM, and optional CHEASE/FreeGS comparison data
    supplied through an external adapter.
15. **Fourier-FEM slab and cylinder.** Mode diagonal operators, retained
    nonlinear triads, real/conjugate packing, and torus-harmonic diagnostics.
16. **Multi-region curl-eigenproblem fixture.** Independent regions, ideal
    interfaces, generic flux and helicity constraints, and a pressure-balance
    residual with manufactured coefficients.
17. **Linear 3D perturbation blocks.** Generic mode response with vacuum, wall,
    singular-layer, and response-matrix toy operators. GPEC or MARS data enter
    only through an external sampler.
18. **Energy-principle toy spectrum.** Radial spline FE, Fourier modes, inertia
    count, eigenpair derivatives, and external-data interchange tests.
19. **Resonant interface sheet.** Ideal current-sheet limit, finite resistive
    layer, XFEM enrichment, fitted sheet, and convergence to a manufactured
    singular solution.
20. **Skew bracket fixture.** Energy-skew nonlinear bracket, Fourier convolution,
    analytical JVP, and long-time structure-preserving integration for a
    caller-defined state.
21. **Resistive diffusion and tearing layer.** Explicit layer, adaptive layer,
    asymptotic enrichment, DG, and ideal-limit comparison.
22. **Generic coupled-field path.** Independently testable magnetic, scalar,
    vector, tensor, interface, and constraint residual blocks. A JOREK-style
    client can map its fields into this path, but FortFEM contains no JOREK
    state or closure implementation.

### Waves, elasticity, and anisotropy

23. **Mixed symplectic acoustics.** Pressure and particle velocity in a first-
    order compatible pair, with energy-preserving symplectic stepping and a
    comparison against the second-order pressure reduction.
24. **General wave family.** The same mixed residual for acoustics, Maxwell,
    elastodynamics, and a scalar wave, with common energy, dispersion, and
    boundary-port plots.
25. **Structure-preserving linear elasticity.** Displacement, velocity, and
    tensor stress with weak symmetry, traction interfaces, and a mixed
    Hellinger--Reissner oracle.
26. **Tensor-pressure wave.** A caller-supplied anisotropic tensor with
    parallel, perpendicular, and gyrotropic parts, including force balance and
    energy diagnostics. The tensor is a generic constitutive fixture.
The current `cgl_pressure_tensor` gallery example provides the first
manufactured constitutive/force-divergence profile and CSV/1D FortPlot
outputs; the coupled wave and higher-dimensional cases remain active.
The `field_aligned_flux` gallery example now provides a generated
parallel/perpendicular profile and (k_\parallel/k_\perp=100) flux plot; a
full assembled anisotropic diffusion gallery case remains active.
27. **Field-aligned diffusion.** A slab, cylinder, and torus with extreme
    \(k_\parallel/k_\perp\), comparing aligned coordinates, FCI maps, Fourier-
    FEM, and an isotropic control case.
28. **Field-aligned edge operator.** A generic anisotropic transport system
    with caller-supplied coefficients, material traces, and a reproducible FCI
    field-line map. This is an operator fixture, not a GRILLIX, BOUT++, or
    Braginskii implementation.

The gallery must show the same case in 1D, 2D, and 3D where a dimensional
reduction exists. Plot scripts use FortPlot and must include mesh edges,
element orientation, patch boundaries, and internal surfaces without dropping
parts of a mesh.

## 12. Verification and benchmark hierarchy

### 12.1 FortFEM-internal levels

1. **Algebraic:** partition of unity, polynomial reproduction, orientation,
   trace signs, exact sequences, rank, nullspaces, and surface measures.
2. **Analytical patch:** constants, affine fields, exact jumps, singular
   enrichments, special functions, torus harmonics, and exact delta integrals.
3. **Manufactured:** FortSym-generated source, boundary data, interface data,
   residual, Jacobian, JVP, VJP, shape derivative, and convergence rates.
4. **Cross-formulation:** fitted CG, XFEM/XIGA, DG/HDG, Nitsche, mortar,
   explicit surface current, regularized layer, FEM/BEM, DtN, and PML on the
   same physical case.
5. **Conservation and structure:** energy, helicity, cross-helicity, flux,
   divergence, charge, reciprocity, passivity, symplectic form, reversibility,
   and dissipative balance.
6. **Performance:** assembly, factorization, solve, matrix-free action, memory,
   iteration count, conditioning, and derivative cost as functions of order,
   mesh size, mode count, and interface complexity.

Report \(L^2\), H1, H(curl), H(div), trace, jump, flux, energy, conservation,
and derivative errors. Report error versus degrees of freedom and versus wall
time. Use a controlled thread count and record compiler, flags, CPU, BLAS,
FortNum/FortSparse revision, and random seed.

### 12.2 External oracle policy

The external comparison matrix covers FreeFEM, MFEM, FEniCSx, deal.II,
GeoPDEs, PetIGA, CHEASE, FreeGS, VMEC/PARVMEC, GVEC, DESC, SPEC, GPEC,
MARS-F, GLISS, STARWALL, and JOREK where a license-safe fixture is available.

- Keep only tiny analytical inputs, adapters, expected tolerances, and a
  provenance manifest in FortFEM.
- Run lightweight FreeFEM, MFEM, and FEniCSx cases in optional jobs. Use `uv`
  for Python environments where practical.
- Do not require heavy, proprietary, or cluster-only packages in the GitHub
  Pages job. Their pre-generated numerical data and performance summaries go
  to a separate `fortfem-benchmarks` data repository, referenced by commit,
  checksum, license, and executable version.
- Never copy external source code or proprietary test data into FortFEM.
- Compare fields and invariants on a common physical sampling set, not by
  assuming that two codes use the same basis or numbering.

## 13. Performance and reproducibility

The performance path is Fortran, FortNum, and FortSparse only. Profile with
the platform tools before changing an algorithm. Track:

- element and cut-cell kernel throughput;
- Fourier transform and mode-convolution cost;
- sparse assembly and matrix-free products;
- direct factorization reuse;
- Krylov iterations and preconditioner setup;
- BEM near-field, far-field, and quadrature cost;
- PML layer overhead;
- JVP/VJP cost relative to primal;
- memory peak and allocation counts.

Small deterministic benchmarks run in focused tests. Larger benchmarks run in
scheduled or manually dispatched Actions, publish JSON and SVG/PNG artifacts,
and update the Pages gallery without committing images. A benchmark result is
not considered comparable unless hardware, compiler, thread count, source
revision, and external data revision are recorded.

Performance gates are workload-specific. The benchmark manifest records the
operator, discretization, order, mesh or mode count, tolerance, thread count,
and reference implementation. It reports time to assemble, time to apply or
solve, peak memory, allocation count, iteration count, derivative overhead,
and error at a common physical sampling set. Reference runs may use MFEM,
Firedrake, NGSolve, deal.II, PETSc, Bempp-cl, or a small analytical kernel,
depending on the operator. The target is a measured improvement on the
declared workload, with accuracy and structure-preservation constraints held
fixed. The project does not use an unqualified speed claim as a substitute for
the benchmark record.

The execution plan is staged:

1. Keep all element and residual kernels local, pure where possible, and
   independent of MPI. Use batched quadrature, generated contractions,
   static condensation, matrix-free actions, and reusable sparse factors on a
   single node.
2. Add the serial ownership and ghost interfaces, deterministic reductions,
   partition-independent numbering, and a no-op communicator backend to the
   core data structures.
3. Add optional MPI halo exchange, distributed assembly, distributed vectors
   and matrices, checkpoint/restart, and solver adapters after the serial
   contracts have independent analytical and cross-code oracles.
4. Add distributed BEM acceleration, domain decomposition, GPU kernels, and
   large-scale strong and weak scaling only when a representative workload
   demonstrates the need.

The single-node path remains a complete supported path at every stage. The
MPI path must reproduce the serial residual, derivative, conservation, and
structure diagnostics within declared reduction tolerances.

## 14. Implementation phases

The phases are ordered by dependency. Each phase ends with a public API,
generated or analytical oracle, focused tests, documentation, and at least one
gallery example.

### Phase 0: Contracts and inventory: **active**

- Complete the primal/JVP/VJP inventory for FEM, BEM, DtN, PML, geometry,
  Fourier, special functions, sparse products, and all iterative solvers.
- GMRES and BiCGSTAB implicit state derivatives now pass finite-difference and
  adjoint-identity tests. Continue the inventory for remaining public
  operators and block solver compositions.
- Publish the complex-adjoint and shape-derivative conventions.
- Add seeded random generators and shrinking to the `fo` test path. The first
  generators cover oriented cells, positive quadrature weights, tensor
  coefficients, mode sets, and fixed-topology interface data.
- Define the partitionable data model before implementing a full MPI backend:
  global entity IDs, owner and ghost metadata, local-to-global maps,
  deterministic reductions, communicator-free local kernels, and serial
  no-op implementations of halo and reduction operations. Add MPI assembly,
  checkpointing, and distributed solver adapters after the serial residual,
  derivative, and invariant contracts are stable.
- Keep FortSym revision pins and generated-kernel checks green.

### Phase 1: Interface calculus: **active**

- The public orientation-safe scalar/vector trace contract now provides
  plus/minus averages and jumps, normal/tangential projections, and the
  rotated Ampere surface-current jump with an independent sign oracle.
- The explicit `assemble_surface_delta_load` primitive now assembles
  trace-basis transpose times positive surface quadrature/source weights,
  providing the fitted δ_\Gamma weak-load contract.
- `assemble_surface_vector_delta_load` adds the tangential trace/surface-
  current pairing needed for an explicit Ampère sheet. Scalar and vector
  delta-load actions now also expose analytic fixed-topology JVP/VJP products
  for basis, quadrature weights, and source/current values; this keeps explicit
  sheet residuals differentiable for fitted, cut, DG, and IGA callers.
- `assemble_interface_surface_current` now evaluates the oriented Ampere
  jump, its integrated current ledger, and fixed-topology JVP/VJP actions.
  Independent analytical, orientation-reversal, finite-difference, and
  real-adjoint tests cover the generic trace algebra; conservation at
  interface edges and material laws remain higher-level contracts.
- `assemble_surface_current_trace_residual` now pairs caller-owned,
  independently sized test and trial tangential trace bases with a target
  current. Its full product-rule JVP and real VJP cover basis, quadrature,
  coefficients, and target data, with a direct vector oracle. This is the
  neutral ownership block for fitted duplicated traces, cut/XFEM or XIGA,
  DG/HDG skeletons, and IGA patches; constitutive pressure laws and
  flux/helicity constraints remain client composition.
- `assemble_interface_jump_penalty` assembles the symmetric positive-
  semidefinite plus/minus jump block used by SIPG and Nitsche penalty terms.
- `assemble_symmetric_nitsche_interface` adds the average-flux consistency
  terms under the same orientation convention. Its value/JVP/VJP actions
  differentiate traces, fluxes, surface weights, and penalty parameters;
  scalar SIPG and non-symmetric consistency variants remain.
- `assemble_mortar_trace_coupling` supplies a weighted cross-mass block for
  independently discretized trace spaces.
  `assemble_scalar_sipg_interface` now composes independently sized test and
  trial trace/flux spaces with symmetric, incomplete, or nonsymmetric
  consistency (`theta=1,0,-1`) and value/JVP/VJP actions; vector FEEC
  consistency, HDG trace blocks, local static condensation, and signed-map
  sparse CSC accumulation are now public. The compatible flux elimination
  primitive now provides a differentiable local recovery map and condensed
  state block; global signed sparse flux ownership remains client-owned.
  `assemble_vector_jump_penalty` now supplies the tensor-valued counterpart
  for caller-owned tangential/normal projectors and anisotropic metrics, with
  value/JVP/VJP actions for vector traces, metrics, weights, and penalties.
  `assemble_vector_sipg_interface` now adds independent test/trial vector
  average-flux consistency for `theta=1,0,-1`, including metric-aware
  value/JVP/VJP actions; higher-order HDG trace composition and concrete
  enriched-space construction remain.
  `assemble_hdg_static_condensation` now exposes the local Schur complement
  and condensed load with implicit-solve value/JVP/VJP actions, so global
  skeleton assembly can be layered without differentiating pivot choices.
  `assemble_compatible_flux_elimination` now exposes the corresponding
  flux-specific recovery map, condensed state block, and full real JVP/VJP;
  `assemble_feec_commuting_projection` now audits all three projected
  differentials with complete value/JVP/VJP actions; actual higher-order
  enriched-space construction and client-owned projection maps remain.
  `assemble_scalar_numerical_flux` now provides conservative central, upwind,
  and Lax--Friedrichs choices with a quadratic-entropy diagnostic and complete
  fixed-topology value/JVP/VJP actions; `assemble_vector_numerical_flux` now
  lifts this to metric-weighted vector states, including strongly anisotropic
  constitutive tensors and complete value/JVP/VJP actions. Entropy-stable
  systems remain.
- The neutral `cell_complex_t` contract now stores oriented integer chain
  boundary maps, checks both boundary-of-boundary identities exactly, and
  reports Euler characteristic and compact Betti diagnostics. Independent
  interval, loop, sphere-CW, torus-CW, and malformed-orientation tests cover
  the primitive; quotient boundary maps, cycle/cocycle kernels, and metric
  harmonic one-forms are now public. `normalize_harmonic_one_forms` now
  supplies fixed-topology period/flux normalization with a central-difference
  JVP and real-adjoint oracle; cycle labels, physical units, gauge constraints,
  and cycle ledgers remain higher graph layers.
- The neutral `region_interface_graph_t` contract now adds oriented plus/minus
  region incidence, periodic self-identifications, and compact connectivity
  labels plus a spanning-forest cycle basis satisfying the exact integer
  incidence nullspace. Independent chain, disconnected, reversed, periodic,
  triangle-cycle, and malformed-endpoint tests cover the graph. Surface laws,
  sheet-current balances, and region physics remain application-owned.
- The neutral `cell_identification_t` contract now validates idempotent
  canonical representatives and signed orientations and reports compact
  quotient classes. `identify_boundary_matrix` pushes oriented incidence to
  quotient classes and rejects inconsistent identified columns. Independent
  identity, signed-periodic, interval-to-circle, signed-column,
  inconsistency-rejection, cycle-rejection, canonical-sign, zero-sign, and
  shape-mismatch tests cover the metadata; quotient geometry remains
  higher-level work.
- The neutral tree--cotree gauge contract now selects a spanning forest on an
  oriented graph, exposes cotree restriction/prolongation, extracts the
  fixed-gauge dense direct system, and composes with real/complex CSC direct
  solves through the existing constrained reduction. Independent triangle,
  disconnected-forest, reduction, sparse solve, fixed-map derivative, and
  malformed-incidence tests cover the selector. `build_tree_cotree_dof_map`
  now lifts the frozen control-edge selector to arbitrary high-order or IGA
  global DOF numbering, retaining extra edge/face/cell moments for the caller's
  sparse direct solver; duplicate and out-of-range maps have independent
  rejection tests. Period constraints remain a separate composition layer.
- The neutral `internal_manifold_graph_t` contract now records oriented
  plus/minus region sides, open or boundaryless manifold endpoints, periodic
  self-identifications, junction incidence, closed-manifold flags, and
  manifold connectivity. Independent slab, cylinder, sphere, torus, and
  malformed-endpoint fixtures cover the metadata; geometry, trace spaces,
  surface divergence, sheet-current balance, and pressure laws remain later
  composition layers.
- The integrated surface-current junction ledger now consumes the fixed
  manifold boundary incidence and distributes each manifold current to open
  junctions, with a zero global balance oracle for conservative columns and
  fixed-topology JVP/VJP actions. Geometric edge flux, constitutive current
  laws, and physical pressure/flux laws remain application composition.
- The fixed-topology loop-current constraint now applies an integer cycle basis
  to integrated manifold currents and subtracts caller-owned target values.
  Its residual, JVP, and VJP have an independent cycle oracle. Open-edge
  junction balance and closed-loop constraints are therefore separate linear
  ledgers; geometric surface divergence, constitutive current laws, and
  physical pressure/flux laws remain composition layers.
- The topology-only surface edge-flux ledger now applies an oriented
  vertex-edge boundary map to caller-integrated tangential edge fluxes. It
  exposes vertex divergence, a global conservation scalar, and fixed-topology
  JVP/VJP actions with open-chain and closed-cycle oracles. Conormal geometry,
  edge quadrature, and the physical tangential-current law remain composition
  layers. `assemble_surface_edge_flux` now supplies the missing neutral
  geometry-to-edge contraction from oriented conormal quadrature and surface
  current, with product-rule JVP/VJP actions and independent orientation,
  finite-difference, and adjoint tests. Surface-basis ownership and physical
  current laws remain composition layers; the neutral test/trial trace
  residual is now public for that ownership layer.
- The neutral normal-traction jump now projects caller-supplied plus/minus
  traction vectors onto a validated unit normal and subtracts a caller-owned
  target. Its product-rule JVP/VJP oracle composes generated CGL, elastic, or
  Maxwell-stress blocks without selecting a constitutive or plasma law.
- Oriented triangle surface measures (area plus unit normal) now have a public
  JVP/VJP API with shared-vertex accumulation and independent finite-difference
  and dot-product oracles. A linear 2D triangle level-set cut primitive now
  returns edge intersections, physical segment length, and gradient normal with
  an affine independent oracle. Exact positive/negative subcell areas and
  interface-length consistency are now available for the same linear cut. A
  degree-one clipped-polygon quadrature primitive now adds exact positive and
  negative centroids, oriented normal data, and an affine-integration oracle.
  The fixed-topology 2D level-set interface JVP now differentiates edge
  intersections, segment length, and the normalized physical normal with an
  independent central-difference and topology-event oracle. The matching
  fixed-topology cut-quadrature JVP propagates edge intersections through
  positive/negative areas and centroids and composes the length/normal
  derivative, with central-difference and area/first-moment conservation
  oracles. A 3D tetrahedral level-set interface now returns ordered triangular
  or quadrilateral cut polygons, area, and gradient-oriented normal with
  independent plane/intersection tests. Exact positive/negative tetra cut
  volumes and centroids now close the degree-one volume/first-moment contract
  with analytic and conservation oracles. Fixed-topology tetra cut JVPs now
  propagate moving vertices and level values through clipped-face moments and
  interface area/normal with central-difference and conservation oracles. The
  internal-manifold graph now provides the topology attachment for this
  surface-measure contract; geometry-to-graph construction remains client
  composition work. The existing vector current pairing consumes the
  surface-measure contract.
- Broken H1, H(curl), H(div), and L2 spaces plus skeleton spaces.
- Explicit delta-source and surface-current weak terms.
- Fitted duplicated spaces, Nitsche, mortar, multipliers, and block constraints.
- Build the region and cell-complex graph around the validated chain maps,
  adding periodic identifications, harmonic representatives, gauge constraints,
  and stable cycle IDs.
- Add closed-loop and open-edge sheet-current constraints, surface divergence,
  pressure balance, and current-ledger oracles on slab, cylinder, sphere, and
  torus fixtures.

### Phase 2: Cut geometry and XFEM/XIGA: **active**

- The shifted Heaviside partition-of-unity primitive is now public. It
  returns the sign-shifted enrichment, has fixed-sign zero JVP/VJP actions,
  and rejects a zero level value as a topology event. Independent sign and
  derivative oracles cover the piecewise-smooth contract. The matching
  scalar product composition `N_i*(H(phi)-H(phi_i))` now has value/JVP/VJP
  actions and a real adjoint oracle. The same activation is now public for
  componentwise vector bases, so H(curl)/H(div) values can be composed while
  retaining the fixed-topology JVP/VJP contract. The corrected-XFEM blending
  ramp `r=sum_i a_i*N_i` and `r*Psi` value/JVP/VJP composition are now public
  with a full-enrichment reproduction oracle. The same ramp and reverse
  contraction are public for vector enrichment arrays, ready for covariant or
  contravariant Piola values and IGA coefficient blocks. Cut-cell geometry,
  support activation, rank/conditioning diagnostics, and Piola-aware vector
  differential operators remain.
- `evaluate_vector_enrichment_differential_3d` now exposes the explicit
  product-rule terms `curl(psi*b)=psi*curl(b)+grad(psi) x b` and
  `div(psi*b)=psi*div(b)+grad(psi) . b`, with full input JVP/VJP actions and
  independent curl/divergence and adjoint oracles. It is a neutral diagnostic
  for H(curl), H(div), Piola, and IGA composition; exact-sequence preservation
  or intentional jump reporting remains a higher-level space contract.
- Cut-cell classification and high-order quadrature. Exact degree-one triangle
  and tetrahedron rules plus exact degree-two raw-moment tensors with
  fixed-topology JVPs are now public. The 2-D linear level-set path now also
  exposes exact degree-three raw-moment tensors and fixed-topology JVPs for
  positive and negative clipped polygons, with an independent simplex-moment
  and central-difference oracle. The matching 3-D tetrahedral path now exposes
  exact degree-three raw-moment tensors and fixed-topology JVPs through its
  oriented tetrahedral fan, with an independent simplex-moment and
  central-difference oracle; higher-order cut rules remain. A weighted
  enrichment-support Gram contract now exposes value/JVP/VJP actions with a
  fixed activation mask, and a symmetric-Jacobi rank/condition diagnostic
  reports the active enrichment rank and singularity with an independent
  weighted, finite-difference, and real-adjoint oracle; support activation maps
  now expose CSR sign activation, unique extrema, signed margins, and
  fixed-topology JVP/VJP actions with explicit owner/tie/topology rejection;
  the physical vector/metric support Gram contraction now composes with
  covariant or contravariant Piola and IGA values and exposes value/JVP/VJP
  actions with SPD-metric rejection. The batched 2D/3D covariant and
  contravariant Piola-enrichment value/JVP/VJP composition now uses FortNum’s
  FortSym-derived inverse/determinant products and rejects singular or
  orientation-reversing maps; exact-sequence cut spaces, higher-order curved
  maps, and commuting conditioning remain.
- Heaviside, kink, singular, helical, and resonant enrichments.
- Shifted bases, corrected blending elements, pruning, conditioning, and
  connected-component activation.
- Trimmed B-spline stabilization and cut H(curl)/H(div) extensions.
- Verify the commuting diagram for every enriched vector space and document
  the exact sequence identity that a physical jump intentionally changes.

### Phase 3: DG and HDG: **active**

- The public mortar trace cross-mass block now has value/JVP/VJP actions for
  independently discretized skeleton traces. The scalar SIPG interface block
  now has symmetric, incomplete, and nonsymmetric consistency variants with
  value/JVP/VJP actions; vector FEEC consistency, HDG traces, local
  hybridization, signed-map sparse CSC accumulation, and compatible flux
  elimination are now public. The latter returns the local compatible flux
  recovery map plus condensed state system with full real JVP/VJP actions;
  global sparse flux ownership remains client-owned.
  A tensor-weighted vector jump penalty now covers component-valued traces and
  anisotropic metric directions; vector consistency fluxes and FEEC commuting
  projections are now public with value/JVP/VJP actions; hybridization and
  static condensation now have a local differentiable Schur primitive;
  `assemble_hdg_global_skeleton` now supplies a signed-map dense reference
  assembler with value/JVP/VJP actions; sparse accumulation remains.
  `assemble_feec_exact_sequence` now exposes the
  metric-independent `curl(grad)` and `div(curl)` defects with independent
  value/JVP/VJP actions for simplicial, IGA, multipatch, and periodic maps.
  `assemble_feec_commuting_projection` now checks projected discrete and
  continuous differentials, including all projection directions in its
  JVP/VJP; concrete generated enriched-space constructors remain a later
  layer.
  The symmetric jump-penalty block also has value/JVP/VJP actions, including
  penalty and surface-weight directions.
- Conservative/upwind vector flux interfaces now include scalar and
  metric-weighted vector value/JVP/VJP paths; entropy-stable system fluxes
  now have an explicit SPD-metric wrapper with fixed-topology value/JVP/VJP
  actions; nonlinear system entropy variables and characteristic HLL/HLLC
  laws remain application-owned.
- Broken vector FEEC, hybridization, static condensation, and mixed CG-DG.

### Phase 4: Open boundaries and vector operators: **active**

- Complete scalar and Maxwell FEM/BEM, DtN, and PML parity on slab, circle,
  sphere, cylinder, and torus.
- Add exact-curved-torus external surfaces and vector DtN.
- Add larger-domain, BEM, DtN, and PML comparisons for Poisson, Ampere, and
  Helmholtz.
- Verify FortPlot mesh and surface rendering for every mesh-bearing plot.

### Phase 5: Fourier-FEM and torus harmonics: **active**

- The public `fourier_mode_registry_t` records field-period phase,
  normalization, real/conjugate packing, retained triads, and caller-selected
  radial regularity. Its analytical mode value, fixed-topology JVP, and
  complex real-part VJP are independently tested. This metadata layer is
  neutral and does not import equilibrium profiles, COCOS/GEQDSK conventions,
  or plasma closures; mode-coupled operators and torus-harmonic radial
  branches remain composition work.
- `assemble_fourier_vector_product` now contracts a caller-owned complex
  coupling tensor over every retained triad for scalar, vector, tensor,
  H(curl), or H(div) coefficient arrays. Its JVP differentiates coupling and
  both factors, while its complex real-part VJP has an independent dot-product
  oracle; truncated triads and de-aliasing remain explicit caller policy.
- FortNum `legendre_q` and `legendre_q_derivative` are now public on `main`,
  with closed-form Q0--Q3, centered-difference derivative, and invalid-domain
  tests. Its toroidal P/Q branch is likewise covered by independent value and
  ODE checks. Standard orthonormal complex spherical harmonics and angular
  derivatives are now public with closed-form degree-one, conjugacy, pole,
  and invalid-angle tests. Complete the remaining stable high-order and
  near-cut continuation before using these functions as production
  DtN/torus-harmonic building blocks.
- The FortFEM public API now re-exports the pinned toroidal P/Q values and
  derivatives through `fortfem_toroidal_harmonics`; an independent adapter
  oracle covers both Hobson branches, derivative values, and invalid-domain
  NaN behavior. Stable continuation and mode-coupled toroidal operators remain.
- Define mode normalization, phase, field-period, and real packing.
- Add mode-coupled scalar, H(curl), H(div), and caller-defined nonlinear
  operators.

### Phase 6: Structure-preserving time evolution: **active**

- The public `advance_mixed_wave_midpoint` step now provides the common
  first-order pressure/velocity, displacement/momentum, and port-Hamiltonian
  Cayley contract. Its independent test checks the oscillator map, energy, and
  signed-step reversibility.
- The public `advance_mixed_wave_symplectic_euler` step now provides a
  partitioned first-order symplectic update for the same mixed state. Its
  independent test checks the two-stage mass-solve oracle and the canonical
  two-state symplectic-form identity; dissipative terms remain separate.
- The public `advance_dissipative_cayley` step now provides the separate
  `(M+hD/2)^(-1)(M-hD/2)` update for positive-time dissipative blocks. Its
  JVP/VJP differentiate mass, damping, step size, and state, while an
  independent two-by-two oracle checks non-increasing SPD-mass energy. It is
  intentionally not labeled symplectic; splitting and absorbing clients own
  the composition.
- Variational/symplectic and Poisson building blocks for ideal terms.
- Energy-dissipative integrators for resistive and viscous terms.
- Symmetric splitting, implicit midpoint/Cayley, discrete-gradient or
  average-vector-field options, and long-time invariant tests.
- Mixed first-order pressure-velocity, displacement-velocity-stress, and
  electromagnetic wave states with a common port-Hamiltonian interface.
- Tensor-valued coefficients and anisotropic constitutive blocks with exact power,
  momentum, and stress-work diagnostics.
- The generated CGL-shaped tensor and product-rule force-divergence blocks are
  public and tested as generic constitutive fixtures. The tensor now also
  feeds compositional traction and `P:grad(v)` stress-work value/JVP/VJP
  diagnostics. `assemble_tensor_volume_work` now provides the neutral full
  volume `stress:grad(test)` residual with tensor, geometry, and quadrature
  JVP/VJP actions and an independent contraction oracle; constitutive laws and
  caller-supplied correction tensors remain external. No plasma closure is
  implemented.

### Phase 7: Equilibrium and linear-response ingredients: **active**

- The neutral `equilibrium_interchange_t` record now defines the external
  adapter boundary for mapped/physical samples, named coefficient and profile
  arrays, explicit scalar/vector/tensor component ranks and offsets, segmented
  boundaries, provenance, units, and normalization. It is validated with
  analytical toroidal and 2D manufactured fixtures plus deep-copy and
  rejection checks. FortFEM still does not implement GEQDSK, COCOS, or any
  plasma-specific reader, profile law, coil model, or equilibrium solver.
- Generic axisymmetric elliptic and Fourier variational fixtures that can be
  sampled by external CHEASE or FreeGS adapters. FortFEM will not implement
  GEQDSK or COCOS readers, profile laws, coil models, or equilibrium solvers.
- Generic multi-region curl-eigenproblem, interface, and constraint blocks that
  external VMEC, GVEC, DESC, or SPEC clients can exercise.
- Generic linear ideal and resistive response blocks, singular layers, vacuum,
  conducting-wall traces, and response matrices for external GPEC, MARS, and
  GLISS oracle data.
- Compose a small manufactured multi-field state from independently verified
  field, tensor, surface-current, vacuum, wall, and constraint blocks. Plasma
  state selection and closure remain outside FortFEM.

### Phase 7a: Field-aligned edge and SOL ingredients: **active**

- FCI RK4 field-line tracing, support-operator gradient/divergence algebra, and
  matrix-free (P K_\parallel Q) diffusion (including FortSym-generated
  gradient and full diffusion JVP/VJP products) are public and independently
  tested. The
  `fci_parallel_diffusion` gallery fixture provides a reproducible open-line
  mass-conservation and negative-energy profile. A 1D linear map builder and
  its fixed-topology JVP/VJP provide independent partition/affine and
  dot-product oracles, and a 2D Cartesian bilinear builder now covers a
  genuine poloidal slice. A generated quadratic Lagrange map covers explicit
  three-node nonuniform stencils on 1D slices, including fixed-stencil JVP/VJP
  products. The RK4 tracer also has a tangent callback path with an exponential
  endpoint oracle. The positive staggered
  flux-box volume constructor
  now covers the traced expansion/area/(B_\varphi) product with pinned
  FortSym-generated value/JVP/VJP kernels. Fixed-cell 2D map JVP/VJP products
  now cover source and target coordinate motion; generated quadratic and
  cubic and quartic fixed-stencil maps now expose value/JVP/VJP actions with
  independent polynomial, finite-difference, and real-adjoint oracles. A
  generated quintic fixed-stencil map now exposes the same value/JVP/VJP
  actions with independent quintic, finite-difference, and real-adjoint
  oracles. A generated fixed-topology quadrilateral area map now supplies
  positive unstructured plane-cell measures with independent shoelace,
  finite-difference, and real-adjoint oracles plus a gallery fixture. Higher
  interpolation derivatives beyond quintic and curved support-volume measures
  beyond quadratic remain separate planned components. The generic polygon map
  covers fixed-topology cells with more than four vertices, and its generated
  quadratic Bezier-edge extension covers arbitrary curved polygon boundaries. A
  generated quadratic Bezier-edge area map now
  supplies a fixed-topology curved-cell value/JVP/VJP contract with an
  independent Gauss--Green oracle and sampled-boundary gallery output.
- The batched 2D bilinear endpoint-to-map adapter now connects traced
  forward/backward endpoints to the support-operator tensor contract and
  carries fixed-topology source-grid and endpoint JVP/VJP actions. Moving
  connectivity, higher-order curved measures beyond quadratic, and stencil rebuilds at
  topology events remain planned; the generic straight polygon map and
  quadratic Bezier-edge quadrilateral map are now the unstructured-cell
  baselines.
- The fixed-stencil quadratic 1D map now has a batched segment adapter with
  generated-kernel-backed JVP/VJP accumulation and independent polynomial,
  finite-difference, and adjoint tests. Degenerate local source coordinates
  are rejected before evaluating the Lagrange products. Curved or moving-cell
  connectivity remains a topology-rebuild concern.
- A positive FCI diffusion diagonal is public as the first anisotropy-aware
  preconditioner contract, with a matching Jacobi apply and positivity test;
  PARALLAX-compatible plane multigrid and additive field splitting are now
  public baselines, while coupled stronger blocks remain planned.
- A PARALLAX-style anisotropic split action now combines independent per-plane
  CSC elliptic blocks with the conservative FCI parallel operator; tensor
  coefficient assembly and coupled stronger field splitting remain planned.
  `compute_fci_anisotropic_diffusion_diagonal` now combines the
  positive support diagonal with each plane block's oriented diagonal and
  rejects non-positive results, providing an independently tested scalar
  preconditioning oracle. `apply_fci_anisotropic_jacobi_preconditioner` applies
  that diagonal directly for small matrix-free solves; cached diagonal use and
  `apply_fci_plane_two_level_vcycle` now provides the next plane-solver layer;
  a public recursive multilevel V(1,1) path now covers arbitrary level sizes;
  the retained-factor path and additive field-split composition are tested
  separately. A public recursive W(1,1) variant now repeats coarse visits;
  coupled blocks remain active.
- The split action now also has a field-only VJP that composes the conservative
  FCI transpose with an explicit transpose of every plane CSC block. An
  independent nonsymmetric-plane oracle and real dot-product test guard this
  adjoint contract; coefficient and geometry derivatives remain separate.
- A batched `apply_fci_plane_two_level_vcycles` adapter now applies the
  independently tested two-level plane cycle to a homogeneous plane stack,
  preserving the PARALLAX field-split boundary between per-plane elliptic
  solves and the global FCI line action. A ragged-offset companion now covers
  variable plane sizes without padding. The public additive field-split
  preconditioner combines cached parallel Jacobi and ragged plane cycles with
  explicit nonnegative weights. A public retained coarse-factor path now
  reuses FortSparse factorizations across right-hand sides, and the recursive
  multilevel V-cycle accepts flat level offsets for nonuniform hierarchies, and
  the W-cycle variant repeats coarse corrections; coupled blocks remain active
  work.
- A fixed-cell barycentric triangle interpolation path now covers logically
  unstructured poloidal targets, including geometry and target JVP/VJP actions;
  its batched endpoint-to-map path now feeds the support-operator tensor
  contract; moving-cell connectivity and higher-order stencils remain planned.
- FCI field-line maps, higher-dimensional interpolation Jacobians, and parallel
  derivative JVPs.
- A geometry-only 2D FCI first-hit search now returns the nearest transverse
  oriented wall/target segment, exact hit parameter/point, and facet normal;
  valid no-hit traces and malformed facets have explicit status contracts. Its
  fixed-topology JVP differentiates hit point, parameter, and facet normal
  against central differences. `assemble_fci_terminal_boundary_flux` now
  provides the volume-weighted conservative contribution with fixed-owner JVP
  and VJP dot-product oracles; owner remaps and physical material laws remain
  separate contracts. The matching 3D triangle search and fixed-topology JVP
  now cover triangulated toroidal vessel/divertor surfaces without copying a
  PARALLAX implementation. A traced-polyline wrapper now scans those facets in
  field-line order and returns the first segment/triangle, endpoint, oriented
  normal, connection length, and fixed-event JVP; no-hit paths return their
  final endpoint and total length.
- The generic nonlinear material-surface flux contract is now public: an
  application callback supplies the local tagged law while FortFEM assembles
  the oriented trace residual, integrated per-tag ledger, and fixed-tag
  value/JVP/VJP products with independent central-difference and dot-product
  oracles. Sheath, Bohm, recycling, radiation, and material databases remain
  application-layer clients.
- Strongly anisotropic diffusion, conduction, resistivity, and viscosity.
- Treat material boundaries as explicit geometry rather than only diffuse
  penalization. [#58](https://github.com/lazy-fortran/fortfem/issues/58)
  owns exact first-hit FCI terminal segments on oriented wall/target facets;
  [#59](https://github.com/lazy-fortran/fortfem/issues/59) owns the matching
  conservative terminal boundary-flux contribution and fixed-topology
  derivatives.
- Keep physical boundary laws outside the numerical library while making them
  first-class residual terms. [#60](https://github.com/lazy-fortran/fortfem/issues/60)
  owns the generic nonlinear material-surface flux value/JVP/VJP contract;
  Bohm, sheath, recycling, sputtering, radiation, and heat-transmission models
  are application-layer clients.
- Support an optional fitted boundary-layer patch without forcing the bulk FCI
  mesh to conform to the vessel. [#61](https://github.com/lazy-fortran/fortfem/issues/61)
  now has the neutral `assemble_fci_boundary_patch_mortar` cross-mass,
  constant-preserving transfer, weighted-adjoint, overlap-measure, and
  ownership-multiplicity contract with fixed-topology JVP/VJP products and
  independent matching, reversed, duplicate, zero-measure, rank-deficiency,
  finite-difference, and dot-product tests. Geometry construction, general
  moving Chimera meshes, and application boundary laws remain outside FortFEM.
- GORILLA owns reusable characteristic stepping and material events, not
  FortFEM: [GORILLA #80](https://github.com/itpplasma/GORILLA/issues/80)
  introduces a common characteristic event contract and
  [GORILLA #81](https://github.com/itpplasma/GORILLA/issues/81) adds neutral
  free flight without misusing guiding-centre equations, while
  [GORILLA #82](https://github.com/itpplasma/GORILLA/issues/82) covers only
  those fluid-advection characteristics that genuinely admit marker tracing.
  Marker collisions, weights, and conservative cell/surface tallies belong
  above the pusher in
  [GORILLA_APPLETS #47](https://github.com/itpplasma/GORILLA_APPLETS/issues/47).
- Small MMS and performance cases aligned with PARALLAX, GRILLIX, and BOUT++
  concepts, without copying their implementations. The FCI gallery fixture
  now records the measured matrix-free action cost alongside its conservation
  and dissipation diagnostics; larger problem-size scaling remains active.
- Document the generic equation-as-data callback ABI for caller-owned fields,
  coefficients, sources, boundary laws, and target ledgers. Keep species,
  closures, sheath, and material physics in client code while testing units,
  signs, residuals, and balances through the public callback contract.

### Phase 8: Cross-code oracles and gallery: **active**

- Ordered examples with FortPlot plots, numerical data, convergence, and
  performance.
- Every example carries the smallest applicable complete contract: a
  manufactured or analytical solution, an independent oracle, primal and
  derivative checks, a conservation or structure diagnostic, a mesh or
  geometry view, and a timing and memory record. Examples progress from a
  one-dimensional Poisson patch test through two-dimensional interfaces and
  waves to three-dimensional torus, sphere, FEM/BEM/DtN, PML, anisotropic,
  and mixed structure-preserving cases.
- Gallery examples use the same generated kernels and public APIs as library
  clients. The gallery does not contain a faster special implementation that
  bypasses the tested residual or derivative path.
- Optional lightweight FEniCSx, FreeFEM, and MFEM runners.
- Sister-repository data for heavy or licensed references.
- GitHub Pages generation, link checks, and periodic deployment monitoring.
- The test workflow now continues through example generation after a failed
  verification stage, uploads its plot artifact, and reports the failure only
  after publication. The Pages workflow consumes the triggering test-run
  artifact even when the test conclusion is red, so known oracle failures do
  not hide the ordered gallery or its plots.
- Seeded property tests and a common sampler compare independent codes without
  assuming matching bases, numbering, or mesh topology.
- The gallery now requires every generated example page to embed a real PNG
  plot; the optional interoperability page retains an explicit SVG provenance
  cover until its external FEniCSx/FreeFEM/MFEM artifact is available. The
  formerly empty circular DtN, circle BEM/CFIE, transmission, mesh, solver,
  and tetrahedral Nédélec pages now emit numerical convergence, response,
  mesh, or performance plots. The test workflow uses the valid
  `fpm run --example --list` compile gate rather than the unsupported
  `fpm build --example` spelling.
- FortPlot mesh-bearing examples have a regression fixture that checks element
  count, boundary edges, patch boundaries, and internal surfaces in the
  rendered output before Pages deployment.
- `doc/examples/primary_plots.txt` is the gallery's explicit visual contract:
  each example names its physical solution, field, geometry, or mesh artifact.
  The documentation generator copies that artifact to `primary.png` and emits
  it before convergence, timing, conditioning, or other diagnostic plots. The
  generator now fails when a mapped artifact is absent instead of silently
  substituting a generated cover, and the docs test rejects diagnostic names
  such as convergence, error, timing, residual, accuracy, or DOF plots as the
  first view. FortPlot's 3-D scatter path must preserve per-sample colour maps
  so solution samples remain physically interpretable. CI collects nested
  per-example output with an independent media validator and keeps the short
  ten-second limit for focused tests while allowing the more expensive
  physical 3-D gallery products a bounded execution window.

### Phase 9: Future application layer: **reference only**

Use the stabilized ingredients to support future research applications that
may resemble GVEC, VMEC, SPEC, GPEC, MARS, GLISS, or JOREK. Do not start a
production code in this repository until the corresponding ingredient,
analytical example, and external oracle have passed the previous phases.
Steady edge/SOL/divertor transport, Braginskii closures, species and reaction
networks, neutral backends, impurity charge states, sheath and material laws,
GORILLA/EIRENE coupling, and target heat-load accounting belong in a separate
application repository. FortFEM owns only their reusable geometry, trace,
operator, residual/JVP/VJP, field-solver, and preconditioner foundations.

## 15. Definition of done

An ingredient is complete only when all applicable checks below pass:

- public API and units are documented;
- the residual and active constraints are explicit;
- FortSym-generated or independently derived kernel is provenance mapped;
- primal, JVP, VJP, and implicit solve derivatives have dot-product or
  finite-difference checks;
- geometry, trace orientation, and topology events are tested;
- an independent analytical or external oracle agrees within a declared
  tolerance;
- convergence, conditioning, conservation, and performance are reported;
- performance is compared with a named reference on a fixed workload, with
  hardware, compiler, thread count, memory, accuracy, and derivative cost
  recorded;
- structure-preserving properties are tested rather than asserted;
- seeded random properties pass with their seed, generated case, and shrink
  record retained in the test log;
- external-code comparisons carry a license, executable version, source
  revision, sampler, tolerance, and data checksum;
- focused `fo` tests meet the short feedback target;
- full CI passes on supported compilers;
- the example and documentation are generated reproducibly;
- transient output is ignored and no unrelated downstream fix is duplicated.

## 16. Non-goals and explicit boundaries

- No full VMEC, GVEC, SPEC, GPEC, MARS, GLISS, JOREK, CHEASE, or FreeGS
  replacement is planned in FortFEM.
- No GEQDSK or COCOS parser, equilibrium profile library, coil or vessel
  physics, plasma closure, species or reaction model, sheath, neutral,
  impurity, or divertor application is planned in FortFEM. These are external
  clients of the generic contracts.
- No C rewrite of numerical kernels is planned.
- No proprietary source, binary, or license-restricted benchmark data is
  checked into this repository.
- No claim of differentiability is made across a changing interface topology,
  adaptive mesh decision, enrichment activation, or eigenvalue crossing
  without a dedicated event or subgradient contract.
- No generic time integrator is accepted as structure preserving without an
  invariant or geometric test.
- No spacetime FEM implementation is in the current scope.

## 17. Literature and provenance

The following references motivate the architecture. Links are kept to papers,
official documentation, or official repositories where possible.

### Interfaces, enrichment, and IGA

- [Arnold, Falk, and Winther, finite element exterior calculus](https://arxiv.org/abs/0906.4325)
- [Fries, corrected XFEM blending elements](https://doi.org/10.1002/nme.2259)
- [Stable isogeometric analysis of trimmed geometries](https://www.sciencedirect.com/science/article/abs/pii/S0045782516308222)
- [XIGA for multi-material problems](https://doi.org/10.1007/s00466-022-02200-y)
- [Independent-field isogeometric boundary elements](https://arxiv.org/abs/1406.0306)
- [DefElement reference](https://defelement.org/)
- [MFEM features and finite-element reference implementation](https://mfem.org/features/)

### Special functions and toroidal harmonics

- [DLMF 14.3: definitions and hypergeometric forms](https://dlmf.nist.gov/14.3)
- [DLMF 14.10: Legendre recurrences and derivatives](https://dlmf.nist.gov/14.10)
- [DLMF 14.19: toroidal (half-integer) specialization](https://dlmf.nist.gov/14.19)
- [FortNum special-function API at the pinned revision](https://github.com/lazy-fortran/fortnum/blob/b5afac7/docs/api.md)

### MHD, equilibria, and linear response

- [CHEASE official page](https://crppwww.epfl.ch/~sauter/chease/)
- [CHEASE toroidal-equilibrium paper](https://doi.org/10.1016/0010-4655(96)00046-X)
- [FreeGS documentation](https://freegs.readthedocs.io/en/stable/creating_equilibria.html)
- [PARVMEC official repository](https://github.com/ORNL-Fusion/PARVMEC)
- [GVEC documentation](https://gvec.readthedocs.io/develop/index.html)
- [DESC official repository](https://github.com/PlasmaControl/DESC)
- [SPEC official documentation](https://princetonuniversity.github.io/SPEC/)
- [GPEC documentation and examples](https://princetonuniversity.github.io/GPEC/)
- [MARS-F response-model paper](https://doi.org/10.1016/j.cpc.2006.09.003)
- [GLISS repository](https://github.com/itpplasma/GLISS)
- [JOREK overview paper](https://arxiv.org/abs/2011.09120)
- [BOUT++ documentation](https://bout-dev.readthedocs.io/en/stable/user_docs/introduction.html)
- [GRILLIX official project page](https://physik.uni-greifswald.de/ag-manz/forschung/codes/grillix/)
- [GRILLIX FCI formulation](https://doi.org/10.1088/1361-6587/aaa373)
- [PARALLAX official repository](https://gitlab.mpcdf.mpg.de/phoenix-public/parallax)
- [PARALLAX/PAccX elliptic solver](https://arxiv.org/abs/2509.11831)
- [Mixed elasticity with weak stress symmetry](https://arxiv.org/abs/math/0701506)
- [JKU Linz TDNNS elasticity formulation](https://www.numa.uni-linz.ac.at/Talks/abstract/154/)
- [Anisotropic mixed finite elements for elasticity](https://doi.org/10.1002/nme.3319)
- [Symplectic mixed finite elements for acoustics](https://doi.org/10.1007/s00211-014-0667-4)
- [Port-Hamiltonian mixed finite elements for linear thermoelasticity](https://doi.org/10.1080/01495739.2021.1917322)
- [Structure-preserving linear elasticity example](https://publiweb.femto-st.fr/tntnet/entries/21775/documents/author/data)
- [HKT resonant current-sheet study](https://collaborate.princeton.edu/en/publications/numerical-study-of-%CE%B4-function-current-sheets-arising-from-resonan/)
- [Matched-asymptotic resistive-layer formulation](https://collaborate.princeton.edu/en/publications/computation-of-resistive-instabilities-by-matched-asymptotic-expa/)
- [STARWALL and linear MHD stability](https://arxiv.org/abs/1508.04911)

### Structure-preserving discretization

- [Splitting-based structure-preserving MHD discretization](https://doi.org/10.5802/smai-jcm.34)
- [Structure-preserving and helicity-conserving FEM for MHD](https://doi.org/10.1016/j.jcp.2021.110894)
- [Structure-preserving Hall-MHD finite elements](https://arxiv.org/abs/2202.11586)
- [Variational integrators in plasma physics](https://arxiv.org/abs/1307.5665)
- [Discrete variational integration for ideal MHD](https://www.osti.gov/servlets/purl/1179782/)
- [Structure-preserving transport-stabilized compatible FEM for MHD](https://doi.org/10.1016/j.jcp.2024.112737)

### FortFEM provenance

- [FortFEM differentiation design](doc/design/differentiation.md)
- [FortFEM FEEC, MEPHIT, and open-boundary audit](doc/roadmap_mephit_feec_bem.md)
- [FortFEM interoperability benchmark protocol](doc/interoperability_benchmarks.md)
- [FortFEM generated-kernel checker](tools/codegen/check_generated.sh)

When a future implementation uses a formula, convention, or compatibility
fixture from another code, add an equation-level entry to the relevant
provenance file and record the external revision. A literature citation alone
does not document a software dependency or justify source reuse.
