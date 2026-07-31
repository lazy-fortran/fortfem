# FortFEM roadmap

Status: living architecture and verification plan, 2026-07-31

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
  Run full `fo` suites and expensive examples asynchronously in CI or a
  controlled background job. Do not wait for a complete gallery build before
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

## 2. Current baseline

The following capabilities are already on FortFEM `main` or in the verified
documentation baseline. The list is intentionally conservative.

| Area | Current state | Next gate |
| --- | --- | --- |
| Scalar FEM | P1/P2/Q1 and arbitrary-order triangular scalar paths, Poisson and diffusion forms, boundary conditions, plotting | General field and coefficient callbacks in the symbolic form compiler |
| FEEC | Oriented triangular and tetrahedral H1, H(curl), H(div), and DG families, Piola maps, commuting tests, sparse assembly, mixed RT-DG Poisson | General multi-field block composition and arbitrary multipatch assembly |
| IGA | Nonuniform B-splines, rational maps, two- and three-dimensional de Rham incidence complexes, cylindrical and toroidal Fourier blocks, initial JOREK magnetic-flux residual/JVP | General patch graphs, trimming, enrichment, and the remaining coupled JOREK variables |
| Special functions | FortNum quadrature, Legendre and spherical functions, Bessel/Hankel paths, toroidal analytical utilities in active development | A documented, independently tested torus-harmonic API and stable half-integer continuation |
| Sparse algebra | FortSparse CSC assembly, retained factors, real and complex solves, sparse products, and CG, PCG, GMRES, and BiCGSTAB converged-state derivative contracts | Preconditioners with measured scaling, flexible Krylov products, and block solver derivatives |
| Open boundaries | Planar, circular, and spherical scalar Helmholtz DtN paths, scalar BEM, Maxwell trace and PML components | General curved Maxwell DtN, toroidal exterior maps, and robust FEM/BEM/DtN comparison fixtures |
| PML | Scalar and curl-curl Cartesian complex-stretching tensors with slab, triangular, and tetrahedral examples | Automated curved-object layers, reflection/error metrics, and derivative coverage for all geometry parameters |
| Differentiation | Analytical FortSym paths, selected Enzyme checks, sparse matrix products, converged CG/PCG/GMRES/BiCGSTAB solves, toroidal coordinate and DtN products | Complete operator inventory, JVP/VJP parity for all public operators, and shape derivatives |
| Examples | Generated documentation pages for Poisson, Maxwell, Helmholtz, BEM, IGA, torus, PML, and solver examples | Ordered gallery beginning with simple Poisson and adding 1D, 2D, and 3D plasma-oriented toy models |

The aggregate full test suite can expose resource-sensitive failures when many
independent numerical targets run in parallel. Focused repetitions are the
first diagnostic and remain the merge gate for a local increment. CI remains
the final cross-platform gate.

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
- tensor pressure in anisotropic plasma closures and reduced-MHD models;
- symmetric and nonsymmetric stress with explicit angular-momentum balance;
- traction and normal-flux traces on fitted, cut, and BEM interfaces;
- constitutive tensors with spatial, parameter, and magnetic-field dependence.

For a magnetized fluid, the CGL-style pressure tensor is a first test case,

\[
\mathbf P=p_\perp I+(p_\parallel-p_\perp)\,\mathbf b\mathbf b^T,
\qquad \mathbf b=\mathbf B/|\mathbf B|,
\]

with optional gyrotropic and Braginskii corrections. The parallel and
perpendicular pressure laws, their force divergence, and their work pairing
are separate residual terms.

The first constitutive slice is now public on `main`: FortSym generates the
six independent symmetric CGL components and their JVP/VJP, while the
`fortfem_cgl_pressure_tensor` wrapper validates the unit magnetic direction,
packs a full symmetric tensor, and combines full-matrix off-diagonal
cotangents. An independent test covers the closed-form oracle, central
differences, the adjoint identity, and invalid directions. Force divergence,
volume assembly, traction traces, Braginskii corrections, and field-aligned
assembly remain separate active work; these blocks are not a claim that the
full anisotropic MHD operator is complete.

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

FortNum must provide tested ordinary Legendre functions, associated and
half-integer Legendre functions, spherical harmonics, cylindrical Bessel and
Hankel functions, and toroidal harmonics with stable recurrences and
continuation rules. Analytical solutions and DtN maps are tested on slabs,
cylinders, spheres, and exact toroidal surfaces. The torus API must state
whether a function is a toroidal harmonic, a Fourier mode in toroidal
coordinates, or a numerical continuation of a special function. These are not
interchangeable names.

## 8. PDE and boundary-operator ingredients

The same residual and derivative contracts cover scalar and vector-valued
equations.

### 8.1 Core equations

- Poisson, diffusion-reaction, and Grad-Shafranov operators;
- scalar Helmholtz with FEM, BEM, DtN, and PML boundaries;
- Ampere, curl-curl, eddy-current, Maxwell, and anisotropic H(curl) forms;
- H(div) flux, mixed Poisson, magnetic induction, and divergence constraints;
- linear ideal and resistive MHD blocks, including singular layers and wall
  response as ingredient problems;
- reduced-MHD brackets, energy functionals, and compatible time operators.

Pressure and stress are represented as scalar, vector, or tensor fields as
required by the closure. A tensor-valued pressure is not silently projected to
its trace. Its symmetry, divergence, traction, and work pairing are tested.

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
an independent dot-product identity; map/coefficient/volume sensitivities
remain separate follow-up contracts. The public 1D linear interpolation-map
builder now checks partition of unity, affine reproduction, fixed-topology
JVP/VJP dot products, and Cartesian bilinear affine reproduction. Higher-order
or unstructured interpolation Jacobians, support-volume construction, and
anisotropy-aware preconditioning remain active work.

### 8.3 FEM/BEM, DtN, and PML

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

## 10. Plasma-oriented ingredient matrix

The following table defines parity targets. It does not authorize copying
source code or shipping a replacement for the named code.

| Reference target | Physics class | FortFEM ingredient parity target |
| --- | --- | --- |
| [CHEASE](https://crppwww.epfl.ch/~sauter/chease/), [paper](https://doi.org/10.1016/0010-4655(96)00046-X) | 2D fixed-boundary axisymmetric toroidal equilibrium | Grad-Shafranov residual, bicubic or spline/FEM geometry, axis regularity, profile and boundary constraints, COCOS conversion, convergence and G-EQDSK oracle |
| [FreeGS](https://freegs.readthedocs.io/en/stable/creating_equilibria.html) | 2D free-boundary Grad-Shafranov equilibrium | Coil and wall geometry, Picard/Newton residual, X/O-point constraints, free-boundary coupling, GEQDSK I/O, analytical Solovev and manufactured profiles |
| [VMEC/PARVMEC](https://github.com/ORNL-Fusion/PARVMEC), [VMEC++ numerics](https://arxiv.org/abs/2502.04374) | 3D nested-surface variational ideal equilibrium | Fourier angles plus radial FE/IGA, energy functional, flux constraints, axis treatment, free boundary, shape JVP/VJP, VMEC-format comparison |
| [GVEC](https://gvec.readthedocs.io/develop/index.html), [DESC](https://github.com/PlasmaControl/DESC) | Flexible 3D variational equilibrium and optimization | General coordinate maps, radial B-splines, Fourier modes, multiple interfaces, exact residual derivatives, continuation and optimization fixtures |
| [SPEC](https://princetonuniversity.github.io/SPEC/) | Multi-region relaxed MHD with ideal interfaces | Region graph, independent fields, Beltrami curl eigenproblem, flux/helicity constraints, total-pressure balance, interface shape derivatives, islands as topology fixtures |
| [GPEC](https://princetonuniversity.github.io/GPEC/), [references](https://princetonuniversity.github.io/GPEC/references.html) | Linear ideal, kinetic, and resistive perturbed response | Equilibrium import, Fourier coupling, singular outer/inner layer contracts, vacuum and wall response, response matrices, normalization and reciprocity tests |
| [MARS-F response work](https://doi.org/10.1016/j.cpc.2006.09.003) | Linear toroidal ideal/resistive MHD and wall response | Linearized residual blocks, complex frequency, plasma-vacuum-wall coupling, Fourier-FEM assembly, resistive layer matching, benchmark outputs without MARS source |
| [GLISS](https://github.com/itpplasma/GLISS) | Global linear ideal-MHD stability in 3D toroidal equilibria | Energy-principle residual, compatible radial spline FE, Fourier mode topology, GVEC/VMEC input adapters, eigenvalue and inertia oracles, Enzyme-compatible kernel boundaries |
| [JOREK](https://www.jorek.eu/), [overview paper](https://arxiv.org/abs/2011.09120) | Nonlinear extended and resistive MHD | 2D FE plus toroidal Fourier blocks, coupled residuals, anisotropic transport, implicit structure-aware stepping, wall and free-boundary traces, operator-level parity tests |
| MEPHIT and STARWALL | Electromagnetic response and resistive-wall coupling | H(curl)/H(div) FEEC, surface traces, FEM/BEM/DtN wall blocks, retained complex factors, interface currents, reciprocity and energy tests |
| [BOUT++](https://bout-dev.readthedocs.io/en/stable/user_docs/introduction.html) | General 3D curvilinear plasma fluid framework, with reduced edge models and implicit time integration | Equation-as-data residuals, curvilinear metric and boundary contracts, field-aligned operators, mixed conservative fluxes, model-level JVP/VJP, and a small transport benchmark |
| [GRILLIX](https://physik.uni-greifswald.de/ag-manz/forschung/codes/grillix/), [FCI paper](https://doi.org/10.1088/1361-6587/aaa373) | 3D edge and scrape-off-layer turbulence using flux-coordinate-independent operators and drift-reduced Braginskii physics | FCI field-line tracing and interpolation, parallel/perpendicular operator split, immersed boundaries, anisotropic transport, sheath interfaces, and manufactured MMS fixtures |
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

### Plasma and MHD ingredients

14. **Cylindrical Grad-Shafranov.** Solovev and CHEASE-style profiles, axis
    regularity, Fourier-FEM, and FreeGS/CHEASE comparison data.
15. **Fourier-FEM slab and cylinder.** Mode diagonal operators, retained
    nonlinear triads, real/conjugate packing, and torus-harmonic diagnostics.
16. **Multi-region Beltrami equilibrium.** SPEC-like regions, ideal interfaces,
    flux/helicity constraints, and a pressure-balance residual.
17. **Linear 3D perturbed equilibrium.** GPEC/MARS-like mode response with
    vacuum, wall, singular-layer, and response-matrix toy operators.
18. **GLISS energy-principle toy spectrum.** Radial spline FE, Fourier modes,
    inertia count, eigenpair derivatives, and GVEC/VMEC-shaped input.
19. **HKT or resonant shielding sheet.** Ideal current-sheet limit, finite
    resistive layer, XFEM enrichment, fitted sheet, and convergence to the
    analytical singular solution.
20. **Reduced-MHD bracket.** Energy-skew nonlinear bracket, Fourier convolution,
    analytical JVP, and long-time structure-preserving integration.
21. **Resistive diffusion and tearing layer.** Explicit layer, adaptive layer,
    asymptotic enrichment, DG, and ideal-limit comparison.
22. **Ingredient-only JOREK path.** Coupled magnetic flux, vorticity, density,
    temperature, parallel velocity, and electric-potential residual stubs as
    separately testable operators. This is not a JOREK implementation.

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
26. **Tensor-pressure plasma wave.** An anisotropic pressure tensor with
    parallel, perpendicular, and gyrotropic parts, including force balance and
    energy diagnostics.
The current `cgl_pressure_tensor` gallery example provides the first
manufactured constitutive/force-divergence profile and CSV/1D FortPlot
outputs; the coupled wave and higher-dimensional cases remain active.
The `field_aligned_flux` gallery example now provides a generated
parallel/perpendicular profile and (k_\parallel/k_\perp=100) flux plot; a
full assembled anisotropic diffusion gallery case remains active.
27. **Field-aligned diffusion.** A slab, cylinder, and torus with extreme
    \(k_\parallel/k_\perp\), comparing aligned coordinates, FCI maps, Fourier-
    FEM, and an isotropic control case.
28. **Edge/SOL ingredient model.** A small drift-reduced Braginskii or heat
    conduction system with parallel transport, sheath or wall traces, and a
    reproducible FCI field-line map. This is a numerical ingredient example,
    not a GRILLIX or BOUT++ replacement.

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
- Keep FortSym revision pins and generated-kernel checks green.

### Phase 1: Interface calculus: **active**

- The public orientation-safe scalar/vector trace contract now provides
  plus/minus averages and jumps, normal/tangential projections, and the
  rotated Ampere surface-current jump with an independent sign oracle.
- Internal manifolds, level sets, and surface measures remain next.
- Broken H1, H(curl), H(div), and L2 spaces plus skeleton spaces.
- Explicit delta-source and surface-current weak terms.
- Fitted duplicated spaces, Nitsche, mortar, multipliers, and block constraints.

### Phase 2: Cut geometry and XFEM/XIGA: **planned**

- Cut-cell classification and high-order quadrature.
- Heaviside, kink, singular, helical, and resonant enrichments.
- Shifted bases, corrected blending elements, pruning, conditioning, and
  connected-component activation.
- Trimmed B-spline stabilization and cut H(curl)/H(div) extensions.

### Phase 3: DG and HDG: **planned**

- Scalar SIPG and non-symmetric variants.
- Conservative/upwind/entropy-aware flux interface.
- Broken vector FEEC, hybridization, static condensation, and mixed CG-DG.

### Phase 4: Open boundaries and vector operators: **active**

- Complete scalar and Maxwell FEM/BEM, DtN, and PML parity on slab, circle,
  sphere, cylinder, and torus.
- Add exact-curved-torus external surfaces and vector DtN.
- Add larger-domain, BEM, DtN, and PML comparisons for Poisson, Ampere, and
  Helmholtz.
- Verify FortPlot mesh and surface rendering for every mesh-bearing plot.

### Phase 5: Fourier-FEM and torus harmonics: **active**

- Stabilize ordinary, associated, and half-integer Legendre and toroidal
  harmonic APIs in FortNum.
- Define mode normalization, phase, field-period, and real packing.
- Add mode-coupled scalar, H(curl), H(div), and reduced-MHD operators.

### Phase 6: Structure-preserving time evolution: **active**

- The public `advance_mixed_wave_midpoint` step now provides the common
  first-order pressure/velocity, displacement/momentum, and port-Hamiltonian
  Cayley contract. Its independent test checks the oscillator map, energy, and
  signed-step reversibility.
- Variational/symplectic and Poisson building blocks for ideal terms.
- Energy-dissipative integrators for resistive and viscous terms.
- Symmetric splitting, implicit midpoint/Cayley, discrete-gradient or
  average-vector-field options, and long-time invariant tests.
- Mixed first-order pressure-velocity, displacement-velocity-stress, and
  electromagnetic wave states with a common port-Hamiltonian interface.
- Tensor-valued pressure and anisotropic constitutive blocks with exact power,
  momentum, and stress-work diagnostics.
- The generated CGL pressure-tensor and product-rule force-divergence blocks
  are public and tested; volume assembly, traction, and stress-work coupling
  remain active follow-up slices.

### Phase 7: Equilibrium and linear-response ingredients: **planned**

- CHEASE and FreeGS-style 2D Grad-Shafranov fixtures.
- VMEC/GVEC/DESC-style Fourier-variational 3D equilibrium fixtures.
- SPEC-style multi-region Beltrami and interface constraints.
- GPEC/MARS/GLISS-style linear ideal and resistive perturbation blocks,
  singular layers, vacuum, conducting wall, and response matrices.

### Phase 7a: Field-aligned edge and SOL ingredients: **active**

- FCI RK4 field-line tracing, support-operator gradient/divergence algebra, and
  matrix-free (P K_\parallel Q) diffusion (including FortSym-generated
  gradient JVP/VJP products) are public and independently tested. The
  `fci_parallel_diffusion` gallery fixture provides a reproducible open-line
  mass-conservation and negative-energy profile. A 1D linear map builder and
  its fixed-topology JVP/VJP provide independent partition/affine and
  dot-product oracles, and a 2D Cartesian bilinear builder now covers a
  genuine poloidal slice. The RK4 tracer also has a tangent callback path with
  an exponential endpoint oracle. The positive staggered flux-box volume constructor
  now covers the traced expansion/area/(B_\varphi) product with pinned
  FortSym-generated value/JVP/VJP kernels. Higher-order or unstructured
  interpolation, while fixed-cell 2D map JVP/VJP products now cover source
  and target coordinate motion. Curved support-volume measures remain a
  separate planned component.
- FCI field-line maps, higher-dimensional interpolation Jacobians, and parallel
  derivative JVPs.
- Strongly anisotropic diffusion, conduction, resistivity, and viscosity.
- Immersed target plates, sheath or wall traces, and open-field-line boundary
  conditions.
- Small MMS and performance cases aligned with PARALLAX, GRILLIX, and BOUT++
  concepts, without copying their implementations.

### Phase 8: Cross-code oracles and gallery: **active**

- Ordered examples with FortPlot plots, numerical data, convergence, and
  performance.
- Optional lightweight FEniCSx, FreeFEM, and MFEM runners.
- Sister-repository data for heavy or licensed references.
- GitHub Pages generation, link checks, and periodic deployment monitoring.

### Phase 9: Future application layer: **reference only**

Use the stabilized ingredients to support future research applications that
may resemble GVEC, VMEC, SPEC, GPEC, MARS, GLISS, or JOREK. Do not start a
production code in this repository until the corresponding ingredient,
analytical example, and external oracle have passed the previous phases.

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
- structure-preserving properties are tested rather than asserted;
- focused `fo` tests meet the short feedback target;
- full CI passes on supported compilers;
- the example and documentation are generated reproducibly;
- transient output is ignored and no unrelated downstream fix is duplicated.

## 16. Non-goals and explicit boundaries

- No full VMEC, GVEC, SPEC, GPEC, MARS, GLISS, JOREK, CHEASE, or FreeGS
  replacement is planned in FortFEM.
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
