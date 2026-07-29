# FortFEM capability audit and implementation roadmap

Date: 2026-07-29

## Scope and evidence

This audit answers four questions:

1. Which parts of FortFEM can use `fortnum` and `fortsparse`?
2. Can FortFEM replace the finite-element work currently performed for MEPHIT?
3. Does FortFEM cover the methods used by the magnetic and acoustic papers?
4. What implementation and validation sequence leads to arbitrary-order
   H(curl), H(div), DtN, BEM, and FEM-BEM support?

The findings use the following evidence:

- FortFEM `main` at `58dda23`, including a clean full `fo` run with 48 passing
  test targets.
- MEPHIT `main` at `fd2c56b` and its FreeFem++ Maxwell daemon.
- The published magnetic-paper source for
  [CPC 277 (2022) 108401](https://doi.org/10.1016/j.cpc.2022.108401), the
  current `paper_magnetic` repository, and pull request 9.
- `paper_acoustics` `main` at `6300ab0`, including its Fortran FEM+DtN
  prototype and FreeFem++ reference implementation.
- The archived ISNVH 2018 paper, FreeFem++ models, FEM/BEM books, and
  H(div)/H(curl) literature under `/mnt/files/archive`.
- The work brain notes for MEPHIT, finite elements, H(div)/H(curl), the
  magnetic paper, and computational acoustics.
- The finite-element definitions and implementation guidance in
  [DefElement](https://defelement.org/),
  [MFEM](https://mfem.org/features/), and the finite element exterior
  calculus literature
  ([Arnold, Falk, and Winther, 2010](https://arxiv.org/abs/0906.4325)).
- The exact nonreflecting-boundary work of
  [Grote and Keller](https://doi.org/10.1006/jcph.1995.1210) and the
  [Bempp boundary-operator documentation](https://bempp.com/handbook/api/boundary_operators.html).

## Repository state

At the initial audit, FortFEM contained 22,622 lines of Fortran across `src`,
`test`, and `example`, of which tests accounted for 10,442 lines. The full
build, test, and lint pipeline passed. The compiler reported array-temporary
warnings in several numerical and plotting paths.

The implemented scalar path covers P1 and P2 triangles, Q1 quadrilaterals,
Poisson assembly, Dirichlet and Neumann conditions, dense LAPACK solves,
iterative solvers, and plotting. The native constrained Delaunay mesher has
substantial MEPHIT-oriented tests. Triangle-compatible output and a CMake
consumer interface remain on unmerged pull request 56.

The expression API records descriptions such as `inner(grad(u), grad(v))`.
It does not compile a general form into element kernels. The high-level
solvers recognize expression strings and dispatch to problem-specific dense
assemblers.

### Vector finite elements

The initial audit found a linearly dependent edge basis, the wrong Piola map,
incorrect curls, no distinct Raviart--Thomas basis, and no global orientation
map. Those reference-element and assembly defects are now repaired for
lowest-order triangles and covered by analytical tests.

The remaining high-level vector path is not a verified PDE solver:

- `vector_function_space(..., "RT", 1)` still shares the generic edge-space
  container and is not connected to RT-specific interpolation or mixed forms.
- The curl-curl assembler does not evaluate the basis module. It assigns
  element-local indices `3*(triangle_id-1)+i` to a space whose dimension is
  the number of unique mesh edges.
- The mass and load terms are placeholders. Boundary values other than zero
  are ignored.
- Existing vector tests establish that calls complete. They do not compare
  edge moments, assembled matrices, fields, or convergence against an
  independent solution.

FortFEM now has verified lowest-order H(curl) and H(div) kernels and global
operators. Arbitrary polynomial order, tetrahedra, mixed spaces, orientation
transforms for higher moments, a discrete de Rham complex, and integration
with the expression API are absent.

### Sparse algebra

FortFEM owns a real CSR matrix, CSR-to-CSC conversion, matrix-vector product,
and preconditioners. Direct real and complex solves now use `fortsparse` with
the permissive in-process SuperLU backend. Assembly still starts from dense
global matrices in the main PDE solvers.

## Consumer requirements

### MEPHIT

MEPHIT currently uses lowest-order elements. Its Maxwell daemon defines:

- P1 scalar fields;
- `RT0Ortho`, a lowest-order H(curl) space;
- `RT0`, a lowest-order H(div) space;
- oriented global edge degrees of freedom with an explicit MEPHIT-to-FreeFem
  sign map;
- the Fourier-harmonic bilinear form
  \[
  \int_\Omega R\,\operatorname{curl} w\,\operatorname{curl} A
  + \frac{n^2}{R}w\cdot A\,dR\,dZ;
  \]
- essential tangential data on boundary label 2;
- repeated solves with one factorization and separate real and imaginary
  right-hand sides.

MEPHIT also exchanges RT0 and P1 coefficients through a C/pipe interface,
uses Triangle-compatible mesh numbering, evaluates RT0 fields and their
derivatives, computes weighted edge projections and L2 norms, and reconstructs
the toroidal component from the Fourier divergence constraint.

FortFEM cannot replace this daemon yet. The minimum replacement gate is a
verified lowest-order H(curl)/H(div) pair, global orientations matching the
MEPHIT mesh, the weighted Fourier form, factor reuse, coefficient
interpolation, and a compatibility test against MEPHIT's six mesh fixtures and
FreeFem++ matrices. Current MEPHIT inputs use only lowest order. The broader
library goal requires arbitrary order.

### Magnetic paper

The 2022 paper uses scalar and edge spaces for decoupled Fourier harmonics in
general curvilinear coordinates. The transverse problem is a weighted
curl-curl Helmholtz problem. Its coefficients include metric density,
inhomogeneous or anisotropic reluctivity, and the Fourier mode. The paper
treats Dirichlet and Neumann data and validates cylindrical, shielding, and
racetrack-coil cases. Current work in `paper_magnetic` adds a three-dimensional
lowest-order Nedelec box problem with an analytical center value and magnetic
field comparison.

FortFEM lacks the verified edge space, tensor coefficients, Fourier form,
general boundary terms, and three-dimensional mesh and element support.
It cannot reproduce the paper at present.

### Acoustic and ultrasonic papers

The acoustic work uses three related exterior treatments:

- a planar half-space DtN operator diagonalized by a Fourier transform;
- Kirchhoff-Helmholtz propagation using the two-dimensional outgoing Green's
  function;
- iterative FEM/BEM-style pressure and normal-velocity coupling.

`paper_acoustics` already contains a Fortran prototype. It assembles
complex-valued P1 elasticity, factors with complex UMFPACK, applies a periodic
planar DtN multiplier with a hand-written DFT, and solves the interface
equation with restarted GMRES. Its committed tests compare the Fortran result
to FreeFem++ and include plane-wave and FDTD validation.

This code is a useful migration source, but it is not a general BEM. The
interface is planar and periodic, the exterior medium is homogeneous, the
transform is quadratic-cost, and no single-layer, double-layer, adjoint
double-layer, or hypersingular boundary operator is assembled. FortFEM
contains none of this functionality.

## Use of fortnum and fortsparse

### fortnum

The following replacements are aligned with current consumers:

| fortnum facility | Replacement or new use in FortFEM |
| --- | --- |
| `fft_c2c` and reusable FFT plans | Replace the quadratic DFT in the planar DtN operator. |
| `gauss_legendre` and `gauss_legendre_ab` | Edge moments, boundary quadrature, and MEPHIT's current GSL fixed Gauss rule. |
| Gauss-Kronrod and CQUAD integration | Verification integrals and near-singular one-dimensional boundary integrals. Singular self-panels still need transformed or analytical rules. |
| complex Bessel functions | Circular and spherical outgoing DtN eigenvalues after a tested Hankel wrapper is added. |
| `det2`, `det3`, `inv2`, and `inv3` | Shared affine and isoparametric mapping kernels after the element API stabilizes. |
| interpolation routines | Boundary data transfer where the finite-element trace and boundary-element mesh differ. |

The existing triangle quadrature cannot be replaced by fortnum today because
fortnum supplies one-dimensional rules. Triangle and tetrahedron rules need a
simplex-rule module or a Duffy tensor-product construction.

FortFEM now uses `fft_c2c` for a standalone periodic planar Helmholtz DtN map.
Exact positive and negative Fourier modes test the outgoing propagating branch
and the decaying evanescent branch at a prime transform length.

### fortsparse

`fortsparse` already supplies the required real and complex CSC types,
duplicate-compressing triplet construction, factor reuse, and real or complex
solves. SuperLU is the permissive in-process default. UMFPACK is available
through a separate helper process.

FortFEM now uses it for direct real and complex sparse solves, with SuperLU as
the in-process backend. The old UMFPACK C binding and direct SuiteSparse link
set have been deleted. The migration retains the internal CSR type only for
existing iterative preconditioners. Assembly should next emit triplets
directly into `fortsparse` CSC storage. The final solver layer should accept an
abstract matrix-vector operator so sparse FEM blocks and dense or fast BEM
blocks can share Krylov methods without forming a monolithic dense matrix.

## Target architecture

### Reference elements and mappings

Reference-element definitions will use polynomial moments, not tables of
pre-expanded basis coefficients. For each family and degree, FortFEM will
construct the moment matrix in a scaled orthogonal polynomial basis and solve
for the nodal basis once.

The initial two-dimensional sequence on triangles is:

\[
P_{k+1}
\xrightarrow{\operatorname{grad}}
\mathcal P^-_{k+1}\Lambda^1
\xrightarrow{\operatorname{curl}}
P_k^{\mathrm{disc}}.
\]

It will include:

- arbitrary-order Lagrange;
- Nedelec first and second kind;
- Raviart-Thomas;
- Brezzi-Douglas-Marini;
- discontinuous scalar spaces;
- covariant and contravariant Piola maps;
- edge-moment sign and permutation transforms;
- interpolation, trace evaluation, curl, and divergence;
- curved isoparametric geometry after affine elements pass.

The three-dimensional sequence on tetrahedra follows after the two-dimensional
API and tests stabilize. Hexahedra and quadrilaterals use tensor-product
families in a later phase.

### Assembly

Forms will carry executable coefficients and operators. Element assembly will
use:

- local-to-global maps derived from mesh entities;
- explicit orientation transforms;
- quadrature order selected from trial degree, test degree, geometry degree,
  and coefficient degree;
- scalar, vector, tensor, real, and complex coefficients;
- cell, boundary, and interface integrals;
- sparse triplet accumulation.

Problem-specific helpers may remain, but they must call the same element
kernels as the form path.

### Exterior operators

Two exact DtN families precede general BEM:

1. A planar periodic half-space operator for the acoustic paper. It handles
   propagating and evanescent Fourier modes and uses `fortnum` FFT plans.
2. Circular and spherical artificial boundaries. Their modal eigenvalues use
   outgoing Hankel functions and support explicit truncation control.

The general BEM starts with two-dimensional Laplace and Helmholtz kernels on
oriented line panels. It provides the four Calderon operators:

- single layer \(V\);
- double layer \(K\);
- adjoint double layer \(K'\);
- hypersingular operator \(W\).

The dense reference assembler uses analytical logarithmic self-panel terms,
singularity subtraction or Duffy transforms for adjacent panels, and adaptive
near-singular quadrature. Hypersingular terms use a weakly singular
regularization. Constant and continuous linear trace spaces are the first
discretizations. Curved panels and higher-order traces follow.

Exterior Helmholtz solves use a combined-field formulation to avoid interior
resonances. FEM/BEM transmission uses Costabel's symmetric coupling as the
reference formulation. Johnson-Nedelec remains an optional lower-cost
formulation for comparison. Large problems use matrix-free FMM or hierarchical
compression only after the dense operators pass the same tests. The fast path
must preserve the dense API and expose an operator application to GMRES.

Three-dimensional Laplace, Helmholtz, and Maxwell boundary operators on
surface triangles form a separate phase. They must not delay the verified 2D
acoustic and MEPHIT paths.

## Validation matrix

Repository-state checks do not count as tests. Every test below has an
independent mathematical or cross-code oracle.

| Area | Required oracle |
| --- | --- |
| Nedelec and Raviart-Thomas basis | Kronecker edge moments, exact polynomial reproduction, known reference-cell curls/divergences, and continuity across a reflected affine cell. |
| Piola maps | Exact tangential or normal trace preservation under non-axis-aligned affine maps. |
| FEEC | Discrete `curl(grad)=0` or `div(curl)=0`, commuting interpolation, exact rank/nullity, and manufactured convergence rates. |
| H(curl) solve | Manufactured smooth field, Maxwell eigenmodes, and the `paper_magnetic` box series value and magnetic field. |
| H(div) solve | Manufactured mixed Poisson/Darcy flux with local conservation and optimal flux/divergence convergence. |
| MEPHIT | FreeFem++ matrix and solution comparison, coefficient round-trip with the edge sign map, and six Triangle mesh fixtures. |
| Planar DtN | Individual propagating and evanescent Fourier modes, zero-reflection plane wave, and the committed acoustic paper fixture. |
| Circular DtN | Exact outgoing cylindrical modes and truncation convergence. |
| BEM operators | Analytical circle spectra, jump relations, Calderon identities, and boundary traces of off-surface point sources. |
| Exterior Helmholtz BEM | Mie-series scattering by a circle across frequencies that include interior resonances. |
| FEM/BEM | Manufactured transmission on a disk, then the same layer-potential solution on a non-circular polygon or spline boundary. |
| Fast BEM | Dense operator applications as the oracle over a size range, followed by asymptotic memory and runtime checks. |

Convergence tests require at least three refinements and compare the measured
slope to the theoretical order. A single coarse-grid tolerance is a smoke
test, not a convergence test.

## Execution sequence

### Implementation status on 2026-07-29

The Phase 0 foundations and the first Phase 1 interfaces are implemented:

- the triangular lowest-order Nedelec basis, curl, and covariant Piola map
  pass exact edge-moment and affine-map tests;
- a distinct RT0 basis, divergence, and contravariant Piola map pass exact
  normal-flux tests;
- the mesh exposes local-to-global edge signs, and the global Nedelec
  curl-curl plus mass operator reproduces the exact continuum energy of
  constant fields;
- RT0 div-div plus mass assembly reproduces its exact reference matrix and
  the continuum energy of globally conforming constant fluxes;
- `fortsparse` replaces FortFEM's direct UMFPACK binding and is tested with
  exact real and complex systems; a retained-factor interface solves multiple
  right-hand sides without refactorization.
- a periodic planar Helmholtz DtN map uses `fortnum` FFTs and passes exact
  propagating and evanescent Fourier-mode tests.
- coefficient callbacks and a built-in axisymmetric Fourier kernel assemble
  the MEPHIT weights \(R\) and \(n^2/R\); the latter converges to an exact
  logarithmic energy on a shifted annular strip.
- oriented real and complex Nedelec tangential and RT0 normal edge moments use
  `fortnum` Gauss-Legendre rules and pass exact polynomial line-integral
  tests; the RT path also exposes MEPHIT's \(R\,v\cdot n\) cylindrical
  projection directly.
- complex RT0 coefficients reconstruct physical fields and element
  divergences, provide an L2 norm, and recover the toroidal Fourier component
  required by the zero-divergence constraint.
- the first C ABI entry accepts zero-based interleaved triangle meshes and
  returns globally oriented edge degrees of freedom and local signs; a
  persistent mesh-handle variant retains that topology for subsequent field
  operations.
- the C ABI also accepts zero-based complex CSC matrices, retains
  `fortsparse` factors behind opaque process-local handles, reuses them across
  right-hand sides, and rejects released handles.
- complex RT0 coefficients cross the C boundary in global edge-DOF order for
  point evaluation, divergence, L2 integration, and toroidal Fourier
  reconstruction.
- an installable `fortfem::capi` shared CMake target builds the focused
  MEPHIT-facing surface with pinned `fortsparse`; a standalone C consumer
  finds the installed package, links, and exercises a retained mesh.
- MEPHIT `main` optionally finds `fortfem::capi`, retains its production
  triangle mesh through the C ABI, maps MEPHIT edge order and direction to
  FortFEM's global RT0 degrees of freedom, and uses FortFEM for the RT0 L2
  norm. An independent affine-field test covers a clockwise cell and
  reordered, reversed edges. The default build and the remaining FreeFem++
  operations are unchanged.

The Phase 0 exit gate is met at the numerical-kernel level. Connection to
executable weak forms, DtN integration into scalar and elastic boundary
forms, full MEPHIT operator replacement, and consumer validation remain.

### Phase 0: correct claims and numerical dependencies

1. Add analytical tests for lowest-order triangular Nedelec and
   Raviart-Thomas elements.
2. Replace the broken edge basis and Piola map.
3. Replace local triangle-edge indices with global oriented entity maps.
4. Migrate direct real and complex factorization to `fortsparse`.
5. Port the planar DtN operator from `paper_acoustics` and replace its DFT
   with `fortnum` FFT.
6. Change README feature claims to distinguish verified, experimental, and
   planned capabilities.

Exit gate: lowest-order H(curl) and H(div) interpolation and assembly pass
analytical tests. The planar DtN map passes propagating and evanescent mode
tests. No direct UMFPACK binding remains in FortFEM.

### Phase 1: MEPHIT replacement

1. Implement the weighted Fourier curl-curl and mass form.
2. Add RT0 projection, interpolation, L2 integration, and Fourier divergence
   reconstruction.
3. Merge or reproduce the Triangle-compatible CMake consumer interface.
4. Add a C ABI for mesh, factor, solve, and coefficient transfer.
5. Compare matrices and solutions against the FreeFem++ daemon.
6. Integrate FortFEM on a MEPHIT branch and run the six mesh cases plus one
   full 33353 reference case.

Exit gate: MEPHIT can run without FreeFem++ or Triangle for the validated case,
with field and residual tolerances recorded from independent reference data.

### Phase 2: arbitrary-order FEEC

1. Implement moment-generated arbitrary-order triangle families.
2. Add orientation transformations and mixed-space degree-of-freedom maps.
3. Add arbitrary-order quadrature and sparse assembly.
4. Prove behavior through the FEEC validation matrix.
5. Add tetrahedral families and the magnetic-paper box case.

Exit gate: optimal convergence holds for H1, H(curl), H(div), and L2 members
of the discrete sequence for degrees 0 through 4. The 3D box test reproduces
its analytical vector potential and curl.

### Phase 3: acoustic DtN and 2D BEM

1. Integrate the planar FFT DtN operator into scalar and elastic boundary
   forms.
2. Add circular outgoing DtN with truncation diagnostics.
3. Implement dense Laplace and Helmholtz Calderon operators.
4. Add combined-field exterior solves and field evaluation.
5. Implement symmetric FEM/BEM transmission coupling.
6. Reproduce the acoustic paper's plane-wave and layered validation cases.

Exit gate: the analytical DtN, circle-scattering, Calderon, and manufactured
non-circular FEM/BEM tests pass. The paper fixture agrees within its published
discretization error.

### Phase 4: fast and three-dimensional BEM

1. Add FMM or hierarchical compression behind the boundary-operator API.
2. Add adaptive boundary refinement and residual estimators.
3. Add curved high-order panels.
4. Add 3D Laplace and Helmholtz surface elements.
5. Add Maxwell traces and electromagnetic FEM/BEM coupling.

Exit gate: fast and dense operators agree to requested tolerance. Memory grows
subquadratically on the benchmark sequence. Three-dimensional sphere
scattering and magnetostatic transmission match analytical series.

## Commit policy

Each commit contains one test-driven capability or one dependency migration.
The red test is run before implementation. The focused test, affected tests,
and bare `fo` pipeline are run before every commit. Commits go directly to
`main` as requested and are pushed after each passing slice. Cross-repository
integration changes remain in their consumer repository and use FortFEM
commits that already pass the local pipeline.
