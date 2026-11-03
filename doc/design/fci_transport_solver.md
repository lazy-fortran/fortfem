---
title: FCI anisotropic transport solver
---

# FCI anisotropic transport solver

The solver is blocked at its first accuracy gate: the native bilinear
flux-coordinate-independent (FCI) map loses the field-line-constant limit
on a non-aligned nested-surface fixture. Even the dense oracle gives nearly
unit relative solution error while satisfying the discrete equations.
Option B, a two-level Krylov preconditioner, cannot correct this spatial
error for the tested operator. No production solver is delivered or qualified.
Island, plate-intersection, single-null, Krylov and cost gates remain open.
Two reduced multiplier candidates also fail the limit-space preflight below:
both retain antisymmetric checkerboards. Neither is selected for production.

The dense kernel/complement reference is test support in
`test/fci_transport_reference.f90`, linked only into its test consumers.
It is neither installed nor part of the public library API. Its original
algebraic verification is retained below.

This is a first part of the transport work required by
[KIN6D's nonlinear resistive MHD design](https://github.com/itpplasma/kin6d/blob/main/docs/nonlinear-resistive-mhd.md),
particularly its anisotropic transport and preconditioner gates RM-B3/RM-B4.
The existing [FCI support operator](fci_parallel_operator.md) supplies the
parallel gradient; no new diffusion discretization or physical closure is
introduced here.

## Reduced-space preflight

On the sheared fixture defined below, `test_fci_transport_reduced_space`
evaluates two spaces for testing the mapped parallel difference. Let `h`
be the Cartesian section spacing. `P0-2h` contains
independent constants on disjoint 2-by-2 blocks of interior support samples
in each plane interval. `Q1-2h` contains continuous bilinear functions on a
coarse section mesh with spacing `2h`; its last interval has width `h` to
include the final support sample. The native bilinear builder supplies its
basis values. These are discrete L2 projections using the existing diagonal
support pairing, rather than exact integration over continuous elements.
Boundary values remain zero on the physical perimeter; periodic planes are
identified exactly as in the original fixture.

Let `G` be the original interval gradient after eliminating wall nodes, `H`
the chosen test basis evaluated at support samples, and
`L L^T = H^T H`. The normalized constraint is `C = L^(-1) H^T G`.
The uniform nodal and support masses cancel in this expression. The two
identical intervals on two periodic planes split orthogonally into
plane-symmetric and plane-antisymmetric blocks. Both blocks are retained.
No connectivity is supplied independently of the native traced maps.

The dense singular value decomposition of `C` reports its rank and kernel.
Its smallest singular value, including a zero if the multiplier has a
cokernel, is `beta_l2`. The second value, `beta_h1`, uses the discrete
H1 norm of the temperature `u`: mass plus the Cartesian five-point section
energy and the periodic third-direction energy. An orthonormal discrete sine
transform diagonalizes this norm; an independent highest-sine-mode eigenvalue checks
the scaling. The multiplier norm remains the discrete L2 norm. Neither
value uses the mesh-dependent dual norm of the stabilized AP analysis.
The rank threshold is the existing oracle's
`8 * machine_epsilon * max(rows,columns) * largest_sigma`.

The meshes have 21, 41 and 81 cells per section direction. The reference
surface-resolution count is 10, 20 and 40 radial bins from the wall to the
axis, at Cartesian spacing `h`; it is a resolution proxy, not a theorem
counting all representable flux functions. Kernel dimension grows
quadratically with section resolution for both candidates.

For an independent physical discriminator, the initial checkerboard has
nodal values `(-1)^(i+j) sin(pi*x_i) sin(pi*y_j)` and opposite signs on the
two planes, with interior grid indices `i` and `j`. It is normalized before
orthogonal projection onto the reduced kernel. Every sampled flux function
has identical values on both planes, so the projected antisymmetric vector
is exactly orthogonal to that whole
physical subspace. The retained vector also has an order-one defect under
the original, unprojected gradient. On the coarse mesh the existing dense
transport oracle independently confirms the full nullity and reproduces
the projected checkerboard using an identity slow operator and parallel
coefficient `1e13`. This is an algebraic mass-projection check, not a
manufactured transport accuracy qualification.

The [complete measurements](figures/fci_transport_reduced_space.csv) retain
all six rows, both inf-sup norms, mass-projection errors of the exact flux,
constraint defects and diagnostic timings. On the finest mesh, `P0-2h` and
`Q1-2h` have kernel dimensions 9,600 and 9,438, retain 98.91% and 98.51% of
the checkerboard, and have H1 inf-sup values `7.643e-4` and `3.641e-4`.
Those values decrease on both refinements; a uniform lower bound is not
established. Reduced constraint defects remain below `3.7e-15`, whereas
the unprojected gradient defects are between 1.73 and 1.87. These quantities
are dimensionless. The full-rank multiplier maps leave respectively 75%
and approximately 74% of the finest temperature space unconstrained.

The [plot source](figures/plot_fci_transport_reduced_space.py) shows the
inf-sup values, kernel growth and the actual coarse kernel vectors retained
in [P0 data](figures/fci_transport_P0-2h_kernel.csv) and
[Q1 data](figures/fci_transport_Q1-2h_kernel.csv). The latter use the same
grayscale amplitude scale and reverse sign on the second plane. They are
vectors computed by the dense oracle, not illustrative patterns.

With the Release/OpenBLAS configuration below:

```sh
cmake --build build/u3c -j 8 --target test_fci_transport_reduced_space
ctest --test-dir build/u3c -j 8 -R '^test_fci_transport_reduced_space_ladder$' -V
ctest --test-dir build/u3c -j 8 -R '^test_fci_transport_reduced_space_gate$' -V
python3 doc/design/figures/plot_fci_transport_reduced_space.py --output /tmp/fci-reduced-space.png
```

The ladder returns 93/93 diagnostic checks in 444.32 s on shared faepmac1,
gfortran 16.2, eight OpenMP threads and one BLAS thread. Timings cover dense
diagnostics and are not a solver cost comparison. The coarse diagnostic
returns 39/39; its separate acceptance invocation returns 39/41, with one
checkerboard rejection for each candidate.

The step-1 absence-of-spurious-kernel requirement fails already on the
coarse mesh. The registered `test_fci_transport_reduced_space_gate` enforces
zero checkerboard retention to roundoff and fails; it is neither disabled
nor marked as an expected failure. The diagnostic and refinement tests
check the measurements, and their passes do not qualify either space.
Step 2 accuracy, constrained-limit and inf-sup qualification on the four
geometries, and step 3 Krylov/cost gates, remain unopened. In particular,
the spatial rates required for a selected method have not been measured.

The failure is in the reduction itself:
`test/test_fci_transport_reduced_space.f90`, `test_basis`, drops about
three quarters of the interval constraints. Full row rank then leaves a
large antisymmetric kernel; `probe` constructs an explicit member using
the SVD row-space projector. Merely making the multiplier smaller is
insufficient. A candidate with the correct limit space and a controlled
inf-sup pairing is still required before solver work.

## First failing gate: nested surfaces with shear

The dimensionless domain has Cartesian coordinates `(x,y,z)` with section
`[0,1]^2` and a periodic third coordinate of period one.
Define the flux function `psi = sin(pi*x)*sin(pi*y)` and
`B = (a*psi_y, -a*psi_x, 1)`, where `a = 0.15*(1+psi)` supplies shear.
Here subscripts denote partial derivatives and `b = B/|B|` is the unit field.
The field is divergence free, tangent to nested `psi` surfaces, and nowhere
zero. Set perpendicular conductivity to one and let `epsilon` be its ratio
to parallel conductivity. With `I` the identity tensor and `b*b^T` the outer
product, the conductivity is `K = I + (1/epsilon-1)*b*b^T`.
The exact temperature solution of
`-div(K*grad(u)) = 2*pi^2*psi`, with zero outer Dirichlet data, is `u = psi`
for every positive epsilon. This follows from `b.grad(psi) = 0` and
`-Delta(psi) = 2*pi^2*psi`.

This is a translated, sheared, guide-field extension of the source and
temperature in Sovinec et al., [JCP 195 (2004), 355–386,
Section 4.2, pp. 374–375](https://doi.org/10.1016/j.jcp.2003.10.004)
([author manuscript](https://w3.pppl.gov/cemm/Project/NIMROD_JCP_rev1_wfigs.pdf),
pp. 35–36): their unnumbered benchmark equations specify
`Q = 2*pi^2*cos(pi*x)*cos(pi*y)` and
`T = cos(pi*x)*cos(pi*y)/chi_perp` on `[-0.5,0.5]^2`.
It does not reproduce their evolving magnetic field or finite-element tables.
Günter et al., [JCP 209 (2005), 354–370,
equations (1)–(2) and Section 2](https://doi.org/10.1016/j.jcp.2005.03.021),
use the same field-constant manufactured-source principle. Their circular
closed-line and `sin(2*pi*x)*cos(2*pi*y)` benchmarks have not been run here.

`test_fci_transport_nested` traces interior Cartesian nodes with the native
classical fourth-order Runge–Kutta field-line tracer, obtains endpoint
interpolation with the native
bilinear map builder, and calls the native parallel-gradient assembler.
Departure maps evaluate at the starting nodes; endpoint maps evaluate after
half a period. Only duplicate periodic planes are identified algebraically.
Connectivity is obtained from those maps, not prescribed as separate lines.
Two canonical planes use 5, 9 and 17 cells per section direction. Odd cell
counts keep the magnetic axis between nodes on every refinement.
Support and interior mass weights are `h^2/2`, where `h` is section spacing;
line lengths come from integrating `|B|` along the traced segment.
The slow matrix is the standard Cartesian Laplacian, including periodic
third-direction differences. The parallel coefficient passed to the oracle
is the excess conductivity `1/epsilon-1`.

An independent sine eigenvalue checks the slow matrix and the
perpendicular-only baseline. The latter uses excess conductivity `1e-20`
to respect the oracle's positive-coefficient contract and agrees with the
zero-excess eigenvalue solution within `1e-11` in maximum norm.
The norm is the relative diagonal-mass norm
defined below. The deliberately loose acceptance bound is `error < 0.1`
at every epsilon and refinement. The registered
`test_fci_transport_nested_accuracy_gate` enforces it and **fails all 24
cases**. It has no expected-failure property. The separate diagnostic test
checks tracing, assembly, the baseline and residuals; its passing result
does not mean the accuracy gate passed.

| cells per section direction | total unknowns | perpendicular-only error | error at epsilon=1e-6 | error at epsilon=1e-13 | parallel defect |
|---|---:|---:|---:|---:|---:|
| 5 | 72 | 0.03356 | 0.99706 | 1.00000 | 0.08789 |
| 9 | 200 | 0.01022 | 0.97841 | 1.00000 | 0.03660 |
| 17 | 648 | 0.002851 | 0.74913 | 1.00000 | 0.01069 |

The parallel defect is `sqrt(sum(W_s*(Q*psi)^2)/sum(M*psi^2))`, with `Q`
the mapped gradient, `W_s` the support weights and `M` the diagonal mass.
All quantities in this table are dimensionless. The original weak residual,
recomputed through sparse gradient and transpose actions rather than the
reduced equations, stays below `4.1e-10`. Traced flux defects stay below
`1.3e-10`; endpoint step-doubling differences stay below `2.6e-10`.
The manufactured solution is lost despite small algebraic residuals.

The cause is the bilinear interpolation in
`src/operators/fci_interpolation_map.f90:1434` and its use in the gradient
at `src/operators/fci_parallel_operator.f90:114`. Partition of unity
preserves constants but not this nonlinear flux. Positive stencil paths
connect every interior departure node to the anchored perimeter, although
the continuous interior field lines never reach it. The diagnostic checks
this reachability independently. Consequently the homogeneous free-node
parallel kernel is trivial; the oracle's numerical-rank check confirms that
fact. At fixed mesh the solution therefore tends to zero as epsilon tends
to zero. This is not a Krylov convergence failure or a theorem about other
FCI interpolation schemes.

Stop at gate 1. The next step is to establish an accurate limit space with
a different spatial/map construction or an asymptotic-preserving (AP)
discretization, then build
option B against that operator. No option-B iteration counts, option-A
implementation, or new matched-accuracy cost comparison are claimed.
The existing additive field split supplies a parallel Jacobi action, not an
exact line solve. Plane cycles and FortSparse's symmetric-positive-definite
algebraic multigrid (AMG) interface (`amg_setup_csr`, `amg_apply`,
`near_nullspace`) were inspected
but cannot repair this failed accuracy gate by preconditioning alone.

The [unfiltered CSV](figures/fci_transport_nested_accuracy.csv) retains all
24 computed rows, residuals, tracing defects and individual wall timings.
The [plot source](figures/plot_fci_transport_nested_accuracy.py) shows
solution error and original weak residual with dimensionless axis labels:

```sh
python3 doc/design/figures/plot_fci_transport_nested_accuracy.py --output /tmp/fci-transport-nested.png
cmake --build build/u3b -j 8 --target test_fci_transport_nested test_fci_transport_reference
ctest --test-dir build/u3b -j 8 -R '^test_fci_transport_' --output-on-failure
```

Use the Release/OpenBLAS configuration below with `build/u3b` as build root.
The transport scope returns seven passes and one failure, the accuracy gate.
On shared faepmac1, gfortran 16.2 and one OpenMP/BLAS thread, oracle setup
took 0.265 ms, 7.66 ms and 1.28 s on the three meshes. Individual solves
are retained in the CSV; these diagnostic timings are not gate-3 evidence.

## Algebra and boundary contract

The input is a weak, positive-sign system

\[
 [\alpha M+\chi_\perp P+\chi_\parallel Q^T W_s Q]u=f.
\]

Here `Q` is the sparse staggered support gradient, `W_s` contains positive
support weights, `M` is a positive diagonal mass, and `P` is a supplied
symmetric perpendicular weak matrix. `alpha` is a nonnegative mass shift;
`chi_parallel` is the positive parallel transport coefficient and
`chi_perp` is its nonnegative perpendicular counterpart, with
`epsilon = chi_perp/chi_parallel`. The caller supplies these physical
coefficients, a coercive slow operator on the parallel kernel, and boundary
conditions. Positive semidefiniteness of `P` is a caller precondition; the
reference checks its symmetry, not its inertia. This differs from the signed,
mass-inverted diffusion action in the existing FCI API: a caller must multiply
that action by the canonical volumes and reverse its sign to obtain this weak
matrix. `f` is a weak load, not a pointwise source rate.

`setup_fci_transport_reference` removes the prescribed degrees of freedom
before computing a singular value decomposition (SVD) of
`W_s^(1/2) Q_free M_free^(-1/2)`. Its right singular vectors give bases `Z` for
the kernel and `Y` for its complement, both orthonormal in the mass pairing.
Thus `Z Z^T M_free` is the discrete mass-weighted average projector. The
caller supplies `expected_nullity`; setup rejects disagreement with the SVD
rank at threshold `8 * machine_epsilon * max(rows, columns) * largest_sigma`.
This checks numerical rank, not magnetic topology or a certified spectral gap.
Weak but nonzero modes below that threshold cannot be distinguished from an
exact kernel. The reported `kernel_defect()` is the maximum absolute entry of
`W_s^(1/2) Q_free Z`.

Write `u_free = Z p + Y y`, and let
`B = alpha M_free + chi_perp P_free`. The two sets of equations are

\[
 \begin{pmatrix}
 Z^T B Z/s & Z^T B Y/s \\
 Y^T B Z/\chi_\parallel &
 \Sigma^2+Y^T B Y/\chi_\parallel
 \end{pmatrix}
 \begin{pmatrix}p\\y\end{pmatrix}
 =\begin{pmatrix}Z^T f_{\rm slow}/s\\f_{\rm micro}/\chi_\parallel\end{pmatrix},
 \qquad s=\max(\alpha,\chi_\perp).
\]

`Sigma` holds the retained singular values. `f_slow` subtracts the
perpendicular boundary lifting; `f_micro` additionally subtracts the parallel
boundary lifting, using the left singular vectors before scaling. The
parallel term is absent from the kernel rows. Consequently no subtraction of
large parallel stiffness entries is needed to recover the slow equations.
For forcing of order `s` and zero boundary lifting, `y` is of order
`s/chi_parallel`, giving the usual average plus small correction structure.
The macro block still needs coercivity, and the retained singular values still
control conditioning with mesh refinement. There is no uniform-in-mesh theorem
or numerical certificate for this reference.

Boundary values are prescribed at caller-supplied intersection degrees of
freedom. An anchored line has no homogeneous constant mode; closed components
can retain theirs in the same system. Nonzero boundary lifting permits the
parallel two-point solution to vary along an open line. The solver does not
trace intersections, construct cut-cell geometry, or apply a sheath law.
Rebuild setup whenever maps, topology, mass, perpendicular matrix or the set of
anchors changes. Coefficients and boundary values can change between solves.

The returned residual is the larger of the relative Euclidean residuals of
the two scaled blocks; a zero block load uses an absolute residual. It checks
the reduced solve, not the continuous equation or the unmodified input
operator after numerical kernel truncation. `status == FORTSPARSE_OK` means
that validation and the dense solve completed; callers must inspect the
residual and the kernel defect. Storage is quadratic and setup/solve are dense;
`storage_bytes()` counts retained arrays only, excluding scratch, allocation
overhead and BLAS/LAPACK workspace.

## Relation to asymptotic-preserving schemes

### Closed-line literature and reduced constraint spaces

The following locations were checked against the linked full-text preprints.
Numbering in the published versions has not been checked. None of these
papers supplies a coarsened FCI difference space with the required discrete
limit kernel and a uniform inf-sup constant.

- Deluzet and Narski, [*A two field iterated asymptotic-preserving method for
  highly anisotropic elliptic equations*, MMS 17 (2019), 434–459](https://doi.org/10.1137/17M115205X),
  [EUROfusion author preprint](https://scipub.euro-fusion.org/wp-content/uploads/eurofusion/WP14ERPR18_19317_submitted.pdf):
  Section 2.1, equations (3)–(9), assumes a smooth unit field, positive
  parallel coefficient and uniformly positive perpendicular tensor on the
  perpendicular subspace. Section 2.2 explicitly distinguishes locking from
  conditioning. Section 2.3, equations (18)–(19), uses two fields in the same
  space `V`; Lemma 5 preserves orthogonality to field-constant functions in
  the H1 seminorm when the initial auxiliary field has zero perpendicular
  gradient. Section 3.1, equation (45), discretizes both fields in the same
  `V_h`, using linear triangles or biquadratic rectangles. Its parameter
  `epsilon_0` keeps each subproblem mildly anisotropic. This is an iterated
  auxiliary-potential formulation. Its applicability to the present FCI
  kernel is **unverified**.
- Narski and Ottaviani, [*Asymptotic Preserving scheme for strongly
  anisotropic parabolic equations for arbitrary anisotropy direction*, CPC
  185 (2014), 3189–3203](https://arxiv.org/pdf/1303.5219):
  Section 3.1, equation (4), removes the inflow gauge and uses an auxiliary
  potential in `L = {q in L2 : grad_parallel q in L2}`, with a positive mass
  stabilization. Equation (5) sets its coefficient to `h^(k+1)` for degree
  `k` elements; Section 4.1 uses biquadratic elements on a Cartesian mesh.
  This permits closed lines without prescribing a point on every line.
  The stabilization changes the second equation; it is not a reduced
  first-derivative multiplier space. The paper postpones theoretical
  justification of this construction.
- Lozinski, Narski and Negulescu, [*Numerical analysis of an
  asymptotic-preserving scheme for anisotropic elliptic equations*](https://arxiv.org/pdf/1507.00879):
  Section 3, equations (34)–(38), uses the same parallel-Sobolev space and
  analyzes its L2-orthogonal complement to field-constant functions.
  Section 3.2, equations (50)–(52), uses `V_h` for both variables, then
  removes the discrete field-constant component for analysis. Its inf-sup
  result uses a mesh-dependent dual norm, not an ordinary L2 multiplier
  norm. Theorem 22, equation (53), assumes regular-mesh approximation and
  the stated H^(k+1) and H1 solution regularity. Remark 23 only expects
  uniform higher regularity in general and proves it in simple aligned
  geometry. These hypotheses are not established for our sheared FCI map,
  islands or X-points.
- The earlier [Degond–Lozinski–Narski–Negulescu micro-macro preprint](https://arxiv.org/pdf/1102.0904),
  Section 2.2, equation (13), fixes the auxiliary potential on inflow;
  Hypothesis B and Theorem 3 require a smooth tube of boundary-crossing
  field lines. The [Degond–Deluzet–Lozinski–Narski–Negulescu duality
  preprint](https://arxiv.org/pdf/1008.3405), Section 2.3, Hypothesis B,
  equations (2.16)–(2.18), assumes an H1-continuous field-line averaging
  projector and a parallel Poincare–Wirtinger inequality. These assumptions
  do not give a closed-line FCI multiplier construction.
- A checked successor is Yang, Deluzet and Narski,
  [*Preserving the accuracy of numerical methods discretizing anisotropic
  elliptic problems*](https://arxiv.org/pdf/1911.11482), Section 3.2,
  equations (32)–(33). It compares the inflow micro-macro construction with
  the two-field iteration and rescaled parallel dynamics. It does not
  specify a reduced FCI difference test space. The later published version,
  [JCP 521 (2025), 113568](https://doi.org/10.1016/j.jcp.2024.113568),
  was located; its full-text equation numbering and any changes from this
  preprint are **unverified**.

The auxiliary potentials in these papers must not be identified with the
first-derivative flux multiplier in the proposed reduced mixed method.
In particular, epsilon-independent conditioning of an auxiliary system does
not by itself establish that its FCI limit space approximates flux functions.
For the unchanged FCI bilinear form, there is also a direct algebraic
obstruction: at an exact two-field fixed point and epsilon zero, the second
equation enforces the same parallel kernel. Its free-node kernel is trivial
on this fixture. This is an inference about transplanting the construction
to this operator, not a claim that the papers analyze this FCI map.

The kernel/complement split follows the limit-space and fluctuation-space
construction in Degond, Deluzet, Lozinski, Narski and Negulescu,
[*Duality-based Asymptotic-Preserving method for highly anisotropic diffusion
equations*, Communications in Mathematical Sciences 10 (2012), 1–31](https://doi.org/10.4310/CMS.2012.v10.n1.a2),
Section 2.3. The small-correction interpretation follows Degond, Lozinski,
Narski and Negulescu,
[*An Asymptotic-Preserving method for highly anisotropic elliptic equations
based on a micro-macro decomposition*, Journal of Computational Physics 231
(2012), 2724–2740](https://doi.org/10.1016/j.jcp.2011.11.040),
Section 2.2 and Remark 2. This is an explicit finite-dimensional decomposition,
not an implementation of either paper's finite element multiplier system.
In particular, the latter paper's inflow-based gauge and its tube-domain
hypothesis do not establish a closed-line or X-point result here.

The algebraic measurements below justify retaining this decomposition as an
oracle. Option B is the selected solver design, but its spatial accuracy
preflight above fails before implementation. FortSparse's available smoothed
aggregation interface (`amg_setup_csr`, `amg_apply`,
`amg_params_t%near_nullspace`) and FortFEM's plane cycles were inspected; neither
is composed into this reference. Comparing those scalable alternatives,
including a parallel complement solve and perpendicular multigrid, remains the
next part. A dense direct solve has no Krylov iteration count to report.

## Verification and reproduction

`test_fci_transport_reference` constructs gradients through the existing FCI
map assembler, identifies periodic end planes, and checks discrete manufactured
solutions with analytic Fourier eigenvalues. One fixture has three periodic
lines. Another has two periodic lines coupled perpendicularly to an open line
with plate values 2 and 5. A separate nontrivial interpolation-map fixture has
one constant mode and nonuniform support/mass weights. These are algebraic
fixtures; their refinement is not a spatial-convergence study of a non-aligned
magnetic field. They do not reproduce the Sovinec or Guenter benchmarks.

The accuracy norm is `sqrt(sum(M*(u-u_exact)^2)/sum(M*u_exact^2))`.
The periodic and mixed ladders cover every decade from `1e-6` to `1e-13` on
three refinements. Periodic checks include zero and nonzero mass shifts and
both fixed-parallel and fixed-perpendicular normalizations. The enforced bound
is `2e-11`; the observed maxima are recorded below. Invalid rank, duplicate
anchors, negative support weights, nonfinite coefficients and an unregularized
kernel are rejected. A fully anchored pure-parallel solve is checked separately.

On faepmac1, GNU Fortran 16.2.0, Release `-O3`, OpenBLAS with one BLAS thread:

```sh
cmake -S . -B build/u3a -G Ninja -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_Fortran_COMPILER=/opt/homebrew/bin/gfortran \
  -DFORTFEM_BUILD_TESTING=ON -DBLA_VENDOR=OpenBLAS \
  -DCMAKE_PREFIX_PATH=/opt/homebrew/opt/openblas
cmake --build build/u3a -j 8 --target test_fci_transport_reference
ctest --test-dir build/u3a -j 8 -R '^test_fci_transport_' --output-on-failure
ctest --test-dir build/u3a -j 8 -V -R '^test_fci_transport_'
```

The native fpm target also passes (359 assertions) with the full LAPACK provider:

```sh
FO_JOBS=8 OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=1 \
  fpm test test_fci_transport_reference --profile release \
  --link-flag /opt/homebrew/opt/openblas/lib/libopenblas.dylib
```

The default Accelerate provider on this host lacks `dsytrf_rk`, called by the
pinned native sparse backend; use the installed full OpenBLAS provider for the
comparison. The cost tests print setup time, mean solve time over ten loads,
relative mass error and retained bytes. Reference setup includes the SVD;
each reference solve factors its dense scaled matrix anew. The sparse direct
path factors the assembled original matrix once and reuses that factor.
These are different work contracts, and their timings are not a comparison at
matched accuracy. Eight OpenMP threads do not accelerate these small dense
problems; the dominant LAPACK work remains on one BLAS thread.

The broader registered `fci` scope currently has two baseline failures on this
host: `test_fci_quartic_curved_polygon_area.f90` rejects its derivative
finite-difference discrepancy, and `test_fci_hendecic_curved_polygon_area.f90`
rejects its Gauss–Green area discrepancy. Both reproduce with original main
sources in a separate baseline consumer: respectively `7.3615e-8 > 5e-8` and
`8.4572e-12 > 8e-12`. Their checks and generated kernels remain unchanged.
Therefore the full existing-FCI verifier is not green, although the new
transport checks pass. These geometry failures do not establish a transport
regression.

## Measured accuracy and cost

The following relative mass errors are dimensionless. Each column uses the
number of unknowns `n` shown in its heading; the periodic rows shown here use
`alpha = 0`. The mass-shifted and coefficient-rescaled cases also pass the
asserted bound. No error growth with decreasing epsilon is observed in these
fixtures.

The retained plot source uses all 48 printed rows from the zero-mass-shift ladders, without
filtering or smoothing. The [CSV data](figures/fci_transport_reference_accuracy.csv)
and [plot source](figures/plot_fci_transport_reference_accuracy.py) are adjacent;
run the source with Python and Matplotlib and `--output /tmp/fci-reference.png`
to reproduce the figure. Rendered figures are not checked into the repository.

| epsilon | periodic n=12 | n=24 | n=48 | mixed n=13 | n=25 | n=49 |
|---|---:|---:|---:|---:|---:|---:|
| 1e-6 | 1.029e-13 | 1.097e-14 | 5.782e-15 | 4.109e-15 | 2.410e-15 | 4.224e-15 |
| 1e-7 | 1.029e-13 | 1.132e-14 | 6.264e-15 | 3.707e-15 | 2.588e-15 | 3.924e-15 |
| 1e-8 | 1.027e-13 | 1.048e-14 | 5.419e-15 | 4.046e-15 | 2.724e-15 | 4.232e-15 |
| 1e-9 | 1.023e-13 | 1.105e-14 | 5.603e-15 | 3.963e-15 | 2.556e-15 | 3.875e-15 |
| 1e-10 | 1.027e-13 | 1.223e-14 | 5.769e-15 | 3.841e-15 | 2.465e-15 | 4.198e-15 |
| 1e-11 | 1.027e-13 | 1.051e-14 | 5.813e-15 | 4.130e-15 | 1.736e-15 | 3.595e-15 |
| 1e-12 | 1.034e-13 | 1.143e-14 | 5.890e-15 | 3.701e-15 | 2.565e-15 | 4.093e-15 |
| 1e-13 | 1.028e-13 | 1.015e-14 | 7.431e-15 | 4.197e-15 | 2.374e-15 | 4.021e-15 |

The largest printed block-scaled residuals are `4.524e-16` (periodic) and
`8.967e-15` (mixed). The mass error stays below `1.035e-13` in the printed
periodic ladder and `4.232e-15` in the mixed ladder.

Timings below are elapsed milliseconds on the shared faepmac1 host. Each
entry gives one OpenMP thread / eight OpenMP threads; BLAS uses one thread
in both cases. These are individual measurements, not confidence intervals.

| method | n | epsilon | setup (ms), 1 / 8 | mean solve (ms), 1 / 8 | relative mass error |
|---|---:|---|---:|---:|---:|
| reference | 48 | 1e-06 | 0.316 / 0.443 | 0.0219 / 0.0857 | 5.782e-15 |
| direct | 48 | 1e-06 | 0.531 / 0.219 | 0.0017 / 0.0015 | 9.356e-09 |
| reference | 48 | 1e-13 | 0.246 / 0.317 | 0.0179 / 0.0795 | 7.431e-15 |
| direct | 48 | 1e-13 | 0.212 / 0.156 | 0.0016 / 0.0017 | 3.251e-01 |
| reference | 96 | 1e-06 | 2.002 / 2.536 | 0.0597 / 0.1183 | 4.144e-14 |
| direct | 96 | 1e-06 | 0.356 / 0.352 | 0.0021 / 0.0021 | 1.950e-08 |
| reference | 96 | 1e-13 | 1.938 / 2.458 | 0.0579 / 0.1124 | 4.505e-14 |
| direct | 96 | 1e-13 | 0.343 / 0.342 | 0.0024 / 0.0024 | 2.411e-01 |
| reference | 192 | 1e-06 | 17.32 / 21.75 | 0.251 / 0.7322 | 2.035e-14 |
| direct | 192 | 1e-06 | 1.986 / 1.957 | 0.0037 / 0.0039 | 4.082e-08 |
| reference | 192 | 1e-13 | 16.48 / 21.06 | 0.2473 / 0.7528 | 2.153e-14 |
| direct | 192 | 1e-13 | 1.951 / 1.979 | 0.0044 / 0.0043 | 1.453e+00 |

The reference retains 55,848 / 222,312 / 887,016 bytes at `n = 48 / 96 / 192`.
The native sparse backend does not expose factor-storage or workspace counts
for these cases (`-1` in the printed output). As a separate process-level
measurement, `/usr/bin/time -l` reports peak resident set sizes of 8.73 /
11.36 MiB for the reference and 6.62 / 6.66 MiB for the direct path, with one /
eight OpenMP threads. One MiB is 1,048,576 bytes. Each process runs both epsilon
endpoints at all three sizes, including fixture construction, setup and ten
solves; these peaks are not isolated per-solve allocation measurements. Reproduce
with the following command, replacing `reference` with `direct` and `1` with `8`:

```sh
/usr/bin/time -l env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  build/u3a/test_fci_transport_reference reference
```

At `epsilon = 1e-13`, the assembled sparse solve has relative errors from
0.24 to 1.45 on these fixtures, whereas the reference remains below `4.506e-14`
in the cost cases. This supports the kernel separation as an accuracy oracle.
It does not establish a practical crossover point or select the production
algorithm. The direct timing cannot be interpreted as time to the reference
accuracy.
