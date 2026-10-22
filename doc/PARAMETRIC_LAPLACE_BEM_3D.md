# General parametric-surface Laplace BEM

[V] `fortfem_parametric_laplace_bem_3d` provides a scalar Galerkin path for a
fixed, connected, closed, oriented triangular surface. It does not assume
axisymmetry. Potential traces are continuous P1 and normal fluxes are P0.
Geometry is evaluated through a caller-owned map and may be curved to any
order, including an exact analytic map. High-order **trace spaces** are not
implemented by this API. `fortfem_parametric_bem_tensor` provides continuous
Q_p potentials and discontinuous Q_pf fluxes on quadrilateral panels of the
same surface contract, with the same singular treatment and operator
container. Its FEM coupling is described in
[the tensor FEM/BEM design](design/tensor_fem_bem.md).

## Surface and operator API

[V] Extend `bem_surface_t` and implement the pure procedure
`evaluate(self,panel,uv,x,a,b,status)`. The reference triangle is
`uv(1)>=0, uv(2)>=0, sum(uv)<=1`; `a` and `b` are the position derivatives
with respect to its two coordinates. The map must be regular, agree between
neighboring panels, and orient `a cross b` out of the bounded region. Pure
maps must be safe to call concurrently. `faces(3,npanels)` supplies the
matching P1 vertex numbering and oriented topological connectivity.

```fortran
call assemble_parametric_laplace_bem(surface, faces, degree, op, status)
call build_exterior_dtn(op, exterior, status)
call exterior%apply(dirichlet_trace, outward_flux, exterior_load, status)
call exterior%matrix(steklov, status)
```

[V] Assemble and build once per fixed geometry, connectivity and quadrature
choice. Each later apply is a matrix action, without assembly or factorization.
The container has `v`, `k`, `adjoint`, `w`, `mass`, `area` and `genus`.
`v` acts on panel flux coefficients; `k` and `mass` map vertex traces to
panel tests; `w` maps vertex traces to vertex tests. `adjoint=transpose(k)`.
All operators and accepted DtN actions use fp64.

[V] With `n` pointing into the exterior and `G=1/(4*pi*|x-y|)`, the physical
outward derivative satisfies `V q=(K-M/2)u`. The symmetric exterior energy map
is `S=W+(M/2-K)^T V^{-1}(M/2-K)`. `apply` returns this positive energy load
`S u` as well as the separate physical flux `q`. These have different discrete
spaces and must not be interchanged. `matrix` returns a copy suitable for the
FEM restriction `R^T S R`. A caller must add its known exterior drive using
the same normal and trace conventions.

[V] Status 1 rejects inputs or topology; status 2 rejects geometry or
nonfinite assembly values. The DtN builder additionally propagates the dense
solver status. An unsuccessful output must not be applied. Changing a public
operator array requires rebuilding the cached DtN object.

## Quadrature and verification limits

[V] Coincident-panel integration uses a source Duffy fan centered at each
target point. Close-panel integration uses a source fan about a projected
reference point and a squared radial coordinate to cluster samples. The
double layer integrates `d_n G (u(y)-u(x))` and restores the exact
closed-surface term `-u(x)/2`. The hypersingular form uses surface curls,
reducing it to a weakly singular integral. Far pairs use product triangle quadrature. Assembly is OpenMP;
shared hypersingular entries use atomic additions.

[V] `degree` controls quadrature, not approximation order. This implementation
has no tolerance-controlled near-pair error bound or Sauter-Schwab edge/vertex
pair rule. A user must demonstrate quadrature convergence and mesh convergence
for the actual shell, particularly for close but nonadjacent panels.

[V] FortSym generates the Green, normal, Jacobian, curl-density and integrand
leaves, including JVP leaves, from `tools/codegen/app/gen_parametric_bem_products.f90`.
Generated output is checked in. The assembled operator currently has no
geometry-JVP/VJP API; fixed near/far selection and the surface-map derivatives
must be included before claiming differentiable assembly.

[V] `parametric_bem_residual(op,u,q,panel_indicator,status)` reports the
area-scaled discrete first-Calderon equation defect. With the assembly operator
used to solve for `q`, it measures algebraic error only. With an independently
more accurate quadrature operator on the same mesh, it can reveal quadrature
error. It is neither a continuous error bound nor a discretization certificate.
Use the refined-mesh estimator below to measure a change of discrete space.
The independent sphere and torus
oracles are in `test/test_parametric*`.

[V] Manufactured flux tests compare the P0 coefficients with the exact
area-averaged normal derivative on each panel. Unresolved variation inside a
panel is excluded from that coefficient error. A continuous field-error norm
and estimator are still needed for a total discretization-error claim.

## Refined-mesh flux estimator

[V] Solve on a coarse mesh and a nested refinement of the **same physical
surface**, with independently checked quadrature convergence on each mesh.
Supply the resulting P0 fluxes, the fine operator, and `parent(j)`, the coarse
panel containing fine panel `j`:

```fortran
call parametric_bem_refinement_estimate(fine_op, coarse_flux, fine_flux, &
    parent, panel_indicator, energy_estimate, status)
```

[V] The caller owns refinement and the physical parent map. Every coarse
panel must have children; these children must partition that panel without
overlap. The routine validates array shapes, finite values, positive areas,
parent indices and coverage, operator symmetry and nonnegative difference
energy. It cannot verify geometric containment from the operator arrays.

[V] With constant injection `I`, `d = q_f - I q_c`, the returned global
estimate is `sqrt(d^T V_f d)`. It approximates the single-layer energy norm
of the flux change, including changes caused by refined Dirichlet data.
For marking, each coarse panel receives
`eta_T = sqrt(h_T sum_{j:parent(j)=T} area_j d_j^2)`, where
`h_T = sqrt(sum_{j:parent(j)=T} area_j)`. These nonnegative local quantities
use a mesh-weighted L2 norm; their squared sum is **not** the global
single-layer energy. Shape regularity is needed to relate the local mesh
weight to a geometric panel diameter. The product `d,d^2` and its JVP are
generated by FortSym; the square identity is checked symbolically.

[V] This is an h/h2 a posteriori **estimate**, conditional on saturation and
accurate quadrature, rather than an unconditional upper bound. To see the
condition, let `E_c = ||q-I q_c||_V`, `E_f = ||q-q_f||_V` and assume
`E_f <= rho E_c` for a specified `rho < 1` in the continuous single-layer
norm. Triangle inequalities give
`eta/(1+rho) <= E_c <= eta/(1-rho)`, where `eta=||d||_V`.
No Galerkin orthogonality is assumed. The routine does not establish `rho`
or bound the replacement of continuous `V` by quadrature `V_f`.
It also does not estimate a separate potential-trace norm, geometry error,
or the error of nonexact harmonic sectors. Near-singular quadrature still
needs independent refinement; an algebraic residual cannot replace this check.

[V] `test_parametric_bem_convergence_slow` audits the exterior sphere dipole
`u=z/r^3`, whose exact outward flux on the unit sphere is `q=-2z`.
Its independent continuous L2 integral includes unresolved within-panel
variation. It requires decreasing true flux error, global energy estimates
and local marking estimates on 8, 32 and 128 panels, while every solved
equation defect remains small. It reports marking effectivity against the
coarse mesh-weighted continuous error; those observations do not prove
saturation for another surface or data set.

## Topology and remaining integration

[V] Connectivity checks reject open edges, nonmanifold edges, inconsistent
edge orientation and disconnected surfaces, and report the orientable genus.
The scalar hypersingular form annihilates the constant trace; the decaying
exterior Dirichlet map itself does not require removing that constant.

[V] A single-valued scalar potential has zero periods. A genus-one exterior
can also carry a nonexact harmonic one-form. Its circulation, and the matching
interior harmonic flux data, are separate data, not extra null vectors of the
scalar single-layer matrix. This API does not construct or couple those fields.
FortFEM's existing `cell_complex` and `harmonic_period_normalization` modules
provide topology and period-normalization primitives for a subsequent adapter.

[V] This API is not a complete general magnetic exterior closure. High-order
trace spaces, an unconditional error certificate, harmonic-sector coupling,
assembled geometry derivatives, compressed assembly/solves and controlled mixed
precision remain open. The existing per-mode axisymmetric API remains the
consumer's established route until these gaps and its accuracy gates close.

[V] Operator and weak-form conventions follow the
[Bempp boundary-operator definitions](https://bempp.com/handbook/api/boundary_operators.html)
and the existing FortFEM `toroidal_mode_boundary_operators_2d` sign conventions.

## Curved P2/P1 traces and bounded pair refinement

[V] `fortfem_parametric_bem_high_order` adds continuous quadratic potentials
and discontinuous linear fluxes on the same `bem_surface_t` geometry. The
geometry connectivity remains `faces(3,npanels)`. The returned
`trace_dofs(6,npanels)` numbers vertices followed by shared edge midpoints
in local edge order (1,2), (2,3), (3,1); `flux_dofs(3,npanels)` numbers
three independent reference-P1 coefficients per panel. Geometry order and
trace order are independent. Maps must agree on shared edges.

```fortran
type(bem_pair_options_t) :: options
type(bem_pair_report_t) :: report
options%relative_tolerance = 1.e-3_dp
options%absolute_tolerance = 1.e-9_dp
call assemble_parametric_bem_p2(surface, faces, options, op, &
    trace_dofs, flux_dofs, report, status)
! Use build_exterior_dtn and exterior%apply as above, with the new DOF counts.
```

[V] Pair refinement increases both target and source quadrature degrees.
Repeated degree labels that select the same Gauss rule are skipped; only
distinct rules count as refinements. An interior closest-point normal peak
must also be sampled within its local distance scale before two rules may
claim convergence; unresolved peaks return status 3 at the caps. Coincident
panels use the product-domain Duffy rule described below. Edge-touching,
vertex-touching and separated pairs retain the projected squared-radial source
fan. For each
local V, K, W and constant-subtraction block, acceptance requires the maximum
entry change to be at most `absolute_tolerance + relative_tolerance *
maxval(abs(current_block))`. The tolerances apply to local blocks, not to
the final solution. Degree and kernel-evaluation caps are per pair;
status 3 reports failure to meet this criterion. Invalid options/topology
return 1; invalid geometry/nonfinite data return 2. An unsuccessful operator
must be discarded. `maximum_degree` is limited to 128.

[V] This is an observed quadrature error estimate, not a rigorous error
bound. Repeated rules can share an unresolved bias; independently refine
quadrature and the mesh for a new geometry. `report` retains work counts,
the largest accepted error ratio and raw V/W symmetry defects measured
before averaging. The raw local-pair API
`integrate_parametric_bem_p2_pair` returns V, K, W and the target-side
constant subtraction without averaging or adding a closed-surface jump.

[V] Assembly batches the target samples before parallel pair integration;
there are no allocations in pair loops. W scatter uses OpenMP atomics.
`op%area` retains geometric panel areas, while V has three rows per panel.
The P0 parent-injection estimator must not be used for these flux DOFs.
Generated P1/P2 basis jets and generated product/JVP/VJP leaves supply the
pointwise mathematics. Assembled geometry derivatives remain unavailable.
The independent pair test uses analytic signed-edge source potentials and
refined target quadrature; the slow sphere audit measures continuous flux
error and exterior energy for degrees 1 and 2, observed order under nominal dyadic mesh refinement, and
1/2/4/8-thread assembly scaling.

Decaying exterior harmonic representatives, period compatibility and their
energy/cross blocks are described in [EXTERIOR_HARMONIC_FIELDS.md](EXTERIOR_HARMONIC_FIELDS.md).

## Rotational block-circulant P2/P1 assembly

[V] `assemble_parametric_bem_p2_circulant` in
`fortfem_parametric_bem_high_order` accepts `nphi` consecutive, identical panel
sectors related by rotation through `2*pi/nphi` about the z axis. Corresponding
panels must use identical local vertex order. The trace rotation is inferred
from connectivity, including shared edges across sector boundaries; global
vertex and edge numbers need not be contiguous within a sector. The cyclic
action on trace DOFs must be free: fixed-axis vertices are unsupported.

```fortran
use fortfem_bem_circulant, only: bem_circulant_t
! Other types and options are the same as for all-pairs P2 assembly.
type(bem_circulant_t) :: circular
call assemble_parametric_bem_p2_circulant(surface, faces, nphi, options, &
    circular, trace_dofs, flux_dofs, report, status)
call circular%expand(op, status) ! Existing bem_operators_t, original DOF order.
```

[V] The assembler evaluates target-sector-zero pairs against every source
sector, using the same adaptive pair rules and tolerance/cap checks as the
all-pairs routine. `accepted_pairs`, `failed_pairs` and `kernel_evaluations`
count actual work. V/W symmetry averaging and the constant-subtracted K jump
are retained. Storage and pair counts are reduced by `nphi`; dense expansion
is optional. Geometry and tangent rotation are checked at each panel's three
vertices, three edge midpoints and centroid. Invalid sector counts,
connectivity/orbits or sampled symmetry return status 1; invalid/nonfinite
geometry returns 2; unresolved pairs return 3. [H] Exact rotational symmetry
between these samples remains a caller-owned geometry contract.

[V] `circular%v`, `%k`, `%w` and `%mass` are `circulant_matrix_t` operators.
Their `blocks(:,:,d+1)` hold the first block row B(d), with
A(r,s)=B(modulo(s-r,nphi)). `row_dofs` and `column_dofs` map sector coordinates
to original global DOFs. `apply(x,y,status)` uses these maps without expanding;
`expand(matrix,status)` produces a dense matrix. `harmonics(modes,status)`
returns the complex independent matrices in sector-frequency order 0 through
`nphi-1`, using FortNum's FFT:

    A_hat(k) = sum_d B(d) exp(+2*pi*i*k*d/nphi).

Sector vectors use the negative-exponent forward DFT; the inverse carries
`1/nphi`. `solve(rhs,x,status)` solves square operators by transforming the
right-hand side, solving each complex system through FortNum, and transforming
back. A nonzero solver status must be checked. The current solve factors anew
on each call; consumers with repeated loads can retain the exposed harmonic
matrices and factor them once. No assembled geometry derivatives are supplied.

[V] `test_bem_circulant_modes` checks nonsymmetric solves for sector counts
1–13 against an independent dense solve. `test_parametric_bem_circulant`
compares all four BEM forms, areas, application, and interior/exterior harmonic
Dirichlet solves against all-pairs assembly on the existing curved torus.
`test_parametric_bem_circulant_scaling_slow` repeats at 12 sectors and reports
all-pairs/sector assembly costs plus 1/2/4/8-thread timings, checking the same
matrix tolerance at every thread count. Timings are measurements, not portable
speed thresholds.

## Coincident-panel product rule

[V] For reference triangle T, the difference y-x lies in a hexagon. Its six
sectors have consecutive corner differences d(s), with 0<s<1. On each sector,
write x=r*t(s)+(1-r)*z and y=r*u(s)+(1-r)*z, where z is in T and t,u are the
piecewise affine negative/positive barycentric parts of d. The four-dimensional
Jacobian is r*(1-r)^2. Thus the weak 1/r singularity is removed in both panel
coordinates, including for regular curved maps, instead of leaving a target
integral with edge logarithms. The surface Jacobians and basis/curl densities
are still evaluated on the actual geometry.

[V] A tangent-metric sinh change of s resolves the narrow angular peak on skinny
panels. Both maps and their Jacobian checks are generated by FortSym. Gauss
refinement, local-block acceptance and the degree/work caps retain their
existing meanings. One coincident rule costs 6*n^4 kernel evaluations, with
n=(degree+3)/2, and it allocates no storage in the pair loop. This change does
not resolve arbitrarily close disjoint panels beyond the existing fan's caps.

[V] The independent planar oracle is the exact self energy for a right triangle
with perpendicular sides a,b. With G=1/(4*pi*|x-y|), its integral is

    (a*b)^2/(12*pi) * ((asinh(a/b)+asinh(b/a))/sqrt(a^2+b^2)
                       + asinh(a/b)/a + asinh(b/a)/b).

The pair test checks this before symmetry averaging, for aspect ratios 1 and
40. The curved-sphere test additionally checks the exact identity D_n G=-G/2
on the unit sphere and raw V/W symmetry at a tighter pair tolerance. Global
sphere-harmonic mesh convergence remains the independent discretization oracle.
