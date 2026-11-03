---
title: Anisotropic heat-transport pollution benchmarks
---

# Anisotropic heat-transport pollution benchmarks

On the two-dimensional benchmarks, high-order H1 elements (continuous elements
with square-integrable first derivatives), the Günter symmetric stencil
and flux-coordinate-independent (FCI) reconstruction control
anisotropic pollution through parallel-to-perpendicular diffusivity ratios of
`1e10`. The two stencils reach the manufactured temperature-error target at
961 free unknowns at all three measured anisotropies; degree-four H1 needs
225, 225 and 961. This is practical evidence over the sampled range, not a
proof of anisotropy-independent error at arbitrary ratios or fixed mesh.

The three-dimensional bilinear FCI operator maps field-line endpoints between mesh
planes by interpolation (plane-to-plane maps). It locks on the tested nested
magnetic flux surfaces (surfaces tangent to the magnetic field, enclosing one
another), which are not aligned with the mesh. On each of three fixed meshes, its
relative solution error tends to one as `epsilon = chi_perp/chi_par`, the
perpendicular-to-parallel diffusivity ratio, tends to zero. Here `chi_perp` and
`chi_par` are the perpendicular and parallel diffusivities. This occurs despite
small algebraic residuals. A preconditioner changes how the discrete equations are
solved, not their solution, so it cannot remove this discretization error. No
production transport solver is qualified by this work. The two-dimensional FCI
fixture traces the exact analytic field with numerical Runge–Kutta steps; it has
no distinct plane coordinate and must not be identified with the three-dimensional
operator.

The problem is numerical perpendicular transport induced by a discrete
parallel operator, which should remain well below physical perpendicular
transport at physical anisotropy. Epsilon-uniform accuracy for general loads
on a fixed mesh is not an accepted target: resonant layers and arbitrarily
small parallel singular values require additional regularity assumptions.
The benchmarks measure signed pollution, independent temperature error and
cost directly. No fitted power law is an accuracy gate or extrapolation rule.

Here `T` is temperature, `Q` the heat source, `b` the unit field direction,
`chi_par` and `chi_perp` the parallel and perpendicular diffusivities, `p` the
element degree, and `h` the Cartesian element width. Coordinates, temperatures and
diffusivities below are normalized as in each source. Write `R = chi_par/chi_perp = 1/epsilon` for the anisotropy, where `epsilon = chi_perp/chi_par` is the
perpendicular-to-parallel diffusivity ratio. The coordinates `x`, `y`, `z` are
Cartesian; `pi` is the circle constant. The L2 norm is the square root of the
integral of the squared error. For the Sovinec test, `chi_perp = 1`, and the
reported paper diagnostic is `|1/T(0,0)-1|`. It includes perpendicular
discretization error. The isolated, signed pollution measure used in the fit is
instead

`E = chi_perp,num/chi_perp = 1/T_R(0,0) - 1/T_1(0,0)`.

Here `T_R` is the discrete solution at anisotropy `R`, and `T_1` its isotropic
baseline on the same mesh. Thus `E` is a signed centre-point diagnostic, evaluated
at `(0,0)`; it does not bound the total solution error. The symbol `chi_perp,num`
denotes the effective numerical perpendicular diffusivity inferred from this
diagnostic. No absolute value is applied to this difference. A fit to `E = C R^b h^a` estimates the mesh exponent `a`, anisotropy exponent `b` (a scalar fit
parameter, distinct from the field direction), and constant `C`; all quantities in
this law are dimensionless in the stated normalization.

## Checked literature

**Verified, Sovinec et al., JCP 195 (2004), 355–386:**
[DOI](https://doi.org/10.1016/j.jcp.2003.10.004),
[author manuscript](https://w3.pppl.gov/cemm/Project/NIMROD_JCP_rev1_wfigs.pdf),
Section 4.2, manuscript pp. 35–36, Fig. 9 on PDF p. 73.
The square is `[-1/2,1/2]^2`, with zero Dirichlet temperature,
`Q = 2*pi^2*cos(pi*x)*cos(pi*y)` and exact
`T = cos(pi*x)*cos(pi*y)/chi_perp`. An out-of-plane electric drive with
the source's spatial dependence induces the magnetic field; large mass
density suppresses motion. The continuous field is tangent to temperature
contours. Tensor elements of degrees 1–5 use `(p+1)^2` Gauss points.
The plotted anisotropies are `1e3`, `1e6`, `1e9`. Fig. 9's bilinear point
at `h = 1/8`, `chi_par = 1e3` is approximately 6. The deliberately broad
reading interval `[5,7]` was fixed before execution. This is a figure
reading, **not a published error bar**. There is no numerical table or
stated graphical accuracy. The journal version's figure coordinates have
not been independently checked against this manuscript.

**Verified, Günter et al., JCP 209 (2005), 354–370:**
[published PDF](https://www.pas.rochester.edu/~shuleli/0629/gunter2007.pdf),
[DOI](https://doi.org/10.1016/j.jcp.2005.03.021).
The PDF filename is misleading; its title page is the 2005 article.
Equations (6a,b), (7a,b), (8), pp. 356–357, define the symmetric
cell-centred parallel gradient and its negative adjoint divergence.
The circular test on p. 357/Fig. 2 uses `[-1/2,1/2]^2`,
`r = sqrt((x-x0)^2+(y-y0)^2)`, the distance from the field centre
`(x0,y0)`, `T = 1-r^3`,
`chi_perp = 1`, source `9*r`, and analytic Dirichlet boundary data.
The centre is `(0,0)` for mostly closed lines and `(-1/2,-1/2)` for
open lines. The normalized central error divides by
`|T(0,0)-T(-1/2,0)|`. The caption reports indistinguishable symmetric
errors for `1 <= chi_par < 1e9`; Fig. 6, p. 361, compares the
cosine-cosine Sovinec problem at `1e3`, `1e6`, `1e9`.
Neither figure supplies a table or an uncertainty. The sine-cosine test
on p. 359 instead uses `sin(2*pi*x)*cos(2*pi*y)` and analytic or periodic
boundary data; it must not be substituted for the cosine-cosine test.

**Verified, Günter et al., JCP 226 (2007), 2306–2316:**
[published PDF](https://www.pas.rochester.edu/~shuleli/0706/gunter2007.pdf),
[DOI](https://doi.org/10.1016/j.jcp.2007.07.016).
Section 2.1, equations (3)–(4), pp. 2307–2308, represents the scalar
parallel flux in discontinuous elements one degree lower than temperature.
The bilinear and biquadratic temperature cases therefore use constant and
bilinear flux spaces. Fig. 1 compares these with conforming elements at
anisotropies `1e3`, `1e6`, `1e9`; the hybrid curves remain almost unchanged.
Section 2.2, equations (5)–(7), pp. 2309–2310, gives the fourth-order
symmetric difference extension. The p. 2308 text prints
`sin(pi*x)*cos(pi*y)` while also specifying `|1/T(0,0)-1|` at the magnetic
axis; these are incompatible. Its reference to the NIMROD test and Fig. 1
are not authority to silently change that printed function. No tabulated
error values or numerical accuracy for figure readings are stated.

**Verified, Stegmeir et al., CPC 198 (2016), 139–153:**
[DOI](https://doi.org/10.1016/j.cpc.2015.09.016),
[checked preprint, arXiv:1505.02040v1](https://arxiv.org/pdf/1505.02040v1).
Section 2.5, equations (29)–(36), uses a uniform inclined field,
perpendicular step `h`, plane spacing `Delta_z`, fractional shift
`f = Delta_z*tan(theta)/h`, where `theta` is the field inclination to the
plane normal, and field-line step length `Delta_s`.
Its linear-interpolation naive diffusivity is
`f*(1-f)*h^2/Delta_s^2` for unit parallel diffusivity. The support scheme
instead has fourth-order modal hyperdiffusion. Section 3.1 uses the
circular shell `0.1 <= rho <= 0.2`, where `rho` is cylindrical radius, and radial sine modes, equation (43).
For zonal modes the measured decay is entirely numerical. Section 3.1.2,
Fig. 11 and equations (47)–(48), reports `h^2`/`h^4` for naive linear/cubic
interpolation, `h^4` for support-linear and an `A*h^8+B*h^4` contribution, with fitted coefficients `A` and `B`,
for support-cubic, each divided by `Delta_z^2`. The cubic lower-order term
is attributed to boundary interpolation. These are modal decay scalings,
not universal scalar perpendicular diffusivities. Published-version
equation numbering is unverified.

**Verified, Held et al., CPC 199 (2016), 29–39:**
[checked preprint, arXiv:1505.05383](https://arxiv.org/pdf/1505.05383).
Section V, equations (46)–(47), uses cylindrical coordinates with
`B = 20*grad(phi) + grad(psi_p) cross grad(phi)`,
`psi_p = cos(pi*(R-10)/2)*cos(pi*Z/2)`,
`R in [9,11]`, `Z in [-1,1]`, `phi in [0,2*pi]`.
Here `R` is cylindrical radius (rather than the anisotropy), `Z` is height,
`phi` is azimuth, `B` is the magnetic field and `psi_p` is poloidal flux.
LDG denotes local discontinuous Galerkin; dG denotes discontinuous Galerkin.
Section V.B uses Neumann walls, `chi_par = 100`, `chi_perp = 0`, and
time step `0.001` in normalized time `t`. Section V.B.2 initializes the stationary profile
`f0 = 0.1*psi_p^2` and measures relative L2 profile loss.
Table VII at `t = 5` fixes `N_phi = 10`, `P_phi = 1`, `P_RZ = 3`;
`N_phi` counts azimuthal cells; the paper's `P` counts polynomial
coefficients, with subscripts `phi` and `RZ` denoting azimuthal and
section directions. It gives:

| section cell width | non-aligned LDG | direct aligned dG | self-adjoint aligned dG |
|---|---:|---:|---:|
| 2/5 | 9.74e-4 | 3.05e-4 | 4.04e-4 |
| 2/10 | 7.22e-5 | 3.11e-5 | 6.22e-5 |
| 2/20 | 6.22e-6 | 3.91e-6 | 8.13e-6 |
| 2/40 | 6.44e-7 | 4.89e-7 | 1.01e-6 |
| 2/80 | 7.31e-8 | 7.75e-8 | 1.36e-7 |

All errors are dimensionless. These tabulated profile errors are not
`chi_perp,num/chi_perp`, since physical perpendicular diffusivity is zero.
Published-version table numbering has not been checked.

## Validation gates and descriptive fits

The figure reading is a reported comparison only. Its interval is not a published
uncertainty, and the prescribed analytic field is not the magnetic field evolved
in NIMROD. Validation requires the original-matrix normwise residual, isotropic L2
convergence of order `p+1`, and convergence of the discrete parallel gradient
applied to an interpolant of an analytic smooth function at order `p` in L2, for
each degree used. The last two refinement rates must be within `0.3` of the design
order. The gradient check uses the analytic scalar test function `u = exp(x+2*y)`,
with independent analytic partial derivatives `u_x = u` and `u_y = 2*u`, and `p+4`
quadrature points per direction. It exercises the same generated basis and
anisotropic integrand as assembly, against an analytic derivative rather than a
second generated operator. The field remains the normalized Sovinec field.

**Verified residual check:** for discrete matrix `A`, solution vector `t`,
and load `f` after boundary elimination, every solve reports the normwise residual `eta`:

`eta = ||A t - f||_inf / (||A||_inf ||t||_inf + ||f||_inf)`.

The matrix norm is the induced infinity norm, the maximum absolute row sum;
vector norms are maximum absolute values. The limit is `1e3 * epsilon_machine`,
where `epsilon_machine` is double-precision machine epsilon, approximately
`2.220446e-13` after multiplication by `1e3`. Residuals use the unfactored
matrix. An exact two-by-two oracle and equation rescaling independently
check this normalization. The six symmetric sine-eigenmode checks pass;
at 16 cells and `R = 1e6`, `eta = 8.24625e-17`. This replaces the
load-normalized residual, whose anisotropy dependence made its fixed limit
inappropriate. It changes the stated error measure, not its tolerance to
the same quantity.

**Verified isotropic convergence:** continuous tensor-product quadrilaterals
use the native GLL (Gauss-Lobatto-Legendre) nodal basis and `p+1` Gauss points
per direction for assembly. The L2 norm is integrated independently with
`p+4` points per direction against the analytic cosine-cosine solution.
All four meshes, with 2, 4, 8 and 16 elements per direction, are retained
in the [isotropic CSV](figures/fci_transport_pollution_isotropic.csv).
Both of the last two refinement rates must be within `0.3` of `p+1`.

| degree p | required L2 order | last measured order | result |
|---|---:|---:|---|
| 1 | 2 | 1.997795 | pass |
| 2 | 3 | 2.994973 | pass |
| 3 | 4 | 3.996267 | pass |
| 4 | 5 | 4.997632 | pass |

**Verified parallel-gradient convergence:** the [parallel CSV](figures/fci_transport_pollution_parallel.csv)
retains all four meshes for each degree. The final measured L2 orders for
`p=1,2,3,4` are `0.999004, 1.999378, 3.000193, 4.000875`; both final rates
pass the stated `p +/- 0.3` check for every degree.

**Descriptive bilinear fit:** conforming bilinear quadrilaterals on the
Sovinec thermal problem, with prescribed analytic
`B = ez cross grad(T_exact)` normalized at each quadrature point, are
measured at `h = 1/8, 1/16, 1/32, 1/64` and `R = 1e7, 1e8, 1e9, 1e10`.
The field is not aligned with the Cartesian mesh. At every point the signed
pollution exceeds ten times the absolute isotropic effective-diffusivity
error and is resolved above the subtraction-roundoff check. All original
normwise residuals pass; the maximum is `1.916e-16`.

A least-squares fit of all 16 logarithmic measurements, without filtering,
gives the following coefficients. An independent NumPy least-squares solve
in the plotting script reproduces the Fortran regression.
The logarithmic root-mean-square fit residual is `0.5685`, so this power law is not used for extrapolation.

| method / problem | a | b | C |
|---|---:|---:|---:|
| H1 p=1 quads / Sovinec, descriptive fit | 1.501624 | 0.841869 | 0.04538546 |

Here H1 denotes the continuous finite-element space with square-integrable
first derivatives.
For fixed `R` in increasing order, the fitted mesh exponents are
`0.580073, 1.412033, 1.946508, 2.067881`. For the four meshes in increasing
resolution, the anisotropy exponents are
`0.998922, 0.982561, 0.841581, 0.544410`. No single law holds in this window. The high-anisotropy mesh exponent
approaches two, while the anisotropy exponent decreases on finer meshes.
This is consistent with a regime change and locking-type saturation of the
conforming bilinear space; it is not evidence of an assembly defect.
The fit and all 16 rows are retained as reports.

## Direct cost to signed pollution threshold

The [complete measurements](figures/fci_transport_pollution_cost.csv) retain
all 67 rows from 24 ladders; the [summary CSV](figures/fci_transport_pollution_cost_summary.csv)
contains the first sampled mesh reaching `E <= 0.1`, or the final mesh at
the cap. Meshes start at four elements per direction and double until the
threshold, with the last mesh clipped to one million free unknowns. The
bilinear cap uses 1001 elements per direction and exactly one million
unknowns; its central value is interpolated between nodes. Higher-degree
quadrilaterals use the existing GLL basis. No fitted law determines a row.

The solver is FortSparse's native real symmetric positive-definite direct
solver, through `sparse_direct_factor_csc` in `fortfem_sparse_direct`; it factors a full,
exactly symmetric compressed sparse column matrix. Small-mesh solutions
at anisotropies `1`, `1e5` and `1e10` agree with the independent LAPACK
banded solver for all four degrees within relative infinity error `1e-6`.
All cost-run normwise residuals are below `5.95e-16`, against the unchanged
limit `2.220446e-13`. No iterative solver or solver tuning is introduced.

The retained test-side Günter symmetric stencil places the parallel flux at
quad centres and the unknowns on the dual cell grid, with its existing
five-point isotropic term. The bilinear FCI candidate uses those same
unknowns and isotropic term. There is **no distinct plane coordinate** in this
2D fixture: both endpoint maps interpolate the same four Cartesian cell
corners. In particular, there are no successive out-of-plane sections; the
Sovinec field has zero out-of-plane component.

At a quad centre `c`, the tracer solves `dX/ds = b(X)` for the position `X` at arc
lengths `+/-h/2`, where `s` is arc length and `c` the starting point, using 64
fourth-order Runge–Kutta (RK4) steps per direction. Its input is the exact
analytic normalized field `b = (ez cross grad(psi))/|grad(psi)|`, where `psi =
cos(pi*x)*cos(pi*y)` labels magnetic flux contours and `ez` is the unit vector
normal to the square, rather than an interpolated discrete field. The endpoint
locations are numerical RK4 approximations, not analytic maps. Each
parallel-gradient row is the difference of the two bilinear interpolation rows
divided by `h`; its `h^2`-weighted adjoint supplies the diffusion form. Thus its
links follow field lines although its unknowns remain on a Cartesian grid. Exact
trajectories preserve `psi`; interpolation need not preserve its nodal samples.
The circular fixture instead traces `(-y,x)/sqrt(x^2+y^2)`. This is a
two-dimensional support reconstruction, not a reproduction of a published
three-dimensional FCI plane discretization. Higher-order FCI is not added.

Times are wall seconds on a shared macOS arm64 machine, GNU Fortran 16.2, Release,
OpenBLAS, with one OpenMP thread and one Basic Linear Algebra Subprograms
(BLAS) thread per case. A CTest resource
lock prevents these timing cases from overlapping each other. Assembly
includes mesh tables, sparsity, tracing where needed, and matrix/load
construction. Solve includes factorization and back substitution; it excludes
the separate isotropic-baseline solve and residual measurement. Each row
is a single measurement, so microsecond differences are not a speed ranking.

The combined table below reports the Sovinec and centred circular Günter
2005 thresholds together with the manufactured temperature-error threshold.
The circular problem uses `T=1-r^3`, source `9*r`, and analytic Dirichlet data.
The two test-side stencils tie at nine unknowns on the signed central-value
pollution criterion for both problems at every sampled anisotropy. Nine
unknowns is the coarsest mesh tried, not a measured minimum resolution.

Negative values are retained, as the signed diagnostic requires. They do
not establish negative physical diffusion. This central-value criterion does
not bound total spatial error, and a nine-unknown success does not establish
accuracy for general loads or magnetic geometry. In particular, the existing
symmetric sine eigenmode is exceptionally well preserved on this uniform
grid. The CSV also retains the isotropic baseline and total central diagnostic.
These tables alone do not qualify a production method; the manufactured
temperature-error comparison below supplies a separate accuracy measure.

The [plotting script](figures/plot_fci_transport_cost.py) displays every row,
uses a signed logarithmic scale for pollution, and regenerates the summary.
It uses no absolute value, smoothing, fitted extrapolation or omitted mesh.

## Parallel consistency and manufactured validation

The original [parallel CSV](figures/fci_transport_pollution_parallel.csv)
checks only the conforming H1 candidates. The additional
[candidate gradient CSV](figures/fci_transport_pollution_candidate_parallel.csv)
applies each actual FCI or symmetric assembly stencil to nodal samples of
`u = exp(x+2*y)` and compares it to the independent analytic derivative
`b.grad(u) = (b_x+2*b_y)*exp(x+2*y)` at every cell centre, where `b_x` and
`b_y` are the Cartesian components of the unit field direction. The L2 error
uses cell areas `h^2`. On 4, 8, 16, 32, 64, 128 and 256 cells per direction,
the final measured orders are `0.993023` for FCI and `1.999740` for the
symmetric stencil; their final errors are `2.078403e-3` and `1.118234e-5`.
Both derivatives are nonzero and converge. The normalized field is singular
at the magnetic axis; no smooth-field global order theorem is inferred from
these measurements.

For the manufactured problem, keep this Sovinec field and `chi_perp=1`, but
set

`v = (cos(pi*x)*sin(pi*y), -sin(pi*x)*cos(pi*y))`,
`s = v_x^2+v_y^2`, and `T_exact = psi + 4*x*psi*s^2`.
Here `v_x` and `v_y` are components of the unnormalized field vector `v`;
`s` now denotes its squared magnitude, rather than tracing arc length, and
`T_exact` is the manufactured exact temperature. The source `Q_R` is the
heat source at anisotropy `R`.

The second term varies along field lines and vanishes on the boundary and
at the magnetic axis. Its L2 norm is `0.12998343504`, independently integrated
by the plotting script with 32, 64 and 128 Gauss points per direction. This
term is odd in `x`, whereas a flux function of `psi` is even in `x`; replacing
it by any flux function cannot remove that L2 component. Thus a flux-only
solution cannot meet the `0.01` error target. The source changes with `R` so that this same exact
temperature solves `-div([I+(R-1)*b*b^T] grad(T)) = Q_R`, where `I` is the
identity tensor and `b*b^T` the outer product of `b` with itself. FortSym generates
the analytic source from the flux; cancellation of the factor `s` before
differentiation gives a finite expression at the field zeros. An independent
fourth-order numerical divergence, using a separate temperature formula and
step `1e-3`, agrees with `Q_R/R` within `1.69e-8` at three off-axis points
for each tested anisotropy, against a fixed `2e-6` check. The axis source is
independently checked against `2*pi^2`.

All candidates use their existing assembly quadrature and the same native
FortSparse symmetric positive-definite (SPD) direct solver. Absolute L2
temperature error is independently overintegrated with `p+4` Gauss points per
direction; the area of the normalized square is one. For the computed temperature
`T_h`, the central error is `|T_h(0,0)-1|`. Nonzero Dirichlet values are included
in the circular error reconstruction; a nodal interpolation oracle verifies
second-order L2 convergence with those boundary values. No residual limit or
solver setting is relaxed. The maximum original-matrix normwise residual across
the 217 new refinement rows is `5.78e-16`, below the unchanged `2.220446e-13`
limit.

The [complete refinement CSV](figures/fci_transport_pollution_validation.csv)
retains every mesh, L2 error, central error, signed pollution where applicable,
residual and timing. Refinements run from 4 to 256 cells per direction; degree
four stops at 250 cells to respect the unknown cap. The unresolved manufactured
H1 degree-one cases continue through 512 and 1001 cells, and degree two at
`R=1e10` continues to 500 cells, so the threshold table does not extrapolate.
The two test-side stencils also retain all 4–256-cell Sovinec and circular
pollution measurements rather than stopping at the first passing mesh.

The following is the single combined cost table. Each entry is free unknowns
followed by assembly plus factorization/solve seconds; `N` denotes the
number of free unknowns. A missing threshold
reports the measured cap, not a fitted prediction. Sovinec and circular
pollution costs retain the preceding direct-cost measurements; manufactured
costs come from the new refinement measurements. The L2 threshold does not
also impose a central-error threshold; both errors remain in the full CSV.

| Candidate | R | Sovinec E <= 0.1: N; time (s) | Circular E <= 0.1: N; time (s) | Manufactured L2 <= 0.01: N; time (s) |
|---|---:|---:|---:|---:|
| H1 p=1 | 1e+06 | not reached at 1,000,000; 2.952 | — | not reached at 1,000,000; 2.622 |
| H1 p=1 | 1e+08 | not reached at 1,000,000; 2.979 | — | not reached at 1,000,000; 2.606 |
| H1 p=1 | 1e+10 | not reached at 1,000,000; 5.193 | — | not reached at 1,000,000; 2.655 |
| H1 p=2 | 1e+06 | 3,969; 0.01563 | — | 16,129; 0.05252 |
| H1 p=2 | 1e+08 | 16,129; 0.06044 | — | 65,025; 0.2048 |
| H1 p=2 | 1e+10 | 261,121; 1.16 | — | 998,001; 3.814 |
| H1 p=3 | 1e+06 | 121; 0.00079 | — | 2,209; 0.0122 |
| H1 p=3 | 1e+08 | 2,209; 0.01907 | — | 9,025; 0.04594 |
| H1 p=3 | 1e+10 | 9,025; 0.06451 | — | 36,481; 0.1955 |
| H1 p=4 | 1e+06 | 225; 0.003128 | — | 225; 0.002149 |
| H1 p=4 | 1e+08 | 225; 0.002829 | — | 225; 0.001902 |
| H1 p=4 | 1e+10 | 961; 0.01424 | — | 961; 0.008907 |
| Gunter symmetric | 1e+06 | 9; 1.5e-05 | 9; 1.3e-05 | 961; 0.001199 |
| Gunter symmetric | 1e+08 | 9; 1.4e-05 | 9; 1.5e-05 | 961; 0.001004 |
| Gunter symmetric | 1e+10 | 9; 1.4e-05 | 9; 1.4e-05 | 961; 0.001002 |
| Bilinear FCI | 1e+06 | 9; 0.000497 | 9; 0.000276 | 961; 0.02634 |
| Bilinear FCI | 1e+08 | 9; 0.00051 | 9; 0.000253 | 961; 0.02658 |
| Bilinear FCI | 1e+10 | 9; 0.000493 | 9; 0.000287 | 961; 0.02663 |

These results reject a completely absent parallel response on this manufactured
case. They do not establish accuracy for every load or field geometry, and the
2D FCI reconstruction does not qualify the production 3D plane operator.
At high anisotropy, the finest-mesh errors and signed pollution can be
nonmonotone even while original-matrix residuals pass. In particular, neither
negative pollution nor a small normwise residual is a forward-error bound.
The [combined threshold CSV](figures/fci_transport_pollution_thresholds.csv)
retains the table at full precision. The
[validation plotting script](figures/plot_fci_transport_validation.py)
plots every signed pollution value and every manufactured L2 and central error,
and regenerates the combined threshold CSV. No data are smoothed or discarded.

## Three-dimensional bilinear FCI locking

The dimensionless nested-surface fixture uses `(x,y)` in `[0,1]^2`, periodic `z`
of period one, `psi = sin(pi*x)*sin(pi*y)`, and magnetic field `B = (a*psi_y,
-a*psi_x, 1)`, with `a = 0.15*(1+psi)` and subscripts denoting partial
derivatives. It is divergence free, nowhere zero and tangent to `psi` surfaces.
Here `a` sets the field-line pitch, and `psi_x`, `psi_y` are partial derivatives
of `psi`. With unit perpendicular conductivity and parallel conductivity
`1/epsilon`, the exact temperature is `psi` for source `2*pi^2*psi` and zero outer
Dirichlet data. This translates and extends the Sovinec thermal source with shear
and a guide field; it does not reproduce NIMROD's evolved field.

Native fourth-order Runge–Kutta tracing, bilinear endpoint interpolation and
parallel-gradient assembly supplied connectivity between two canonical planes half
a period apart, with duplicate periodic planes identified. On 5, 9 and 17 section
cells per direction (72, 200 and 648 free unknowns), relative mass-norm errors
(discrete L2 errors weighted by nodal cell volumes and divided by the
exact-temperature norm) at `epsilon=1e-6` were `0.99706`, `0.97841` and `0.74913`.
On each fixed mesh the error tends to one as `epsilon` decreases toward zero,
reaching near-one errors at `epsilon=1e-13`; all 24 rows fail `error < 0.1`.
Perpendicular-only errors were `0.033558`, `0.010216`, `0.002851`. The [complete
CSV](figures/fci_transport_nested_accuracy.csv) retains every error, residual,
tracing check and timing; the [plot
source](figures/plot_fci_transport_nested_accuracy.py) displays all rows. The
original weak residual stays below `4.1e-10`, traced flux defect below `1.3e-10`,
and endpoint step-doubling difference below `2.6e-10`.

Bilinear interpolation preserves constants but not this nonlinear flux.
Positive interpolation paths connect every free node to the anchored wall,
even though the continuous interior field lines do not reach it. An independent
graph check and singular-value rank check found a trivial free-node parallel
kernel, so the fixed-mesh solution tends to zero as `epsilon` tends to zero.
The responsible interpolation is in `src/operators/fci_interpolation_map.f90`
(`build_fci_bilinear_interpolation_map_2d`) and its gradient use in
`src/operators/fci_parallel_operator.f90`. This is a spatial locking result,
not a Krylov failure or a conclusion about all FCI interpolation schemes.

## Constraint-space and averaging limitations

In a multiplier formulation an additional unknown enforces the weak
parallel-gradient constraint: the parallel gradient is tested against every
function in the multiplier space. A flux function is constant along magnetic field
lines. Coarsening a multiplier space cannot by itself isolate flux functions: its
kernel (the temperature fields whose tested parallel gradient vanishes) has
dimension at least the temperature dimension minus the multiplier dimension. The [reduced-space CSV](figures/fci_transport_reduced_space.csv)
records piecewise-constant and continuous-bilinear multipliers on a section grid
with twice the temperature spacing, at 21, 41 and 81 section cells. At 81 cells
the temperature space has dimension 12800, but the two kernels have dimensions
9600 and 9438, against a radial-resolution proxy of only 40 bins. Checkerboard
modes alternate sign between neighbouring grid nodes; the odd-plane modes also
have opposite signs on the two planes. Projected odd-plane checkerboards retain
98.91% and 98.51% of their norm. The proxy is not an exact flux-space dimension.
These measurements do not qualify this generic reduction route.

An explicit field-line-average projector, which would replace a temperature by its
average along each field line, was considered but not implemented. A small
parallel defect of its proposed range does not yield an epsilon-uniform
consistency estimate without controlling the parallel flux multiplier and the
complement. The energy estimate alone permits multiplier growth proportional to
`epsilon^(-1/2)`. Smooth resonant shear need not give a uniformly bounded
multiplier, and irrational field lines on a torus need not give a positive
parallel spectral gap. Thus neither accurate averaging nor kernel identification
alone proves the proposed projected penalty accurate. No averaging operator or
modified transport discretization was implemented.

## Dense algebraic oracle

The kernel/complement oracle remains only in `test/fci_transport_reference.f90`,
linked into its independent algebraic tests, never the installed library. The
kernel consists of fields annihilated by the parallel gradient; the complement
is its mass-orthogonal remainder.
It eliminates prescribed values and uses a mass-weighted singular value
decomposition of the support gradient to split kernel and complement before
scaling the weak equations. Rank is numerical, not certified; weak nonzero
modes below the rank threshold may be classified as kernel. Setup and solves
are dense, and storage is quadratic. This is not a scalable production solver
or an implementation of a published finite-element multiplier system.

Periodic and mixed closed/open line fixtures use independent Fourier
eigenvalues. The latter includes plate values 2 and 5. The
[accuracy CSV](figures/fci_transport_reference_accuracy.csv) and
[plot source](figures/plot_fci_transport_reference_accuracy.py) retain all
48 zero-mass-shift rows over `epsilon=1e-6` through `1e-13`. Maximum relative
mass errors are `1.035e-13` (periodic) and `4.232e-15` (mixed), below the
`2e-11` gate. Shifted and rescaled equations, nonuniform interpolation weights,
rank rejection and invalid inputs are also tested. These algebraic fixtures
do not establish spatial convergence on non-aligned magnetic surfaces.
At `epsilon=1e-13`, assembled sparse solves of the separate cost fixtures
had relative errors `0.24` to `1.45`, while the split oracle stayed below
`4.506e-14`. Their timings compare different work and accuracy, not a practical
solver crossover; the cost fixtures in
[`test_fci_transport_reference.f90`](../../test/test_fci_transport_reference.f90)
provide the corresponding algebraic comparison.

## Consequence for temperature transport

The tested two-dimensional Günter symmetric and FCI stencils show
anisotropy-independent sampled resolution requirements for the manufactured load.
High-order H1 controls the same error at modest resolution; degree four has a
small but measurable resolution increase with anisotropy, retained in the table,
so a strict uniform-error theorem is not claimed. The three-dimensional bilinear
plane operator locks on the non-aligned nested fixture. These benchmarks do not
qualify a production solver; island, plate, single-null and general
three-dimensional accuracy and preconditioner gates remain unqualified.

## Figure report on quadrilaterals

The existing reconstruction already uses quadrilaterals. The source, square,
temperature boundary, conforming bilinear space and four-point quadrature
match the stated thermal setup. The analytic field does not establish
equivalence to the paper's discrete magnetic field.

| case, h=1/8 and R=1e3 | free unknowns | central T | paper diagnostic |
|---|---:|---:|---:|
| prescribed analytic field, bilinear quads | 49 | 0.2378644074 | 3.204075806 |
| Sovinec manuscript Fig. 9 | -- | -- | about 6; reading [5,7] |

This row is retained without a pass/fail comparison to the figure reading.

## Reproduction

The [figure row](figures/fci_transport_pollution_published.csv),
[sine residuals](figures/fci_transport_pollution_checks.csv),
[isotropic errors](figures/fci_transport_pollution_isotropic.csv),
[scaling rows](figures/fci_transport_pollution_scaling_p1.csv),
[independent fit](figures/fci_transport_pollution_fits.csv), and
[plot source](figures/plot_fci_transport_pollution.py) retain all measured rows.
The plot has normalized units and no smoothed or discarded data.

```sh
export PKG_CONFIG_PATH=/opt/homebrew/opt/openblas/lib/pkgconfig
cmake -S . -B build/transport -G Ninja -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_Fortran_COMPILER=/opt/homebrew/bin/gfortran \
  -DFORTFEM_BUILD_TESTING=ON -DBLA_VENDOR=OpenBLAS \
  -DCMAKE_PREFIX_PATH=/opt/homebrew/opt/openblas
cmake --build build/transport -j 8
ctest --test-dir build/transport -j 8 -L '^fci$' -LE benchmark --output-on-failure
ctest --test-dir build/transport -j 8 -L '^fci$' --output-on-failure
ctest --test-dir build/transport -j 8 -L '^fast$' --output-on-failure
```

The fast CTest selection (`-L '^fast$'`, or FCI with `-LE benchmark`) includes
validation gates and one four-cell cost row for each of the six candidates.
These rows check the original-matrix residual and report cost; they do not
assert a pollution threshold. Full cost, scaling and refinement ladders carry
the `benchmark` label and remain registered and runnable. Bare CTest includes
benchmarks; use explicit labels to select the fast scope. The native fpm default
runs the diagnostic program, not its command-line ladders.

CTest sets one OpenMP and one OpenBLAS thread per FCI case; builds use eight jobs.
The timing cases share a CTest resource lock. Retained CSV times are single
measurements on a shared macOS arm64 machine with GNU Fortran 16.2,
Release/OpenBLAS, without uncertainty estimates. Direct ladders stop at signed `E <= 0.1` or one million unknowns; validation refinements retain every row,
including nonmonotone errors. No solver tolerance or residual gate is relaxed.

Plot sources read the retained CSVs by default, and accept output paths outside
the repository. For example:

```sh
MPLBACKEND=Agg python3 doc/design/figures/plot_fci_transport_cost.py \
  --output /tmp/sovinec-cost.png --circle-output /tmp/circle-cost.png \
  --summary-output /tmp/pollution-cost-summary.csv
MPLBACKEND=Agg python3 doc/design/figures/plot_fci_transport_validation.py \
  --output-dir /tmp/transport-validation \
  --summary-output /tmp/transport-thresholds.csv
```
