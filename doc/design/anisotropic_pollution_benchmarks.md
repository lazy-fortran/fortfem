---
title: Anisotropic transport pollution benchmarks
---

# Anisotropic transport pollution benchmarks

The symmetric scheme and bilinear flux-coordinate-independent (FCI) candidate
meet the signed pollution threshold on the smallest tested mesh, with nine
free unknowns, on both problems at all three anisotropies. Among continuous
finite elements, degree four reaches the threshold at 961 unknowns when the
anisotropy is `1e10`; degree one does not reach it at one million unknowns.
These are minima on the sampled refinement ladders, not optimal mesh searches.

The comparison uses directly measured signed pollution and cost; no fitted
power law is an accuracy gate or an extrapolation rule. The bilinear
quadrilateral measurements cross response regimes, so no single power law
holds throughout the retained window. Method selection remains separate
from these measurements.
The [accuracy target](fci_transport_solver.md#accepted-accuracy-target) is
numerical perpendicular transport well below physical perpendicular transport
at physical anisotropy, rather than epsilon-uniform accuracy for general loads
on a fixed mesh.

Here `T` is temperature, `Q` the heat source, `b` the unit field direction,
`chi_par` and `chi_perp` the parallel and perpendicular diffusivities, `p`
the element degree, and `h` the Cartesian element width. Coordinates,
temperatures and diffusivities below are normalized as in each source.
Write `R = chi_par/chi_perp` for the anisotropy. For the Sovinec test,
`chi_perp = 1`, and the reported paper diagnostic is `|1/T(0,0)-1|`.
It includes perpendicular discretization error. The isolated, signed pollution
measure used in the fit is instead

`E = chi_perp,num/chi_perp = 1/T_R(0,0) - 1/T_1(0,0)`.

Here `T_R` is the discrete solution at anisotropy `R`, and `T_1` its
isotropic baseline on the same mesh. No absolute value is applied to this
difference. A fit to `E = C R^b h^a` estimates the mesh exponent `a`,
anisotropy exponent `b`, and constant `C`; all quantities in this law are
dimensionless in the stated normalization.

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
`r = sqrt((x-x0)^2+(y-y0)^2)`, `T = 1-r^3`,
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
`f = Delta_z*tan(theta)/h`, and line step `Delta_s`.
Its linear-interpolation naive diffusivity is
`f*(1-f)*h^2/Delta_s^2` for unit parallel diffusivity. The support scheme
instead has fourth-order modal hyperdiffusion. Section 3.1 uses the
circular shell `0.1 <= rho <= 0.2` and radial sine modes, equation (43).
For zonal modes the measured decay is entirely numerical. Section 3.1.2,
Fig. 11 and equations (47)–(48), reports `h^2`/`h^4` for naive linear/cubic
interpolation, `h^4` for support-linear and an `A*h^8+B*h^4` contribution
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
Section V.B uses Neumann walls, `chi_par = 100`, `chi_perp = 0`, and
time step `0.001`. Section V.B.2 initializes the stationary profile
`f0 = 0.1*psi_p^2` and measures relative L2 profile loss.
Table VII at `t = 5` fixes `N_phi = 10`, `P_phi = 1`, `P_RZ = 3`;
the paper's `P` counts polynomial coefficients. It gives:

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

The figure reading is a reported comparison only. Its interval is not a
published uncertainty, and the prescribed analytic field is not the magnetic
field evolved in NIMROD. Validation requires the original-matrix normwise residual, isotropic L2
convergence of order `p+1`, and convergence of the discrete parallel gradient
applied to an interpolant of an analytic smooth function at order `p` in L2,
for each degree used. The last two refinement rates must be within `0.3` of
the design order. The gradient check uses `u = exp(x+2*y)`, with independent
analytic derivatives `u_x = u` and `u_y = 2*u`, and `p+4` quadrature points per
direction. It exercises the same generated basis and anisotropic integrand
as assembly, against an analytic derivative rather than a second generated
operator. The field remains the normalized Sovinec field.

**Verified residual check:** for discrete matrix `A`, solution vector `t`,
and load `f` after boundary elimination, every solve reports

`eta = ||A t - f||_inf / (||A||_inf ||t||_inf + ||f||_inf)`.

The matrix norm is the induced infinity norm, the maximum absolute row sum;
vector norms are maximum absolute values. The limit is `1e3 * epsilon`,
where `epsilon` is double-precision machine epsilon, approximately
`2.220446e-13` after multiplication by `1e3`. Residuals use the unfactored
matrix. An exact two-by-two oracle and equation rescaling independently
check this normalization. The six symmetric sine-eigenmode checks pass;
at 16 cells and `R = 1e6`, `eta = 8.24625e-17`. This replaces the former
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
unknowns and isotropic term. At each quad centre the existing RK4 tracer
follows the unit field by `+/-h/2`, with 64 steps in each direction. The
existing bilinear map and generated gradient form the endpoint difference
divided by `h`; its volume-weighted adjoint supplies the diffusion form.
This is a two-dimensional support-operator reconstruction with line spacing
`h`, not a reproduction of a published three-dimensional FCI table.
Higher-order FCI is not added.

Times are wall seconds on shared faepmac1, GNU Fortran 16.2, Release,
OpenBLAS, with one OpenMP and one BLAS thread per case. A CTest resource
lock prevents these timing cases from overlapping each other. Assembly
includes mesh tables, sparsity, tracing where needed, and matrix/load
construction. Solve includes factorization and back substitution; it excludes
the separate isotropic-baseline solve and residual measurement. Each row
is a single measurement, so microsecond differences are not a speed ranking.

| Sovinec candidate | R | free unknowns | signed E | setup (s) | factor + solve (s) |
|---|---:|---:|---:|---:|---:|
| H1 p=1 | 1e+06 | 1,000,000 (not reached) | 0.5088 | 1.185 | 1.767 |
| H1 p=1 | 1e+08 | 1,000,000 (not reached) | 34.627 | 1.19 | 1.789 |
| H1 p=1 | 1e+10 | 1,000,000 (not reached) | 3431.3 | 1.98 | 3.213 |
| H1 p=2 | 1e+06 | 3,969 | 0.04315 | 0.01005 | 0.005578 |
| H1 p=2 | 1e+08 | 16,129 | 0.088998 | 0.04307 | 0.01737 |
| H1 p=2 | 1e+10 | 261,121 | 0.066699 | 0.711 | 0.4492 |
| H1 p=3 | 1e+06 | 121 | 4.6729e-05 | 0.000707 | 8.3e-05 |
| H1 p=3 | 1e+08 | 2,209 | 0.048762 | 0.01418 | 0.004888 |
| H1 p=3 | 1e+10 | 9,025 | 0.063953 | 0.05406 | 0.01045 |
| H1 p=4 | 1e+06 | 225 | 0.00037415 | 0.002921 | 0.000207 |
| H1 p=4 | 1e+08 | 225 | 0.011337 | 0.002642 | 0.000187 |
| H1 p=4 | 1e+10 | 961 | 0.020701 | 0.01271 | 0.001531 |
| Gunter symmetric | 1e+06 | 9 | 4.5535e-11 | 1e-05 | 5e-06 |
| Gunter symmetric | 1e+08 | 9 | 2.1334e-09 | 1e-05 | 4e-06 |
| Gunter symmetric | 1e+10 | 9 | 2.6823e-07 | 1e-05 | 4e-06 |
| Bilinear FCI | 1e+06 | 9 | -0.0024323 | 0.000486 | 1.1e-05 |
| Bilinear FCI | 1e+08 | 9 | -0.0024323 | 0.000506 | 4e-06 |
| Bilinear FCI | 1e+10 | 9 | -0.0024323 | 0.000488 | 5e-06 |

The two candidates with the fewest unknowns on Sovinec are tied at nine.
They are compared below on the centred circular Günter 2005 problem,
with `T=1-r^3`, source `9*r`, and analytic Dirichlet data.

| Circular candidate | R | free unknowns | signed E | setup (s) | factor + solve (s) |
|---|---:|---:|---:|---:|---:|
| Gunter symmetric | 1e+06 | 9 | 3.2589e-11 | 9e-06 | 4e-06 |
| Gunter symmetric | 1e+08 | 9 | -1.079e-08 | 1e-05 | 5e-06 |
| Gunter symmetric | 1e+10 | 9 | -7.6268e-07 | 9e-06 | 5e-06 |
| Bilinear FCI | 1e+06 | 9 | -0.00023449 | 0.000271 | 5e-06 |
| Bilinear FCI | 1e+08 | 9 | -0.0002345 | 0.000248 | 5e-06 |
| Bilinear FCI | 1e+10 | 9 | -0.00023675 | 0.000283 | 4e-06 |

Negative values are retained, as the signed diagnostic requires. They do
not establish negative physical diffusion. This central-value criterion does
not bound total spatial error, and a nine-unknown success does not establish
accuracy for general loads or magnetic geometry. In particular, the existing
symmetric sine eigenmode is exceptionally well preserved on this uniform
grid. The CSV also retains the isotropic baseline and total central diagnostic.
No production method is selected or qualified by these tables.

The [plotting script](figures/plot_fci_transport_cost.py) displays every row,
uses a signed logarithmic scale for pollution, and regenerates the summary.
It uses no absolute value, smoothing, fitted extrapolation or omitted mesh.

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
cmake -S . -B build/u3g -G Ninja -DCMAKE_BUILD_TYPE=Release -DCMAKE_Fortran_COMPILER=/opt/homebrew/bin/gfortran -DFORTFEM_BUILD_TESTING=ON -DBLA_VENDOR=OpenBLAS -DCMAKE_PREFIX_PATH=/opt/homebrew/opt/openblas
cmake --build build/u3g -j 8 --target test_fci_transport_pollution
ctest --test-dir build/u3g -j 8 --stop-on-failure -R '^test_fci_transport_pollution$' -V
ctest --test-dir build/u3g -j 8 --stop-on-failure -R '^test_fci_transport_pollution_isotropic_gate$' -V
ctest --test-dir build/u3g -j 8 --stop-on-failure -R '^test_fci_transport_pollution_parallel_gate$' -V
ctest --test-dir build/u3g -j 8 --stop-on-failure -R '^test_fci_transport_pollution_scaling_p1$' -V
ctest --test-dir build/u3g -j 8 --stop-on-failure -R '^test_fci_transport_pollution_figure_report$' -V
ctest --test-dir build/u3g -j 8 --stop-on-failure -R '^test_fci_transport_pollution_sparse_gate$' -V
ctest --test-dir build/u3g -j 8 --stop-on-failure -R '^test_fci_transport_pollution_(cost_|circle_)' -V
python3 doc/design/figures/plot_fci_transport_cost.py --build-dir build/u3g --data /tmp/fci_transport_pollution_cost.csv --output /tmp/sovinec-cost.png --circle-output /tmp/circle-cost.png --summary-output /tmp/pollution-cost-summary.csv
python3 doc/design/figures/plot_fci_transport_pollution.py --output /tmp/transport-pollution.png --fit-output /tmp/transport-pollution-fits.csv
```

CTest runs with one OpenMP thread and one OpenBLAS thread per case; the native
build uses eight jobs. Timings in the CSV are diagnostics on shared faepmac1
with gfortran 16.2, single-run measurements without uncertainty estimates.
Direct cost ladders refine until `E <= 0.1` or the unknown cap. They retain
the signed diagnostic,
the isotropic baseline on each mesh, assembly time, and factorization plus
solve time. No solver-tolerance adjustment or solver-method research is part
of this comparison.
