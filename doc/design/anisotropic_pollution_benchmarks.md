---
title: Anisotropic transport pollution benchmarks
---

# Anisotropic transport pollution benchmarks

**Verified outcome:** isotropic assembly passes its convergence gate, but the
first conforming pollution-scaling gate fails. The bilinear quadrilateral fit
gives a mesh exponent of 1.502 rather than the required `2 +/- 0.3`.
The candidate comparison stops here; no method is selected or qualified.
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

## Validation gates and first failure

The figure reading is a reported comparison only. Its interval is not a
published uncertainty, and the prescribed analytic field is not the magnetic
field evolved in NIMROD. Validation instead requires the theoretical/published
pollution exponents within `+/- 0.3` over at least four mesh widths and four
anisotropies in a pollution-dominated regime, plus isotropic L2 convergence
of order `p+1` for every degree used. For naive conforming elements the
required exponents are `a = 2p`, `b = 1`; the support-operator expectations
and plane-spacing dependence are recorded above. None of these exponents
or tolerances is fitted to the measurements.

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

**Verified first failing gate:** conforming bilinear quadrilaterals on the
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
The logarithmic root-mean-square fit residual is `0.5685`, so the failed
power law is not used for extrapolation.

| method / problem | a | b | C | unknowns for E <= 0.1 at R=1e10 |
|---|---:|---:|---:|---|
| H1 p=1 quads / Sovinec, failed validation fit | 1.501624 | 0.841869 | 0.04538546 | unavailable: scaling gate fails |
| H1 p=2--4 / either problem | not measured | not measured | not measured | unavailable |
| bilinear FCI / either problem | not measured | not measured | not measured | unavailable |
| higher-order FCI / either problem | not measured | not measured | not measured | unavailable |
| Guenter symmetric / either problem | not measured | not measured | not measured | unavailable |

Here H1 denotes the continuous finite-element space with square-integrable
first derivatives, and FCI denotes flux-coordinate-independent interpolation.
For fixed `R` in increasing order, the fitted mesh exponents are
`0.580073, 1.412033, 1.946508, 2.067881`. For the four meshes in increasing
resolution, the anisotropy exponents are
`0.998922, 0.982561, 0.841581, 0.544410`. Thus even the aggregate mesh
exponent is outside `[1.7, 2.3]`, and the per-anisotropy gate fails as well.

The failure is enforced at
`test/test_fci_transport_pollution.f90:239`; the signed diagnostic is
formed at line 213. The prescribed field comes from
`tools/codegen/app/gen_anisotropic_pollution_products.f90:27` and is used
by the quadrature assembly in that test. This identifies the failing
observable and exercised discretization, not a proved code defect.
Origin/main contains no corresponding pollution test or existing fix.
**Hypothesis:** the chosen parameter window mixes distinct response regimes;
the variation of marginal exponents is consistent with this explanation,
but it does not prove it. Field normalization, quadrature effects and
anisotropic solution-shape changes have not been isolated. The isotropic
order check alone does not establish anisotropic consistency.

No higher-degree pollution gate, Stegmeir support-operator gate, Guenter
circular problem, four-candidate ladder, or solver work follows this failure.
The symmetric scheme has only its retained sine oracle; no circular result
is implied. No valid unknown count at `R = 1e10` or matched-accuracy cost
comparison is established.

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
[failed scaling rows](figures/fci_transport_pollution_scaling_p1.csv),
[independent fit](figures/fci_transport_pollution_fits.csv), and
[plot source](figures/plot_fci_transport_pollution.py) retain all measured rows.
The plot has normalized units and no smoothed or discarded data.

```sh
cmake -S . -B build/u3f -G Ninja -DCMAKE_BUILD_TYPE=Release -DCMAKE_Fortran_COMPILER=/opt/homebrew/bin/gfortran -DFORTFEM_BUILD_TESTING=ON -DBLA_VENDOR=OpenBLAS -DCMAKE_PREFIX_PATH=/opt/homebrew/opt/openblas
cmake --build build/u3f -j 8 --target test_fci_transport_pollution
ctest --test-dir build/u3f -j 8 -R '^test_fci_transport_pollution$' -V
ctest --test-dir build/u3f -j 8 -R '^test_fci_transport_pollution_isotropic_gate$' -V
ctest --test-dir build/u3f -j 8 -R '^test_fci_transport_pollution_scaling_p1$' -V
# Stop here: the scaling command fails. The figure report is not a validation gate.
ctest --test-dir build/u3f -j 8 -R '^test_fci_transport_pollution_figure_report$' -V
python3 doc/design/figures/plot_fci_transport_pollution.py --output /tmp/transport-pollution.png --fit-output /tmp/transport-pollution-fits.csv
```

CTest runs with one OpenMP thread and one OpenBLAS thread per case; the native
build uses eight jobs. Timings in the CSV are diagnostics on shared faepmac1
with gfortran 16.2, not matched-accuracy performance evidence. The next step
is to explain the missing common scaling regime or repair a demonstrated
anisotropic discretization defect before restarting the unchanged gates.
