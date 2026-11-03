---
title: Anisotropic transport pollution benchmarks
---

# Anisotropic transport pollution benchmarks

**Verified outcome:** the first literature-reproduction screen fails.
No discretization is selected, and no candidate refinement or anisotropy
ladder is qualified. The [accuracy target](fci_transport_solver.md#accepted-accuracy-target)
is numerical perpendicular transport well below physical perpendicular
transport at physical anisotropy, rather than epsilon-uniform accuracy for
general loads on a fixed mesh.

Here `T` is temperature, `Q` the heat source, `b` the unit field direction,
`chi_par` and `chi_perp` the parallel and perpendicular diffusivities, and
`h` the Cartesian element width. Coordinates and diffusivities below are
normalized as in each source. For the Sovinec test, `chi_perp = 1`, so
`|1/T(0,0)-1|` is the effective perpendicular diffusivity error. It includes
perpendicular discretization error; it is not by itself an isolated measure
of pollution from the parallel operator. A subsequent comparison must retain
the isotropic baseline and report the signed difference of the effective
diffusivities as well as the paper's absolute diagnostic.

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

## Reproduction failure

**Verified:** the test-side reconstruction uses continuous analytic
`B = ez cross grad(T_exact)`, normalized at each quadrature point, rather
than solving the driven magnetic problem. The field definition is in
`tools/codegen/app/gen_anisotropic_pollution_products.f90:27`; its use is
in `test/test_fci_transport_pollution.f90:79`. The source, square,
temperature boundary, conforming bilinear space and four-point quadrature
match the stated thermal problem. Matching those ingredients does not
establish equivalence to the paper's discrete magnetic field.

| case | free unknowns | central T | effective error | original weak residual |
|---|---:|---:|---:|---:|
| prescribed analytic field, bilinear, h=1/8, ratio=1e3 | 49 | 0.2378644074 | 3.204075806 | 1.3970e-13 |
| Sovinec manuscript Fig. 9, same h and ratio | — | — | approximately 6; reading interval [5,7] | — |

**Hypothesis:** the different magnetic representation contributes to the
discrepancy. No experiment isolates that contribution, and the result is
not a demonstrated defect in FortFEM's tensor basis or FCI operator.
Origin/main has no corresponding pollution-reproduction test; its FCI map
and parallel operator are unchanged on this branch.

**Verified independent check and limitation:** the symmetric-stencil sine
problem has central value
`[pi*h/(2*sin(pi*h/2))]^2`, independent of anisotropy in exact arithmetic.
All six forward-value checks pass. At 16 cells and ratio `1e6`, however,
the double-precision band solve's original residual is `2.3934e-9`, failing
the unchanged `1e-9` check. The residual is recomputed from the unfactored
band matrix in `test/test_fci_transport_pollution.f90:250`.
Roundoff amplified by anisotropy is a hypothesis, not an established cause.
This failure is retained as well. No threshold was relaxed.

| requested candidate | available foundation | qualified measurements |
|---|---|---|
| bilinear FCI | native 2D bilinear maps and support gradient | none |
| higher-order FCI | native cubic 1D maps; no native bicubic 2D builder found | none |
| conforming H1, p=1–4 | arbitrary-degree GLL line tables, Gauss quadrature | first p=1 reconstruction fails reproduction |
| Günter symmetric | test-side cell-centred stencil | sine oracle only; circular benchmark not run |

There is no measured pollution exponent, cost comparison at accepted
accuracy, or unknown count satisfying the ratio-`1e10` selection criterion.
No Krylov or option-B implementation was attempted.

## Reproduction

The [reproduction row](figures/fci_transport_pollution_published.csv),
[sine-oracle residuals](figures/fci_transport_pollution_checks.csv) and
[plot source](figures/plot_fci_transport_pollution.py) retain the failed
checks without filtering. The plot labels every quantity as dimensionless.
The graphical interval is intentionally a conservative rejection screen;
even agreement would not by itself validate the full harness to an accuracy
that the source does not state.

```sh
cmake -S . -B build/u3e -G Ninja -DCMAKE_BUILD_TYPE=Release -DCMAKE_Fortran_COMPILER=/opt/homebrew/bin/gfortran -DFORTFEM_BUILD_TESTING=ON -DBLA_VENDOR=OpenBLAS -DCMAKE_PREFIX_PATH=/opt/homebrew/opt/openblas
cmake --build build/u3e -j 8 --target test_fci_transport_pollution
ctest --test-dir build/u3e -j 8 -R '^test_fci_transport_pollution(_published_gate)?$' -V
python3 doc/design/figures/plot_fci_transport_pollution.py --output /tmp/transport-pollution.png
```

The two registered tests fail as described above. The new fixtures use only
analytic data. Timings in the CSV are diagnostic: on shared faepmac1 with
gfortran 16.2 and one BLAS thread the 49-unknown setup took 45 microseconds
and factorization/solve 16 microseconds, or 0.327 microseconds per unknown.
These are not matched-accuracy performance evidence. The next step is to
reconcile the published discretized field and thermal operator, then repeat
this reproduction gate before attempting any candidate ladder.
