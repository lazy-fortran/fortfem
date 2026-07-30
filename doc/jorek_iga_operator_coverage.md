# Structure-preserving IGA coverage for JOREK-type models

FortFEM's isogeometric path follows the tensor-product spline de Rham complex
and does not reproduce JOREK's historical scalar `G1` Bézier discretization.
The objective is compatible operator coverage with exact discrete differential
identities.

The operator inventory is based on:

- M. Hoelzl et al., *The JOREK non-linear extended MHD code and applications*,
  Nuclear Fusion 61 (2021), 065001,
  [doi:10.1088/1741-4326/abf99f](https://doi.org/10.1088/1741-4326/abf99f).
- E. Franck et al., *Energy conservation and numerical stability for the
  reduced MHD models of the non-linear JOREK code* (2015),
  [arXiv:1408.2099](https://arxiv.org/abs/1408.2099).
- D. Toshniwal and T. J. R. Hughes, *Isogeometric discrete differential
  forms: Non-uniform degrees, Bézier extraction, polar splines and flows on
  surfaces*, CMAME 376 (2021), 113576,
  [doi:10.1016/j.cma.2020.113576](https://doi.org/10.1016/j.cma.2020.113576).

## Implemented primitives

| Published-model requirement | FortFEM IGA primitive |
|---|---|
| Poloidal scalar mass and diffusion | Physical H1 mass/diffusion |
| Grad--Shafranov operator | Cylindrical `1/R` weighted H1 diffusion |
| Toroidal Fourier derivatives | Exact complex `i*n` multiplier |
| Toroidal part of scalar Laplacian | Mode-dependent `n^2/R` mass |
| Nonlinear poloidal Poisson brackets | Energy-skew Galerkin bracket and exact analytical JVP |
| Nonlinear toroidal mode coupling | Exact retained-mode convolution `p+q=n` |
| Scalar coefficient products | Coefficient-weighted H1 Galerkin mass |
| Magnetic-flux evolution | Cylindrical weak residual and analytical JVP |
| Density and pressure transport | No-parallel cylindrical residual and analytical JVP |
| Ideal poloidal flux subflow | Cayley/implicit-midpoint Poisson propagator |
| Magnetic/vector differential sequence | H1--H(curl)--H(div)--L2 incidence |
| Curl--curl and div--div energies | Physical Piola-mapped vector forms |
| Patch interfaces | Orientation-aware 2D and 3D quotient complexes |
| Magnetic-axis regularity | Type-1 polar H1/Hcurl/L2 extraction, exact incidence, and physical Hodge operators |

The nonlinear bracket convolution is deliberately truncated only after exact
integer mode addition. It therefore cannot alias a discarded triad into an
unrelated retained harmonic. Its analytical directional derivative retains
both bilinear product-rule terms. Deterministic random trials verify the
quadratic skew invariant, and an independent central difference verifies the
JVP for complex multi-mode fields.

The no-parallel thermodynamic block implements equations (26)--(27) of Franck
et al. as `R[q,u] + 2 gamma_q q dZ(u)`, with `gamma_q=1` for density and the
physical adiabatic index for pressure. Both products retain exact Fourier
convolution. An affine manufactured field checks the compressibility factor,
an independent cylindrical test function verifies mass conservation, and a
complex randomized central difference verifies the complete product-rule JVP.

The magnetic-axis extraction maps satisfy both commuting diagrams against an
independently assembled periodic tensor-product cell complex. FortSym proves
the generic polar-cap and regular-cell identities exactly. Sparse physical
operators use the Galerkin restriction `E A E^T`; no dense pseudoinverse is
formed. Direct radial--periodic-angular quadrature assembles the physical polar
H1 mass and diffusion, covariant-Piola H(curl) mass and curl--curl, and density
mapped L2 mass. Its constant mass reproduces the exact area of a regular-polygon
polar map, while its stiffness annihilates constants. Deterministic random
trials independently verify the H1-gradient/H(curl)-mass and
H(curl)-curl/L2-mass energy identities. The block layout is also cross-checked
during development against the MIT-licensed STRUPHY polar implementation,
while FortFEM retains its own Fortran implementation and tests.

## Time discretization

The ideal poloidal flux subflow

`dt(psi) = R [psi,u]`

with fixed electric potential is advanced by the implicit midpoint/Cayley map.
Because the assembled bracket is skew, this map preserves the discrete mass
norm and is exactly time reversible. Repeated steps retain one FortSparse
factorization.

FortNum provides symmetric Strang and fourth-order Yoshida compositions for
combining such model propagators. This follows the propagator architecture of
F. Holderied, S. Possanner, and X. Wang,
*MHD-kinetic hybrid code based on structure-preserving finite elements with
particles-in-cell*, JCP 433 (2021), 110143,
[doi:10.1016/j.jcp.2021.110143](https://doi.org/10.1016/j.jcp.2021.110143).

The fourth-order composition has negative substeps and is restricted to
reversible Hamiltonian pieces. Resistive and viscous pieces must be composed
with contractive dissipative propagators instead. Space--time FEM is not part
of the current implementation.

## Remaining coupled-model work

- Complete the residual and Newton/Jacobian block layout beyond the implemented
  magnetic-flux and no-parallel density/pressure equations, for electric
  potential, vorticity, and parallel velocity.
- Add the parallel-flow products to density/pressure transport and the
  parallel-gradient blocks.
- Add resistive, viscous, thermal, and source terms with coefficient fields.
- Generalize pairwise patch quotients to arbitrary patch-interface graphs.
- Replace the remaining dense tensor result workspace in physical magnetic-axis
  H(curl)/L2 assembly with direct CSC accumulation, and publish scaling
  benchmarks. Quadrature already traverses only active local support.
- Verify a nonlinear manufactured equilibrium and an energy-transfer problem
  against the published continuous invariants.
- Benchmark cached harmonic-bracket assembly and nonlinear residual evaluation.

These remaining items are model integration work; the single-patch physical
de Rham forms and pairwise conforming topology are already executable.
