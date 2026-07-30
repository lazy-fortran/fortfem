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

## Implemented primitives

| Published-model requirement | FortFEM IGA primitive |
|---|---|
| Poloidal scalar mass and diffusion | Physical H1 mass/diffusion |
| Grad--Shafranov operator | Cylindrical `1/R` weighted H1 diffusion |
| Toroidal Fourier derivatives | Exact complex `i*n` multiplier |
| Toroidal part of scalar Laplacian | Mode-dependent `n^2/R` mass |
| Nonlinear poloidal Poisson brackets | Energy-skew Galerkin bracket |
| Nonlinear toroidal mode coupling | Exact retained-mode convolution `p+q=n` |
| Magnetic/vector differential sequence | H1--H(curl)--H(div)--L2 incidence |
| Curl--curl and div--div energies | Physical Piola-mapped vector forms |
| Patch interfaces | Orientation-aware 2D and 3D quotient complexes |

The nonlinear bracket convolution is deliberately truncated only after exact
integer mode addition. It therefore cannot alias a discarded triad into an
unrelated retained harmonic.

## Remaining coupled-model work

- Assemble the complete reduced-MHD residual and Newton/Jacobian block layout
  for flux, electric potential, vorticity, density, temperature, and parallel
  velocity.
- Add conservative density/pressure transport and parallel-gradient blocks.
- Add resistive, viscous, thermal, and source terms with coefficient fields.
- Generalize pairwise patch quotients to arbitrary patch-interface graphs.
- Verify a nonlinear manufactured equilibrium and an energy-transfer problem
  against the published continuous invariants.
- Benchmark cached harmonic-bracket assembly and nonlinear residual evaluation.

These remaining items are model integration work; the single-patch physical
de Rham forms and pairwise conforming topology are already executable.
