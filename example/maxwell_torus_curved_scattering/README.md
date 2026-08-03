# Exact-curved torus Maxwell scattering

This example solves time-harmonic perfect-electric-conductor scattering on an
exact parametric torus. The numerical trace space uses surface-Piola RWG
functions, a barycentrically refined Buffa--Christiansen dual space, and the
nondegenerate RWG--RBC pairing. Radial-Duffy quadrature treats coincident
Green kernels, adaptive product quadrature treats adjacent panels, and a
one-sided extrapolation supplies the MFIE trace.

The resonance-safe CFIE composes the real-wavenumber EFIE with a coercive
imaginary-wavenumber BC operator through the discrete dual mass matrix. This
is the range-aware product algebra described by Scroggs, Betcke, Burman,
Smigaj, and van 't Wout, and the purely imaginary regularizer follows the
operator-preconditioned formulation of Le and Cools. The dual trace space
comes from the Buffa--Christiansen complex:

- [Scroggs et al., *Software frameworks for integral equations in
  electromagnetic scattering based on Calderón
  identities*](https://arxiv.org/abs/1703.10900)
- [Le and Cools, *An operator preconditioned combined field integral equation
  for electromagnetic scattering*](https://doi.org/10.1137/23M1581674)
- [Buffa and Christiansen, *A dual finite element complex on the barycentric
  refinement*](https://doi.org/10.1090/S0025-5718-07-01965-5)

Two incident plane waves share one assembled and factored CFIE. Their
cross-observed far fields must satisfy Lorentz reciprocity; this independent
physical identity is the example's accuracy gate. Degree two near-trace
quadrature fails that gate, while degree three passes and is recorded in the
benchmark output.

Generated, uncommitted artifacts are:

- `maxwell_torus_solution_2d.png`: reconstructed scattered magnetic-field
  magnitude in a physical (x)-(z) slice through the torus hole, with arrows
  showing the real instantaneous (H_x,H_z) vector components of the solved
  nonzero scattering field (the first gallery plot is therefore the solution,
  not a convergence diagnostic);
- `scattered_field.csv`: the same physical slice samples, including the
  reconstructed vector components and complex-field magnitude;
- `maxwell_torus_rcs_1d.png`: equatorial bistatic radar-cross-section cut;
- `maxwell_torus_rcs_2d.png`: angular radar-cross-section map;
- `maxwell_torus_rcs_3d.png`: normalized three-dimensional radiation surface;
- `maxwell_torus_dtn_1d.png`: weak RWG/RBC Maxwell DtN response on the same
  exact-curved torus surface;
- `rcs_trace.csv`: reproducible equatorial samples;
- `benchmark.txt`: unknown counts, CFIE, weak-DtN, and field-reconstruction
  timings, quadrature, reciprocity error, DtN response norm, and peak radar
  cross section.

`provenance.json` records the method references and makes clear that this is a
solver-generated torus fixture rather than a redistribution of external
scattering data.  Setting `MAXWELL_TORUS_FAST=1` selects a two-by-two curved
surface mesh for the independent ten-second output/data oracle; the reference
gallery keeps the three-by-three mesh and the same physical operator.

GitHub Actions generates these files and publishes them in the example
gallery. No rendered image is checked into the repository.
