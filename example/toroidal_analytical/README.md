# Toroidal Poisson and Ampère: analytical versus BEM

This example evaluates the real \(n=2,m=1\) toroidal harmonic

\[
\Psi=\sqrt{\cosh\eta-\cos\theta}\,
P_{3/2}^{1}(\cosh\eta)\cos(2\theta)\cos(\phi)
\]

and the current-free magnetostatic field \(\mathbf H=-\nabla\Psi\). The
half-integer Legendre function and derivative come from Fortnum's
Hobson-normalized toroidal special functions.

The numerical path uses exact parametric panels on the constant-\(\eta\)
torus with \(\cosh\eta=2\). It samples only Dirichlet data, solves the dense
Galerkin \(V\phi=(K-\frac12M)u\) system for the unknown P0 Neumann trace, and
then reconstructs the exterior field. The Ampère field is the numerical
gradient of that reconstructed scalar potential. The Helmholtz comparison
uses the corresponding complex curved-torus DtN solve with Laplace
singularity subtraction. Analytical toroidal harmonics and an outgoing point
source are independent oracles.

The program produces:

- `toroidal_trace_1d.png`: analytical and BEM Poisson/Ampère traces;
- `toroidal_surface_2d.png`: the analytical scalar mode in \((\theta,\phi)\);
- `toroidal_bem_error_2d.png`: BEM relative error over the exterior surface;
- `toroidal_helmholtz_1d.png`: outgoing Helmholtz point source, exact vs BEM;
- `toroidal_geometry_3d.png`: the physical toroidal interface;
- `toroidal_trace.csv`: reproducible source data for comparisons.
- `benchmark.txt`: mesh size, separate DtN-solve/evaluation timings, and
  global relative errors.

Rendered outputs are generated in CI and are intentionally not committed.

The separated toroidal harmonic follows
[DLMF §14.19](https://dlmf.nist.gov/14.19). The boundary representation uses
the standard Laplace/Helmholtz Green formulas with outward surface normals.
