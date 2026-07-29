# Toroidal Poisson and Ampère: analytical versus BEM

This example evaluates the real \(n=2,m=1\) toroidal harmonic

\[
\Psi=\sqrt{\cosh\eta-\cos\theta}\,
P_{3/2}^{1}(\cosh\eta)\cos(2\theta)\cos(\phi)
\]

and the current-free magnetostatic field \(\mathbf H=-\nabla\Psi\). The
half-integer Legendre function and derivative come from Fortnum's
Hobson-normalized toroidal special functions.

The numerical path triangulates the constant-\(\eta\) torus with
\(\cosh\eta=2\), samples exact Dirichlet and Neumann data there, and evaluates
the three-dimensional Green representation in the exterior. Thus this is a
boundary-element representation test with an independent analytical oracle;
it does not claim to solve an unknown boundary density.

The program produces:

- `toroidal_trace_1d.png`: analytical and BEM Poisson/Ampère traces;
- `toroidal_surface_2d.png`: the analytical scalar mode in \((\theta,\phi)\);
- `toroidal_bem_error_2d.png`: BEM relative error over the exterior surface;
- `toroidal_helmholtz_1d.png`: outgoing Helmholtz point source, exact vs BEM;
- `toroidal_geometry_3d.png`: the physical toroidal interface;
- `toroidal_trace.csv`: reproducible source data for comparisons.
- `benchmark.txt`: mesh size, timings, and global relative errors.

Rendered outputs are generated in CI and are intentionally not committed.

The separated toroidal harmonic follows
[DLMF §14.19](https://dlmf.nist.gov/14.19). The boundary representation uses
the standard Laplace Green formula with outward surface normals.
