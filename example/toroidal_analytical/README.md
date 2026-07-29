# Toroidal Poisson and Ampère analytical reference

This example evaluates the real \(n=2,m=1\) toroidal harmonic

\[
\Psi=\sqrt{\cosh\eta-\cos\theta}\,
P_{3/2}^{1}(\cosh\eta)\cos(2\theta)\cos(\phi)
\]

and the current-free magnetostatic field \(\mathbf H=-\nabla\Psi\) on the
constant-\(\eta\) torus with \(\cosh\eta=2\). The half-integer Legendre
function and derivative come from Fortnum's Hobson-normalized toroidal
special functions.

The program produces:

- `toroidal_trace_1d.png`: scalar trace and three orthonormal field components;
- `toroidal_surface_2d.png`: the scalar mode in \((\theta,\phi)\);
- `toroidal_geometry_3d.png`: the physical toroidal interface;
- `toroidal_trace.csv`: reproducible source data for comparisons.

Rendered outputs are generated in CI and are intentionally not committed.
The same figures will gain FEM, FEM–BEM, DtN, PML (Helmholtz only), and
far-boundary overlays as those numerical examples are completed.
