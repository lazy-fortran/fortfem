# Spherical-harmonic contract

`fortfem_spherical_harmonics` exposes the pinned FortNum implementation of
standard orthonormal complex spherical harmonics. The contract uses the
Condon--Shortley phase and accepts

\[
  Y_l^m(\theta,\phi),\qquad 0\leq\theta\leq\pi,
\]

with periodic `phi`. The module also exposes analytic derivatives with
respect to `theta` and `phi`. Values at the poles are defined; angular
derivatives at a pole are intentionally undefined and return a quiet NaN,
because the azimuthal coordinate is singular there.

`spherical_harmonic_product_coefficient` returns the Gaunt coefficient for
expanding a product of two harmonics into the orthonormal basis. It applies
the triangle, parity, and azimuthal selection rules and uses the finite
Wigner-3j Racah sum owned by FortNum.

FortNum owns the normalization and recurrence implementation. FortFEM only
provides this stable public adapter so boundary operators, Fourier-FEM, and
toroidal examples consume the same pinned special-function revision. The
independent contract test checks the closed-form \(Y_0^0\) value through
`fortfem_api`, rather than checking repository wiring or a generated file.

The provenance and current degree/continuation limitations are recorded in
the [FortNum API](https://github.com/lazy-fortran/fortnum/blob/9bb6e7b/docs/api.md)
and the [DLMF spherical-harmonic definitions](https://dlmf.nist.gov/14.30).
