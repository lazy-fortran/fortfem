# Toroidal-harmonic contract

`fortfem_toroidal_harmonics` exposes the pinned FortNum analytical functions

\[
P_{n-1/2}^{m}(x),\quad Q_{n-1/2}^{m}(x),\qquad x>1,
\]

and their derivatives with respect to `x`. The functions use the
Hobson normalization and nonnegative degree and order indices. In toroidal
coordinates, `x = cosh(eta)`, so these functions are the radial factors for
separated Laplace solutions on an exact torus. The P/Q names are deliberately
kept distinct from Fourier modes and from a toroidal-coordinate map.

FortNum owns the hypergeometric definitions, order recurrence, derivative
recurrence, and moderate-degree domain policy. Its zero-order Q branch uses a
DLMF zero-balanced continuation near `x = 1` and obtains associated orders from
the analytical derivative before applying the generated order recurrence.
FortFEM only provides the public adapter used by analytical DtN/BEM/FEM
fixtures. The adapter test checks independent high-precision values for both
branches, derivative values, near-cut finiteness, and NaN rejection outside
`x > 1`.

The normalization and provenance are recorded in the [FortNum
API](https://github.com/lazy-fortran/fortnum/blob/d8be030/docs/api.md), with
the half-integer specialization specified by [DLMF
14.19](https://dlmf.nist.gov/14.19).
