# Verified affine RT0 audit

`fortfem_verified_rt0_majorant` independently audits a candidate recovered flux
for scalar positive diffusion on a conforming, nonoverlapping affine triangle
mesh with homogeneous Dirichlet boundary data. The caller supplies an interval
extension of the coefficient, source and discrete gradient on physical boxes,
and an outward-rounded Friedrichs constant. A sampled callback is insufficient.

Only the recovered flux's zeroth edge moments are used. The first occurrence
of each global edge selects one exactly represented real moment; neighbouring
cells use opposite signs of that same number. The canonical reference RT0
fields are `(xi, eta-1)`, `(xi, eta)` and `(xi-1, eta)`. Their normal integrals
are the unit vectors for the three counterclockwise edges. The exact
contravariant Piola map therefore gives an H(div) field: neighbouring normal
traces cancel, including when input RTp coefficients originally disagree due
to floating-point solves or accumulation. This construction bypasses the
floating-point RTp moment-matrix inverse.

For this field sigma and any conforming candidate u_h in H1_0, integration by
parts, Cauchy--Schwarz, and `||v|| <= C_F ||v||_A` prove

```
||u-u_h||_A <= ||A^(-1/2)(sigma-A grad u_h)|| + C_F ||f+div sigma||.
```

The code integrates both squared norms with interval range enclosure. At level
L the reference triangle is partitioned into `4**L` triangles of exactly
representable dyadic area `1/(2*4**L)`. A child's coordinate bounding box
contains it; the generated interval Piola and integrand leaves enclose the
integrands throughout that box. Multiplying each range by its area encloses
the child's exact integral, without a quadrature-exactness assumption. Interval
summation and square root enclose arithmetic as well. Refinement narrows these
enclosures, but no requested width is guaranteed. There is no source-projection
or local-solve accuracy assumption.

This is a separate RT0 audit. It does not certify the narrower sampled RTp
estimate. It does not certify geometry mappings, curved cells, nonzero
Dirichlet data, discontinuous discrete states, arbitrary callbacks, or a
pointwise error band. Source and geometry modelling error stay separate.
The strict interval arithmetic execution requirements of FortNum apply:
IEEE arithmetic, no unsafe reassociation or unaccounted contraction, and
supported finite intervals. Invalid callback results or nonpositive
coefficient/determinant enclosures return failure rather than a certificate.

`test_verified_rt0_majorant` checks polynomial and rational integral
containment against independent antiderivatives, shared-edge repair, refinement
and rejection of a coefficient interval crossing zero.
