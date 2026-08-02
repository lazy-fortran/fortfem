---
title: Piola-enriched vector composition
---

# Piola-enriched vector composition

`fortfem_piola_enriched_vector` composes a standard affine or mapped-element
Piola transform with a pointwise vector enrichment.  It is geometry- and
application-neutral, so the same batched contract serves simplex FEEC,
multipatch IGA, and cut XFEM/XIGA callers.

For a reference vector `b_hat`, Jacobian `J`, and activation `a`, the two
supported transforms are

```
covariant H(curl):       b = J^(-T) b_hat,
contravariant H(div):    b = J b_hat / det(J),
enriched value:          e = a b.
```

The arrays use `(component, point, basis)` layout for vector values and
`(point, basis)` for activation.  Jacobians and their directions use
`(component, component, point)` layout.  Only orientation-preserving,
finite, nonsingular 2D and 3D Jacobians are accepted.

```fortran
call evaluate_piola_enriched_vector_values( &
    PIOLA_COVARIANT, jacobians, reference_values, activation, physical_values, &
    status)
call evaluate_piola_enriched_vector_values_jvp( &
    map_kind, jacobians, reference_values, activation, jacobians_dot, &
    reference_values_dot, activation_dot, physical_values_dot, status)
call evaluate_piola_enriched_vector_values_vjp( &
    map_kind, jacobians, reference_values, activation, physical_values_bar, &
    jacobians_bar, reference_values_bar, activation_bar, status)
call evaluate_piola_enriched_vector_differential_3d(...)
call evaluate_piola_enriched_vector_differential_3d_jvp(...)
call evaluate_piola_enriched_vector_differential_3d_vjp(...)
call evaluate_piola_enriched_vector_differential_2d(...)
call evaluate_piola_enriched_vector_differential_2d_jvp(...)
call evaluate_piola_enriched_vector_differential_2d_vjp(...)
```

The JVP differentiates the inverse/determinant Piola factors, reference
coefficients, and activation.  The VJP uses FortNum’s FortSym-derived inverse
and determinant reverse kernels and returns the geometry, reference, and
activation cotangents.  This makes shape sensitivity and enrichment
conditioning composable without duplicating covariant/contravariant formulas
in IGA or XFEM clients.

For 2D and 3D affine cells, the differential contract makes the de Rham effect
explicit.  A covariant H(curl) map reports

```text
curl(a b) = a J curl_hat(b)/det(J) + grad(a) x b,
```

in 3D (with the corresponding scalar rotated-gradient term in 2D), while a
contravariant H(div) map reports

```text
div(a b) = a div_hat(b)/det(J) + grad(a) . b.
```

The value, JVP, and VJP APIs differentiate the Jacobian, reference values,
reference curl/divergence, activation, and physical activation gradient.  The
unused differential is returned as zero for the selected FEEC family, so a
caller can inspect an intentional exact-sequence break without silently
mixing H(curl) and H(div) laws.  The affine restriction is deliberate:
curved-cell derivative maps and cut quadrature remain separate compositions.

`test_piola_enriched_vector` is an independent behavioral oracle: it computes
2D and 3D inverse/determinant formulas locally, checks both map kinds against
the primal result, checks central-difference JVPs and real dot-product VJPs,
and rejects singular and unknown map choices.  Exact-sequence preservation
across a cut, support activation, and higher-order curved-map construction
remain separate contracts.

`test_piola_enriched_differential_2d` and
`test_piola_enriched_differential_3d` independently evaluate the affine
Piola pullbacks and enrichment product terms, checks both FEEC families by
central differences, and verifies the complete real dot-product identity for
geometry and enrichment-gradient directions.
