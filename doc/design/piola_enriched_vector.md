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
```

The JVP differentiates the inverse/determinant Piola factors, reference
coefficients, and activation.  The VJP uses FortNum’s FortSym-derived inverse
and determinant reverse kernels and returns the geometry, reference, and
activation cotangents.  This makes shape sensitivity and enrichment
conditioning composable without duplicating covariant/contravariant formulas
in IGA or XFEM clients.

`test_piola_enriched_vector` is an independent behavioral oracle: it computes
2D and 3D inverse/determinant formulas locally, checks both map kinds against
the primal result, checks central-difference JVPs and real dot-product VJPs,
and rejects singular and unknown map choices.  Exact-sequence preservation
across a cut, support activation, and higher-order curved-map construction
remain separate contracts.
