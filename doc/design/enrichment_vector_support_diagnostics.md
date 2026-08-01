---
title: Vector enrichment-support conditioning
---

# Vector enrichment-support conditioning

`fortfem_enrichment_support_vector_diagnostics` provides the metric-aware
support Gram contraction needed to detect rank loss and conditioning problems
in vector XFEM/XIGA enrichments.  The caller supplies physical vector values
at quadrature samples.  These may be covariant Nédélec values, contravariant
Raviart--Thomas/BDM values, mapped IGA fields, or another compatible vector
representation; the upstream Piola or spline map is therefore not duplicated
inside the conditioning primitive.

For vector basis values `b_i(x_p)`, symmetric positive-definite metric
`M(x_p)`, and nonnegative quadrature weights `w_p`, the active-support Gram
matrix is

```
G_ij = sum_p w_p b_i(x_p)^T M(x_p) b_j(x_p).
```

Only the caller's fixed activation mask contributes.  The metric is required
to be symmetric positive definite in two or three physical components, so a
negative or nonsymmetric material/geometry metric is rejected rather than
silently producing a misleading rank estimate.  The resulting `gram` can be
passed directly to `evaluate_enrichment_support_rank_condition`.

```fortran
call evaluate_enrichment_support_vector_gram( &
    basis_values, metric, sample_weights, active_mask, gram, status)
call evaluate_enrichment_support_vector_gram_jvp( &
    basis_values, metric, sample_weights, active_mask, basis_dot, metric_dot, &
    weights_dot, gram_dot, status)
call evaluate_enrichment_support_vector_gram_vjp( &
    basis_values, metric, sample_weights, active_mask, gram_bar, basis_bar, &
    metric_bar, weights_bar, status)
```

The JVP differentiates physical vector values, metric coefficients, and
quadrature weights while retaining the discrete activation space.  This lets
callers propagate Piola Jacobian or IGA control-point derivatives into the
conditioning diagnostic without hiding geometry in a second implementation.
The VJP returns all three cotangents and satisfies the real dot-product
identity, including nonsymmetric Gram cotangents while metric directions are
restricted to the declared symmetric metric manifold.

The independent test constructs physical values from a hand-written affine
covariant-Piola formula, evaluates the contraction with an explicit metric
oracle, checks central-difference and real-adjoint identities, and rejects a
nonsymmetric metric.  Cut-cell activation, map construction, and the final
rank policy remain separate contracts and can therefore be reused by IGA,
fitted FEEC, DG, or HDG callers.
