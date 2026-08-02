---
title: Field-aligned constitutive tensor
---

# Field-aligned symmetric and gyrotropic tensor

`evaluate_field_aligned_constitutive_tensor` provides the neutral pointwise
constitutive tensor

\[
  K = k_\perp I + (k_\parallel-k_\perp)bb^T + k_H[b]_\times,
  \qquad [b]_\times v=b\times v,
\]

for a caller-owned unit direction `b`.  The trailing `hall_coefficient`
argument is optional and defaults to zero.  Thus the same API covers a
symmetric field-aligned tensor and a nonsymmetric Hall/gyrotropic extension;
it does not select a plasma, transport, or material closure.

```fortran
call evaluate_field_aligned_constitutive_tensor( &
    k_parallel, k_perpendicular, b, tensor, status, hall_coefficient)
call evaluate_field_aligned_constitutive_tensor_jvp( &
    k_parallel, k_perpendicular, b, k_parallel_dot, k_perpendicular_dot, &
    b_dot, tensor_dot, status, hall_coefficient, hall_dot)
call evaluate_field_aligned_constitutive_tensor_vjp( &
    k_parallel, k_perpendicular, b, tensor_bar, k_parallel_bar, &
    k_perpendicular_bar, b_bar, status, hall_coefficient, hall_bar)
```

The symmetric projector part composes with the generated CGL tensor.  The
cross-product matrix is assembled explicitly so that its Hall coefficient and
direction derivatives are visible in the JVP/VJP.  The wrapper rejects
non-finite coefficients, cotangents, or increments, non-three-dimensional
arrays, and directions whose norm differs from one by more than the fixed
unit-direction tolerance.  Direction increments are ambient derivatives; a
caller imposing a unit-vector chart should supply its tangent projection.

The focused test has an independent matrix oracle, central-difference JVP,
real transpose identity including the Hall coefficient, optional-Hall limit,
and invalid-input cases.  Geometry, quadrature, and constitutive laws remain
outside this pointwise contract.
