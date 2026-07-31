---
title: Tensor volume work
---

# Tensor volume work

`assemble_tensor_volume_work` is the neutral volume residual shared by
anisotropic transport, CGL pressure, Maxwell stress, and linear elasticity.
For quadrature point `q` and test function `i`, it evaluates

\[
  r_i=\sum_q w_q\,\nabla v_i(q):P(q),
\]

with no assumption that `P` is symmetric or that it came from a particular
constitutive model. `test_gradient(q,i,a,b)` and `stress(q,a,b)` are caller-
owned physical arrays; positive quadrature weights and fixed topology are
validated by the primitive.

```fortran
call assemble_tensor_volume_work( &
    test_gradient, stress, weights, residual, status)
call assemble_tensor_volume_work_jvp( &
    test_gradient, stress, weights, test_gradient_dot, stress_dot, &
    weights_dot, residual_dot, status)
call assemble_tensor_volume_work_vjp( &
    test_gradient, stress, weights, residual_bar, test_gradient_bar, &
    stress_bar, weights_bar, status)
```

The JVP applies the product rule to test gradients, stress, and quadrature
weights. The real VJP is the transpose contraction, so geometry, tensor
coefficients, and test-space data can all participate in an adjoint or shape
derivative. A generated CGL tensor or an elasticity/Maxwell client can supply
`stress`; FortFEM does not select its law or add plasma-specific closures.

`test_tensor_volume_work` compares against a separately written contraction,
central differences, the real dot-product identity, and an incompatible
residual-shape rejection.
