---
title: Tensor power and dissipation split
---

# Tensor power and dissipation split

`evaluate_tensor_power_split` is a neutral pointwise contract for a real
constitutive tensor (K) and vector (x). It reports

\[
p_s=x^T\operatorname{Sym}(K)x,\qquad
p_a=x^T\operatorname{Skew}(K)x,\qquad
p=x^TKx.
\]

The antisymmetric contribution is exactly zero for a real vector. Thus a
caller can use `p_s` as the dissipative/energy contribution and `p_a` as an
independent Hall or gyrotropic power diagnostic without FortFEM selecting a
plasma closure. The tensor may be a field-aligned pullback, an anisotropic
elastic tangent, a wave impedance, or any other caller-owned constitutive
block.

The accompanying JVP and VJP cover the tensor and vector arguments. The
contract is pointwise: quadrature weights, geometry maps, and integration
over elements or field lines remain caller-owned. In particular, no
regularization is applied to a nearly singular anisotropic coefficient.

`test_tensor_power_split` uses an independently assembled symmetric/skew
oracle, central differences for all three reported powers, and the real
transpose identity for the complete VJP. It also checks that a non-finite
tensor is rejected.
