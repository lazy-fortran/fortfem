---
title: Weighted wave reflection diagnostics
---

# Weighted wave reflection diagnostics

`fortfem_wave_reflection_diagnostics` is a solver- and geometry-neutral
diagnostic layer for FEM, BEM, DtN, PML, and external-code samples. Given a
reference field (u_r), a candidate field (u_c), and positive physical
sample weights (w_s), it reports

\[
 e = \left(\sum_s w_s\lvert u_{c,s}-u_{r,s}\rvert^2\right)^{1/2},
 \qquad
 e_{\rm rel}=e/\left(\sum_s w_s\lvert u_{r,s}\rvert^2\right)^{1/2}.
\]

The component dimension is included in the sum. For an incident field (u_i)
and total field (u_t), the named reflection coefficient is

\[
 R=\|u_t-u_i\|_w/\|u_i\|_w.
\]

The public routines provide primal, JVP, and real-part complex VJP actions for
the absolute and relative error, plus the reflection-coefficient wrapper:

```fortran
call evaluate_weighted_complex_error(..., absolute_error, relative_error, status)
call evaluate_weighted_complex_error_jvp(..., absolute_error_dot, &
    relative_error_dot, status)
call evaluate_weighted_complex_error_vjp(..., absolute_error_bar, &
    relative_error_bar, reference_bar, candidate_bar, weights_bar, status)

call evaluate_weighted_reflection_coefficient( &
    incident, total, weights, coefficient, status)
```

The reverse convention is

\[
 \operatorname{Re}(\bar e^H\dot e)
 =\operatorname{Re}(\bar u_r^H\dot u_r)
 +\operatorname{Re}(\bar u_c^H\dot u_c)
 +\bar w^T\dot w.
\]

Samples, weights, and all active directions must be finite and shape-compatible
with positive weights. Derivatives reject a zero error or zero reference norm,
where the norm ratio is not classically differentiable. Coordinates and mesh
motion are intentionally outside this layer: a caller composes its physical
sample and quadrature geometry JVP/VJP before invoking the metric.

`test_wave_reflection_diagnostics_ad` uses an incident plane wave plus a known
counter-propagating reflected wave. It checks the exact coefficient, both JVPs
against central differences, the VJP dot identity, and invalid-weight rejection.
