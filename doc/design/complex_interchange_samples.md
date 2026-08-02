---
title: Complex physical interchange samples
---

# Complex physical interchange samples

`complex_interchange_sample_set_t` is the physical-grid boundary for
frequency-domain and linear-response comparisons. It carries real physical
coordinates, positive quadrature weights, complex scalar/vector/tensor
components, producer, and provenance. Component ordering and normalization
remain adapter-owned; FortFEM does not parse GEQDSK, COCOS, GPEC, MARS-F,
GLISS, or STARWALL data.

`compare_complex_interchange_samples` requires matching coordinates and
weights within an explicit tolerance and reports weighted complex L2 absolute
and relative errors together with the componentwise maximum error. The
independent real-valued sample contract remains available for equilibrium
comparisons; this type is its complex counterpart for FEM/BEM/DtN/PML and
external ideal/resistive response fields.

The weighted L2 and relative metrics expose fixed-coordinate derivative
actions as well:

```fortran
call compare_complex_interchange_samples_jvp( &
    reference, candidate, coordinate_tolerance, reference_values_dot, &
    candidate_values_dot, weights_dot, absolute_error_dot, relative_error_dot, status)
call compare_complex_interchange_samples_vjp( &
    reference, candidate, coordinate_tolerance, absolute_error_bar, &
    relative_error_bar, reference_values_bar, candidate_values_bar, weights_bar, status)
```

The derivative uses the real-part convention
`Re(sum(conjg(values_bar)*values_dot))`; zero error and zero reference norm
are rejected, while the nonsmooth maximum norm remains diagnostic-only.
`test_complex_interchange_sample_derivatives` checks central differences and
the complex real-part adjoint identity.
