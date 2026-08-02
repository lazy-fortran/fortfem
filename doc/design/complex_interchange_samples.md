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
