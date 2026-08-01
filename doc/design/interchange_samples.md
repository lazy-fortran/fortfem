---
title: Physical interchange sample sets
---

# Physical interchange sample sets

`interchange_sample_set_t` is the small, solver-independent boundary for
comparing fields produced by different codes. An adapter supplies one common
set of physical coordinates, one value vector at every point, positive
quadrature weights, and provenance. The contract does not interpret an
equilibrium file, a coordinate convention, a PDE, or a plasma closure.

The values may be scalar, vector, or tensor components. Component ordering is
owned by the adapter and recorded in its provenance; FortFEM only checks that
the candidate and reference have the same declared shape and physical points.
`compare_interchange_samples` then reports a weighted absolute (L^2) error,
relative error against the reference norm, and componentwise (L^\infty)
error. A coordinate and weight tolerance is required, so comparing two
different meshes silently is impossible.

This is intended for license-safe samplers for CHEASE, FreeGS, GPEC, MARS,
GLISS, VMEC, SPEC, or any other external producer. The external code and its
reader remain outside FortFEM; only the sampled numerical arrays and a
revision/provenance record cross the boundary.

The independent test covers initialization and validation, deep-copy
assignment, weighted vector errors, and rejection of mismatched physical
coordinates. It is deliberately separate from any external solver run.
