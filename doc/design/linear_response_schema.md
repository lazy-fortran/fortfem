---
title: Linear-response interchange schema
---

# Linear-response interchange schema

`write_linear_response_interchange` and
`read_linear_response_interchange` provide a small, versioned, neutral text
boundary for reusable response matrices. The record contains only modal
metadata, provenance, complex equilibrium/inertia/resistive/vacuum/wall
blocks, and named response channels. It does not parse GPEC, MARS-F, GLISS,
STARWALL, or any other application format.

The first line is the magic string
`FORTFEM_LINEAR_RESPONSE_TEXT 1`; the following metadata and dimension lines
are followed by labels and column-major complex matrix entries. Complex
values are stored as real/imaginary pairs, so the file remains inspectable
with ordinary text tools. A reader validates the schema version, dimensions,
finite values, labels, and block shapes before publishing a record.

Reads apply a bounded payload guard before allocation. This is deliberate:
an interchange file is an adapter boundary and must not be able to request an
unbounded response matrix. Larger or compressed data belongs in a separate
caller-owned backend; the neutral record remains a small dense oracle for
reciprocity, passivity, Schur, FEM/BEM/DtN, and wall-response tests.

Both routines return `status == 0` on success and a nonzero status for invalid
records, open failures, malformed data, or failed validation. The schema is
metadata-only and therefore remains usable by equilibrium, linear-response,
free-boundary, and structure-preserving time clients without importing their
physics or normalization conventions.
