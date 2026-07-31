---
title: Quotient cell-complex contract
---

# Quotient cell-complex contract

`quotient_cell_complex` composes signed representative maps with the integer
boundary matrices of an oriented cell complex. It is the topology-only layer
between periodic or multipatch metadata and compatible finite-element
incidence matrices.

## Contract

For each chain group (C_k), a `cell_identification_t` supplies a canonical
representative and an orientation sign. The lower boundary of every
identified cell must equal the declared sign times the boundary of its
representative:

\[
  D_{k-1} S_k = S_{k-1} D_{k-1}^{\mathrm{quot}}.
\]

The routine applies this relation successively to (D_1), (D_2), and
(D_3), then validates the resulting quotient complex with exact integer
boundary-of-boundary products. A map that has inconsistent signs or does not
define a chain map is rejected; no averaging or tolerance is used.

## API

```fortran
call quotient_cell_complex( &
    complex, vertex_identification, edge_identification, quotient, status, &
    face_identification=face_identification, &
    volume_identification=volume_identification)
```

Face and volume identifications are required only when the source complex has
faces or volumes. Zero-cell chain groups are represented by zero-width or
zero-height allocatable matrices, so an interval can be closed into a circle
without special mesh code.

The returned complex owns copied integer incidence matrices. It does not own
coordinates, metric Hodge matrices, basis functions, trace ownership, or
application-specific periodic laws. Those belong to a higher geometry or
space layer.

## Independent oracle

The focused test closes the endpoints of an oriented interval and checks the
circle Betti numbers ((1,1,0,0)). It also supplies a malformed source
complex and verifies that quotienting rejects it before constructing output.
The existing signed-column tests independently exercise the chain-map sign
condition used by each boundary level.

