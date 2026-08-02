---
title: Boundary-region graph contract
---

# Boundary-region graph contract

`boundary_region_graph_t` is the neutral free-boundary geometry boundary.  It
delegates oriented plus/minus region incidence, connectivity, and cycle-basis
metadata to the existing fixed-topology region graph, then attaches
caller-owned interface metadata:

- nonnegative genus per interface;
- an exterior-interface flag;
- contiguous physical samples with `sample_offsets(i):sample_offsets(i+1)-1`
  belonging to interface `i`;
- three-dimensional points, oriented normals, and positive quadrature weights.

The wrapper validates finite arrays, positive weights, nonzero normals, offset
ownership, and topology without deciding whether the samples came from
triangles, high-order panels, Fourier surfaces, NURBS, or IGA patches.  The
`boundary_region_graph_interface_samples` accessor returns an independent copy
of one interface's samples, so a BEM, DtN, PML, virtual-casing, or external
free-boundary client can consume the same physical sampler contract.
`boundary_region_graph_interface_metadata` returns the genus and exterior flag
without exposing the graph's private storage.

Geometry values and shape derivatives remain caller-owned operator data.  This
type is metadata, not a mesh or an equilibrium representation; it has no coil,
profile, GEQDSK/COCOS, VMEC, SPEC, GPEC, or MARS logic.
