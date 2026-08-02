---
title: Geometry-generated tetrahedral Nedelec PML
---

# Geometry-generated tetrahedral Nédélec PML

The `assemble_tetra_nedelec_geometry_pml_csc` family composes the neutral
curvilinear PML geometry contract with the existing order-one-or-higher
tetrahedral H(curl) CSC assembly:

1. each cell receives a caller-owned layer origin, unit normal, width, and
   attenuation parameters;
2. the fixed-active-set geometry builder creates the full complex stretch;
3. the existing covariant-Piola curl--curl/mass assembly consumes that tensor.

The JVP and VJP propagate through both layers, including mesh vertices,
origins, normals, widths, wave number, and attenuation. The reverse path
accumulates mesh cotangents from the finite-element geometry and from the PML
distance map. No closest-point search, active-cell update, or global matrix
allocation beyond the caller-owned CSC assembly is hidden here. A cell
crossing its layer plane is a fixed-topology event and must be handled by the
caller.

The independent test compares the full CSC values with perturbed
geometry-generated reassemblies and checks the complete real-part complex
adjoint identity. This is the geometry chain needed before comparing vector
FEM/PML fields with BEM or DtN on curved external surfaces.
