---
title: Geometry-generated scalar tetrahedral PML
---

# Geometry-generated scalar tetrahedral PML

The `assemble_tetra_lagrange_geometry_pml_csc` family is the scalar
Helmholtz counterpart of the geometry-generated Nédélec PML chain. It
composes caller-owned per-cell layer origins, unit normals, widths, wave
number, and attenuation with the existing curvilinear scalar CSC assembly.

Value, JVP, and VJP products propagate through the normal-frame stretch and
the P1 or arbitrary-degree finite-element geometry. The reverse path returns
mesh, layer, width, wave-number, and attenuation cotangents. Active-cell
classification and closest-point construction are intentionally outside the
API; crossing the layer plane is a fixed-topology event.

The independent test checks CSC pattern preservation, geometry-generated
central reassembly, and the complete real-part complex adjoint identity.
