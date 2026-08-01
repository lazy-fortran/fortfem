---
title: Curvilinear PML geometry builder
---

# Curvilinear PML geometry builder

`fortfem_curvilinear_pml_geometry` supplies a small, solver-neutral geometry
layer for curved PML cells.  Each cell receives a point on the local layer and
a unit normal from the caller.  For the cell centroid (c), layer point (o),
and normal (n), the signed distance is

\[
d = n\mathbin{\cdot}(c-o).
\]

Active cells ((d>0)) receive the full tensor

\[
S = I + i\,a\,n n^T,
\qquad
a = \frac{\sigma_{\max}}{k}\left(\frac{d}{w}\right)^p,
\]

while cells on the physical side retain (S=I).  Since (n n^T) is formed in
physical coordinates, a curved layer produces off-diagonal entries and can be
passed directly to the full curvilinear Helmholtz or curl--curl coefficient
contract.  No closest-point search, surface representation, or mesh policy is
imposed; those remain owned by the FEM, IGA, or boundary client.

The public routines are:

```fortran
call build_curvilinear_normal_pml_element_stretch( &
    vertices, cells, layer_origins, layer_normals, layer_width, &
    wave_number, sigma_max, polynomial_degree, stretch, status)
call build_curvilinear_normal_pml_element_stretch_jvp( &
    vertices, cells, layer_origins, layer_normals, layer_width, &
    wave_number, sigma_max, polynomial_degree, vertices_dot, origins_dot, &
    normals_dot, width_dot, wave_number_dot, sigma_max_dot, stretch_dot, &
    status)
call build_curvilinear_normal_pml_element_stretch_vjp( &
    vertices, cells, layer_origins, layer_normals, layer_width, &
    wave_number, sigma_max, polynomial_degree, stretch_bar, vertices_bar, &
    origins_bar, normals_bar, width_bar, wave_number_bar, sigma_max_bar, &
    status)
```

The attenuation and its scalar derivatives use the FortSym-generated Cartesian
PML products.  The JVP and VJP are valid for a fixed active-cell set; a layer
client should keep cell centroids away from (d=0) while differentiating.
Normals are required to be finite and unit length.  The VJP uses the real-part
complex pairing and returns real geometry and layer-parameter cotangents.  The
independent `test_curvilinear_pml_geometry_ad` test checks the full tensor
oracle, central differences, the complex dot identity, and invalid-normal
rejection.
