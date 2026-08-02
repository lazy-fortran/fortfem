---
title: Linear level-set triangle interfaces
---

# Linear level-set triangle interfaces

`evaluate_level_set_triangle_interface_2d` is the first internal-manifold
geometry primitive for unfitted interface work. It linearly interpolates three
nodal level-set values over a physical triangle and returns the two edge
intersections, their physical segment length, and the unit normal in the
direction of the physical level-set gradient.

`evaluate_level_set_triangle_cut_areas_2d` uses the same linear cut to clip the
positive and negative sub-polygons exactly and returns both physical areas plus
the interface length. Its focused test checks the affine `x+y-0.4` case against
an independent clipped-triangle oracle and verifies parent-area conservation.

`evaluate_level_set_triangle_cut_quadrature_2d` extends that primitive with a
centroid for each non-empty side and the oriented interface normal. The
centroid/area pair is an exact degree-one quadrature rule, so constants and
affine manufactured fields integrate without an additional approximation. The
zero-area side is returned with a zero centroid, and an identically zero level
set is rejected because it has no defined positive/negative partition.

`evaluate_level_set_triangle_cut_moments_2d` adds the symmetric raw tensor

\[
M_{ab}^{\pm}=\int_{\Omega^{\pm}}x_a x_b\,dA,
\]

using the same clipped polygon and orientation. This is an exact degree-two
moment contract: the positive and negative tensors sum to the parent triangle
moment, including for a cut cell. Its fixed-topology JVP differentiates the
edge intersections and polygon moments directly, so quadratic manufactured
loads can use the cut geometry without a finite-difference kernel.

`evaluate_level_set_triangle_cut_fourth_moments_2d` extends the same contract
to the symmetric rank-four tensor

\[
M_{abcd}^{\pm}=\int_{\Omega^{\pm}}x_a x_b x_c x_d\,dA.
\]

The implementation uses the exact Green-theorem primitive
\(\int_K x^p y^q\,dA=(p+1)^{-1}\oint_{\partial K}x^{p+1}y^q\,dy\)
and a finite binomial edge sum, so it remains exact for every clipped
triangle or quadrilateral. Its fixed-topology JVP differentiates the same
primitive term by term. The focused test checks an independent simplex
monomial oracle, conservation of the parent quartic tensor, and a central-
difference JVP oracle. This is the first degree-four cut rule; the matching
tetrahedral degree-four rule and general curved-cell moment fitting remain
open.

The edge topology is intentionally explicit: a triangle with no proper cut, a
degenerate physical map, or a zero-gradient level set is rejected. A cut that
passes through a vertex is deduplicated, but derivative paths must still treat
that event as a connectivity change and rebuild the cut stencil. The focused
interface test uses the affine level set `x+y-0.4` and an independent
 segment/normal oracle; the quadrature test additionally checks polygon
centroids and affine integration. The cut-moment test checks an independent
one-corner quadratic oracle, parent-moment conservation, and a central-
difference JVP oracle. Higher-dimensional higher-order cut-cell quadrature
remains Phase 2 work; the neutral
internal-manifold graph is documented separately.

`evaluate_level_set_triangle_interface_2d_jvp` provides the fixed-topology
forward derivative of the same primitive. It differentiates both the physical
triangle vertices and the nodal level-set values, returning derivatives of the
two ordered intersection points, segment length, and normalized level-set
normal. The implementation differentiates the edge fractions and the affine
physical gradient directly, so it does not finite-difference geometry inside a
kernel. Nodal level-set zeros are rejected as topology events; callers must
rebuild the cut stencil after an interface crosses a vertex or an edge-crossing
pair changes. Its focused test compares an independently derived intersection
and normal formula with a central-difference oracle while checking that a
topology event is rejected.

`evaluate_level_set_triangle_cut_quadrature_2d_jvp` extends the same contract
to the exact degree-one positive/negative subcell areas and centroids. It
propagates the edge-intersection derivatives through the clipped-polygon
shoelace moments, and composes the interface JVP for length and normal. It
supports fixed uncut states as well as proper cuts, rejects nodal zeros, and is
checked against central differences plus independent total-area and
first-moment conservation identities.

`evaluate_level_set_tetra_interface_3d` is the first higher-dimensional cut
primitive. It computes the affine physical level-set gradient, collects the
three or four edge intersections, orders the convex polygon in its physical
plane, and returns its area and gradient-oriented normal. Its test covers both
triangular and quadrilateral cuts against independent edge points and plane
geometry, as well as uncut, nodal-topology, and degenerate-tetrahedron
rejection. Fixed-topology tetrahedral JVPs and the neutral internal-manifold
graph now form the next 3D geometry attachment layer; geometry-to-graph
construction remains a client composition concern.

`evaluate_level_set_tetra_cut_quadrature_3d` now supplies that volume layer for
linear tetrahedra. It clips each oriented parent face, closes the positive and
negative polyhedra with the interface polygon, and accumulates exact volumes
and first moments from oriented tetrahedral fans. The resulting degree-one
quadrature data is checked against analytic one-corner cuts, uncut limits,
quadrilateral cuts, and volume/first-moment conservation. Fixed-topology
tetrahedral JVPs now propagate moving vertices, level values, clipped-face
intersections, volumes, centroids, interface area, and gradient-oriented normal.
The JVP is checked against central differences, side-volume/first-moment
conservation, and the independent interface-JVP contract. The neutral
internal-manifold graph supplies the topology attachment; geometry-to-graph
construction remains a client composition layer.

`evaluate_level_set_tetra_cut_fourth_moments_3d` extends the same oriented
tetrahedral fan to the symmetric rank-four raw tensor. Each fan tetrahedron is
integrated by the exact barycentric multinomial sum, and its fixed-topology JVP
differentiates the vertex products directly. The focused test checks an
independent simplex-moment oracle, conservation of the parent quartic tensor,
and a central-difference JVP oracle. Curved-cell moment fitting and adaptive
cut integration remain separate higher-level layers.

`evaluate_level_set_tetra_cut_moments_3d` extends the same oriented fan with
the symmetric raw second-moment tensor

\[
M_{ab}^{\pm}=\int_{\Omega^{\pm}}x_a x_b\,dV.
\]

The tetrahedral moment identity is exact for quadratic manufactured fields and
the positive/negative tensors conserve the parent tensor. Its fixed-topology
JVP differentiates both clipped-face and interface-fan contributions. The
focused test covers an independent one-corner tetrahedron, conservation, and a
central-difference derivative oracle.
