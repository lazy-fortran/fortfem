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

The edge topology is intentionally explicit: a triangle with no proper cut, a
degenerate physical map, or a zero-gradient level set is rejected. A cut that
passes through a vertex is deduplicated, but derivative paths must still treat
that event as a connectivity change and rebuild the cut stencil. The focused
interface test uses the affine level set `x+y-0.4` and an independent
segment/normal oracle; the quadrature test additionally checks polygon
centroids and affine integration. Higher-dimensional level sets, higher-order
cut-cell quadrature, and internal-manifold graphs remain Phase 2 work.

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
rejection. Fixed-topology tetrahedral JVPs and internal-manifold graphs remain
the next 3D geometry layer.

`evaluate_level_set_tetra_cut_quadrature_3d` now supplies that volume layer for
linear tetrahedra. It clips each oriented parent face, closes the positive and
negative polyhedra with the interface polygon, and accumulates exact volumes
and first moments from oriented tetrahedral fans. The resulting degree-one
quadrature data is checked against analytic one-corner cuts, uncut limits,
quadrilateral cuts, and volume/first-moment conservation. Fixed-topology
tetrahedral JVPs and internal-manifold graphs remain separate follow-up
contracts.
