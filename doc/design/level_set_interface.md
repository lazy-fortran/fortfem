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
cut-cell quadrature, and fixed-topology shape derivatives remain Phase 2 work.
