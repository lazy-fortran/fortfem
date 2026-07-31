# Linear level-set triangle interfaces

`evaluate_level_set_triangle_interface_2d` is the first internal-manifold
geometry primitive for unfitted interface work. It linearly interpolates three
nodal level-set values over a physical triangle and returns the two edge
intersections, their physical segment length, and the unit normal in the
direction of the physical level-set gradient.

The edge topology is intentionally explicit: a triangle with no proper cut, a
degenerate physical map, or a zero-gradient level set is rejected. A cut that
passes through a vertex is deduplicated, but derivative paths must still treat
that event as a connectivity change and rebuild the cut stencil. The focused
test uses the affine level set `x+y-0.4` and an independent segment/normal
oracle. Higher-dimensional level sets, cut-cell quadrature, and fixed-topology
shape derivatives remain Phase 2 work.
