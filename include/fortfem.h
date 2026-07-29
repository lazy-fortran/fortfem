#ifndef FORTFEM_H
#define FORTFEM_H

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Build the globally oriented edges of a zero-based triangular C mesh.
 *
 * vertices stores [x0, y0, x1, y1, ...], and triangles stores
 * [v0, v1, v2, ...]. Edges are returned in degree-of-freedom order as
 * [start0, end0, start1, end1, ...]. triangle_edge_dofs and
 * triangle_edge_signs contain three entries per triangle in local edge order
 * (v0,v1), (v1,v2), (v2,v0).
 *
 * status is 0 on success, -1 for invalid input, or -2 when edge_capacity is
 * too small. In the latter case, n_edges returns the required capacity.
 */
void fortfem_triangle_edge_map(
    int n_vertices,
    const double *vertices,
    int n_triangles,
    const int *triangles,
    int edge_capacity,
    int *n_edges,
    int *edges,
    int *triangle_edge_dofs,
    int *triangle_edge_signs,
    int *status);

#ifdef __cplusplus
}
#endif

#endif
