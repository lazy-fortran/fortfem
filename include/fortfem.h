#ifndef FORTFEM_H
#define FORTFEM_H

#include <complex.h>

#ifdef __cplusplus
typedef __complex__ double fortfem_complex;
extern "C" {
#else
typedef double complex fortfem_complex;
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

/*
 * Persistent variant of fortfem_triangle_edge_map. Mesh handles are positive,
 * process-local integers. Call fortfem_triangle_mesh_free when finished.
 */
void fortfem_triangle_mesh_create(
    int n_vertices,
    const double *vertices,
    int n_triangles,
    const int *triangles,
    int *handle,
    int *n_edges,
    int *status);

void fortfem_triangle_mesh_edges(
    int handle,
    int edge_capacity,
    int *n_edges,
    int *edges,
    int *triangle_edge_dofs,
    int *triangle_edge_signs,
    int *status);

void fortfem_triangle_mesh_free(int handle, int *status);

/*
 * RT0 coefficients use the global degree-of-freedom order returned by
 * fortfem_triangle_mesh_edges. Triangles are zero-based at this boundary.
 */
void fortfem_rt0_evaluate(
    int handle,
    int triangle,
    double x,
    double y,
    int n_dofs,
    const fortfem_complex *dofs,
    fortfem_complex *value,
    fortfem_complex *divergence,
    int *status);

void fortfem_rt0_l2_norm(
    int handle,
    int n_dofs,
    const fortfem_complex *dofs,
    double *norm,
    int *status);

void fortfem_rt0_toroidal(
    int handle,
    int mode,
    int n_dofs,
    const fortfem_complex *dofs,
    fortfem_complex *toroidal,
    int *status);

/*
 * Assemble the real axisymmetric Nedelec operator
 *
 *   integral R curl(w) curl(A) + mode^2 / R w dot A
 *
 * on a retained counter-clockwise triangle mesh with positive radial
 * coordinates. The returned CSC arrays are zero-based. n_dofs and nnz are
 * always set for valid input; status is -2 when nnz_capacity is too small.
 */
void fortfem_nedelec_axisymmetric_fourier_csc(
    int handle,
    int mode,
    int quadrature_degree,
    int nnz_capacity,
    int *n_dofs,
    int *nnz,
    int *col_ptr,
    int *row_ind,
    double *values,
    int *status);

/*
 * Assemble the mixed integral of a Nedelec test field with an RT0 trial
 * field. Storage and capacity conventions match the Fourier operator above.
 */
void fortfem_nedelec_rt0_mass_csc(
    int handle,
    int quadrature_degree,
    int nnz_capacity,
    int *n_dofs,
    int *nnz,
    int *col_ptr,
    int *row_ind,
    double *values,
    int *status);

/*
 * Factor a zero-based complex CSC matrix and return a positive opaque handle.
 * Handles retain the factorization across solves and must be released with
 * fortfem_factor_free. The handle registry is process-local.
 */
void fortfem_complex_factor_csc(
    int n,
    int nnz,
    const int *col_ptr,
    const int *row_ind,
    const fortfem_complex *values,
    int *handle,
    int *status);

void fortfem_complex_solve(
    int handle,
    int n,
    const fortfem_complex *rhs,
    fortfem_complex *solution,
    int *status);

void fortfem_factor_free(int handle, int *status);

#ifdef __cplusplus
}
#endif

#endif
