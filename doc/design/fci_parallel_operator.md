# FCI parallel support operator

FortFEM now has a small, dependency-light algebraic contract for a
field-coordinate-independent (FCI) parallel derivative. It is the part of the
PARALLAX design that belongs in a reusable FEM/operator library: interpolation
maps are supplied by a geometry or field-line-tracing service, while FortFEM
assembles the sparse staggered gradient and its conservative support
divergence.

The geometry side now also exposes `trace_fci_field_line_rk4`, a fixed-step
classical RK4 service with a callback for `d(point)/dphi`, and
`trace_fci_field_line_rk4_jvp`, which advances a supplied tangent through the
same stages. Both are intentionally agnostic about magnetic-field storage and
mesh lookup. A constant-velocity oracle is exact, and exponential tests check
the primal and tangent endpoints before the resulting points are passed to a
separate interpolation/map builder.
The first such builder is `build_fci_linear_interpolation_map_1d`, which
provides partition-of-unity and affine-reproduction contracts for coordinate
slices without coupling the operator to a mesh lookup implementation. Its
fixed-topology JVP and VJP are also available for differentiable geometry
updates away from stencil-cell crossings. A Cartesian bilinear map builder is
available for a first 2D poloidal slice, with x-fast source-column ordering.
Its source-grid and target-point JVP/VJP products are available on a fixed
cell topology as well.
`compute_fci_staggered_flux_box_volumes` combines traced forward/backward flux
expansion with plane-cell area and (B_\varphi) for the support weights.

For segment `k`, let `Q_plus(k)` map the upper poloidal plane to staggered flux
boxes and `Q_minus(k)` map the lower plane. With line lengths `ell`, the
gradient is

\[
  Q_k = \operatorname{diag}(\ell_k^{-1})
        \left(Q_{+}(k)-Q_{-}(k)\right).
\]

For canonical-plane volumes `W_c` and staggered volumes `W_s`, the support
divergence is assembled as

\[
  P = -W_c^{-1}Q^T W_s.
\]

Therefore the discrete power/adjoint identity is exact up to floating-point
roundoff,

\[
  u^T W_c P f = -(Qu)^T W_s f.
\]

With positive cellwise parallel coefficients `K_parallel`, the matrix-free
diffusion action is

\[
  L_\parallel u = -W_c^{-1}Q^T W_s K_\parallel Q u,
  \qquad
  u^T W_cL_\parallel u = -(Qu)^T W_sK_\parallel Qu.
\]

The public entry points are:

```fortran
call assemble_fci_parallel_gradient_csc( &
    forward_map, backward_map, line_lengths, gradient, status)
call assemble_fci_parallel_support_divergence_csc( &
    gradient, canonical_volumes, staggered_volumes, divergence, status)
call apply_fci_parallel_diffusion( &
    forward_map, backward_map, line_lengths, parallel_coefficient, &
    canonical_volumes, staggered_volumes, field, diffusion_field, status)
call apply_fci_parallel_diffusion_jvp( &
    forward_map, backward_map, line_lengths, parallel_coefficient, &
    canonical_volumes, staggered_volumes, field, forward_map_dot, &
    backward_map_dot, line_lengths_dot, parallel_coefficient_dot, &
    canonical_volumes_dot, staggered_volumes_dot, field_dot, &
    diffusion_field_dot, status)
call apply_fci_parallel_diffusion_vjp( &
    forward_map, backward_map, line_lengths, parallel_coefficient, &
    canonical_volumes, staggered_volumes, field, diffusion_field_bar, &
    forward_map_bar, backward_map_bar, line_lengths_bar, &
    parallel_coefficient_bar, canonical_volumes_bar, staggered_volumes_bar, &
    field_bar, status)
call compute_fci_parallel_diffusion_diagonal( &
    forward_map, backward_map, line_lengths, parallel_coefficient, &
    canonical_volumes, staggered_volumes, diagonal, status)
call apply_fci_parallel_jacobi_preconditioner( &
    diagonal, residual, correction, status)
call apply_fci_anisotropic_diffusion( &
    perpendicular_operators, forward_map, backward_map, line_lengths, &
    parallel_coefficient, canonical_volumes, staggered_volumes, field, &
    diffusion_field, status)
call apply_fci_parallel_diffusion_field_vjp( &
    forward_map, backward_map, line_lengths, parallel_coefficient, &
    canonical_volumes, staggered_volumes, diffusion_field_bar, field_bar, &
    status)
call apply_fci_parallel_gradient_jvp( &
    forward_map, backward_map, line_lengths, field, forward_map_dot, &
    backward_map_dot, line_lengths_dot, field_dot, gradient_dot, status)
call apply_fci_parallel_gradient_vjp( &
    forward_map, backward_map, line_lengths, field, gradient_bar, &
    forward_bar, backward_bar, line_bar, field_bar, status)
call build_fci_bilinear_interpolation_maps_2d( &
    source_x, source_y, forward_x, forward_y, backward_x, backward_y, &
    forward_map, backward_map, status)
call build_fci_bilinear_interpolation_maps_2d_jvp( &
    source_x, source_y, forward_x, forward_y, backward_x, backward_y, &
    source_x_dot, source_y_dot, forward_x_dot, forward_y_dot, &
    backward_x_dot, backward_y_dot, forward_map_dot, backward_map_dot, status)
call build_fci_bilinear_interpolation_maps_2d_vjp( &
    source_x, source_y, forward_x, forward_y, backward_x, backward_y, &
    forward_map_bar, backward_map_bar, source_x_bar, source_y_bar, &
    forward_x_bar, forward_y_bar, backward_x_bar, backward_y_bar, status)
```

The canonical unknowns are ordered plane-by-plane. Staggered rows are ordered
segment-by-segment. `forward_map` and `backward_map` have shape
`(n_staggered, n_plane, n_segment)`, and `line_lengths` has shape
`(n_staggered, n_segment)`. Zero interpolation coefficients are omitted from
the CSC matrix, so the resulting action is matrix-free at the solver level
through the existing FortSparse `csc_matvec` interface.

The batched bilinear adapter turns traced endpoint coordinates with shape
`(n_staggered, n_segment)` into the forward and backward map tensors consumed
by the support operator. It applies the single-slice builder independently to
each toroidal segment, preserving the x-fast source-column ordering and
rejecting any segment that leaves the source box. Its JVP/VJP paths reuse the
fixed-cell single-slice contracts and accumulate source-grid cotangents over
all segments. The topology rule is deliberate: a traced endpoint on a grid
line is valid for the primal map but is rejected by the derivative paths until
the stencil cell is rebuilt.

The focused test uses identity maps and an analytically linear field as an
independent oracle. It also checks the weighted adjoint identity and a flux
balance vector. A separate diffusion test checks the explicit `P K_parallel Q`
oracle, the weighted negative-energy identity, and the diffusion field VJP
against an independent dot-product identity. Field-line tracing,
interpolation stencils, support-volume construction, perpendicular terms,
boundary conditions, and multigrid remain separate follow-up ingredients. The
pointwise gradient and local diffusion contribution value, JVP, and VJP
products are generated by the pinned FortSym revision. The full diffusion AD
test checks maps, line lengths, coefficients, both volume weights, and the
field against a central-difference JVP and an independent real dot-product
VJP identity. These are fixed-topology contracts: a field-line stencil or
interpolation-cell crossing is a discrete event and must be handled by a
separate remap/rebuild step. This API deliberately does not copy PARALLAX
implementation code.

`compute_fci_parallel_diffusion_diagonal` returns the positive diagonal of
`-L_parallel`, namely `diag(W_c^{-1} Q^T W_s K_parallel Q)`. It is a cheap,
geometry-aware Jacobi baseline for the field solver and a stable contract for
later plane multigrid or field-split preconditioners. Its per-stencil
`Q^2` contribution is generated by FortSym and is checked independently in
`test_fci_parallel_diagonal`.
`apply_fci_parallel_jacobi_preconditioner` performs the corresponding
positive diagonal solve and rejects zero or negative entries before division.
`apply_fci_anisotropic_diffusion` adds one square CSC perpendicular block per
canonical plane to the parallel support action. Each plane block is validated
against the FCI plane size, so the split remains explicit and compatible with
PARALLAX-style independent 2D elliptic solves plus a matrix-free 3D line
operator. The anisotropic test uses negative symmetric plane blocks and checks
the combined dissipative energy oracle.

## Provenance

The design follows the FCI trace/interpolation and support-operator
descriptions in the local PARALLAX checkout at
`/home/ert/code/genex/parallax`, especially its `pages/description/fci`
material. The published support-operator construction is described by
[Stegmeir et al.](https://doi.org/10.1016/j.cpc.2016.12.014). No PARALLAX
source, binary, or license-sensitive benchmark data is included in FortFEM.
