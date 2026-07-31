---
title: Tree-cotree gauge contract
---

# Tree-cotree gauge contract

`fortfem_tree_cotree_gauge` provides the topology-only reduction needed to
make a direct curl--curl system nonsingular. Given an oriented vertex--edge
incidence matrix, it constructs a spanning forest. Tree-edge degrees of
freedom are fixed to zero; cotree degrees of freedom are retained.

For a fixed restriction (R), a dense curl--curl block is reduced as

\[
    K_c=R^T K R,\qquad b_c=R^T b.
\]

The tree selector is a discrete topology operation. Its choice is frozen for
coefficient, geometry, and state differentiation. JVP/VJP actions of the
reduced operator therefore select the corresponding entries of the full
operator; a graph rebuild or tree change is reported as a topology event and
is not silently differentiated.

## API

```fortran
call build_tree_cotree_gauge(incidence, gauge, status)
call tree_cotree_gauge_edges(gauge, tree_edges, cotree_edges, status)
call apply_tree_cotree_restriction(gauge, full, reduced, status)
call apply_tree_cotree_prolongation(gauge, reduced, full, status)
call reduce_tree_cotree_dense_system( &
    gauge, matrix, rhs, reduced_matrix, reduced_rhs, status)
call reduce_tree_cotree_dense_system_jvp( &
    gauge, matrix_dot, rhs_dot, reduced_matrix_dot, reduced_rhs_dot, status)
call reduce_tree_cotree_dense_system_vjp( &
    gauge, reduced_matrix_bar, reduced_rhs_bar, matrix_bar, rhs_bar, status)
```

The graph must have exactly one `+1` and one `-1` per edge column. Parallel
edges are allowed; self-loops and ambiguous columns are rejected because a
control graph must supply their topology explicitly. Disconnected graphs
produce a spanning forest. On a multiply connected domain, non-tree edges
include cotree loop representatives; period constraints and harmonic cuts are
separate contracts.

## FEEC and IGA scope

For lowest-order edge elements the graph is the mesh one-skeleton. For
high-order Nédélec and spline H(curl) spaces, only the gradient/control-mesh
subspace is selected by the tree; higher-order edge, face, and cell moments
are retained through a caller-owned DOF map. The same rule applies to the IGA
control mesh and to nonmatching multipatch mortars. The primitive never
assumes a plasma geometry or a particular material law.

Tree--cotree gauging is the direct-solve path for curl--curl nullspaces. It is
complementary to compatible nullspace projection and to ILU/IC preconditioners
for iterative solves; no preconditioner is hidden inside this topology module.

## Independent oracle and literature

The triangle test checks the exact spanning-tree/cotree split, restriction and
zero-tree prolongation, and the reduced direct block. A disconnected graph
checks forest behavior and an invalid incidence column is rejected.

The implementation follows the graph-gauge motivation in [Manges and
Cendes (1995)](https://doi.org/10.1109/20.376275), the magnetostatic and
eddy-current analysis by [Bíró, Preis, and Richter
(1996)](https://ieeexplore.ieee.org/document/558631), and the high-order tree
gauge discussion in [the tree-gauge review](https://doi.org/10.3390/j5010004).
For spline and mortared spaces, see [Kapidani, Merkel, Schöps, and Vázquez
(2022)](https://doi.org/10.1016/j.cma.2022.114654) and the
[open preprint](https://arxiv.org/abs/2110.15860). The partial second-order
edge-element reduction is illustrated by [Bíró and Preis
(2015)](https://www.compumag.org/Proceedings/2015_Montreal/papers/OA7-1.pdf).
