---
title: Surface-current trace contract
---

# Surface-current trace contract

`fortfem_surface_current` turns an oriented vector trace jump into a generic
surface-current field and an integrated ledger. The normal points from the
minus side to the plus side:

\[
\mathbf K=\mathbf n\times(\mathbf H^+-\mathbf H^-),
\qquad
\bar{\mathbf K}=\sum_q w_q\mathbf K_q.
\]

This is the orientation-safe Ampere trace algebra. It does not choose a
constitutive law, a plasma closure, a wall model, or a material boundary
condition.

## API

```fortran
call assemble_interface_surface_current( &
    plus_field, minus_field, normals, surface_weights, &
    surface_current, integrated_current, status)
call assemble_interface_surface_current_jvp( &
    plus_field, minus_field, normals, surface_weights, plus_dot, minus_dot, &
    normals_dot, surface_weights_dot, surface_current_dot, &
    integrated_current_dot, status)
call assemble_interface_surface_current_vjp( &
    plus_field, minus_field, normals, surface_weights, surface_current_bar, &
    integrated_current_bar, plus_bar, minus_bar, normals_bar, &
    surface_weights_bar, status)
call assemble_surface_current_junction_balance( &
    junction_incidence, manifold_current, junction_balance, &
    global_balance, status)
call assemble_surface_current_junction_balance_jvp( &
    junction_incidence, manifold_current, manifold_current_dot, &
    junction_balance_dot, global_balance_dot, status)
call assemble_surface_current_junction_balance_vjp( &
    junction_incidence, manifold_current, junction_balance_bar, &
    global_balance_bar, manifold_current_bar, status)
call assemble_surface_current_loop_constraints( &
    loop_basis, manifold_current, target_loop_current, residual, status)
call assemble_surface_current_loop_constraints_jvp( &
    loop_basis, manifold_current, target_loop_current, manifold_current_dot, &
    target_loop_current_dot, residual_dot, status)
call assemble_surface_current_loop_constraints_vjp( &
    loop_basis, manifold_current, target_loop_current, residual_bar, &
    manifold_current_bar, target_loop_current_bar, status)
```

The primal validates positive quadrature weights, finite three-component
traces, and unit normals. The JVP and real VJP keep the quadrature topology and
orientation fixed. The integrated-current cotangent is accumulated into each
pointwise current cotangent before the cross-product reverse rule is applied.

Reversing the interface convention (swap the two traces and negate the normal)
leaves the physical \(\mathbf K\) unchanged. This is an independent sign
oracle, not a property inferred from a caller's boundary law.

The junction ledger consumes the integer boundary incidence supplied by the
neutral internal-manifold graph and an integrated current per manifold. It
returns a vector balance at every open junction and the global sum. Closed or
boundaryless manifolds have zero incidence columns and therefore cannot create
a spurious net current. Its topology is fixed for differentiation; the JVP
and VJP act on the manifold-current values only.

Closed cycles or externally prescribed current loops use an integer loop basis
`loop_basis(loop, manifold)`. The loop residual is

\[
  r_\ell = \sum_m B_{\ell m} K_m - K_{\ell,\mathrm{target}}.
\]

This is a neutral constraint block: the target may be a flux, current, or
other application-owned normalization with declared units. Its fixed-topology
JVP differentiates manifold and target values, while the VJP is the transpose
of the same integer map. The loop basis is not differentiated across a graph
rebuild or cycle-basis change.

## Independent verification

`test_surface_current` checks a two-point analytical jump and weighted ledger,
orientation reversal, central-difference JVP, real dot-product VJP, and
malformed normal/weight rejection. Surface quadrature geometry and current
conservation at interface edges remain separate geometry contracts.
`test_surface_current_balance` checks open, closed, and malformed incidence,
global conservation, and the real ledger adjoint identity.
`test_surface_current_constraints` checks an independent two-cycle residual,
the fixed-topology JVP, the real adjoint identity, and shape rejection.
