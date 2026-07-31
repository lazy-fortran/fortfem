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
```

The primal validates positive quadrature weights, finite three-component
traces, and unit normals. The JVP and real VJP keep the quadrature topology and
orientation fixed. The integrated-current cotangent is accumulated into each
pointwise current cotangent before the cross-product reverse rule is applied.

Reversing the interface convention (swap the two traces and negate the normal)
leaves the physical \(\mathbf K\) unchanged. This is an independent sign
oracle, not a property inferred from a caller's boundary law.

## Independent verification

`test_surface_current` checks a two-point analytical jump and weighted ledger,
orientation reversal, central-difference JVP, real dot-product VJP, and
malformed normal/weight rejection. Surface quadrature geometry and current
conservation at interface edges remain separate topology contracts.
