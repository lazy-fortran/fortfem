---
title: Enrichment-support activation
---

# Enrichment-support activation

`fortfem_enrichment_support_activation` turns a caller-owned level-set sample
and CSR support map into the fixed-topology activation data used by XFEM and
XIGA assembly.  The map is deliberately independent of a mesh or a basis
family: `support_offsets` is one-based and half-open, and each entry in
`support_nodes` refers to a level-set sample.

For support `s`, let

```
m_s = min(phi_i),  M_s = max(phi_i),  d_s = min(M_s, -m_s).
```

The support is active exactly when `m_s < 0 < M_s`.  The routine also returns
the unique owners of `m_s` and `M_s`, and a margin branch (`+1` for `M_s`,
`-1` for `-m_s`).  A zero sample, duplicate support entry, invalid CSR map,
or non-finite level set is rejected as a topology or input error.  The
activation mask itself is discrete; this is intentional because changing it
changes the algebraic space.

```fortran
call evaluate_enrichment_support_activation( &
    level_values, support_offsets, support_nodes, active_mask, support_min, &
    support_max, activation_margin, min_owner, max_owner, margin_branch, &
    status)
call evaluate_enrichment_support_activation_jvp( &
    level_values, support_offsets, support_nodes, min_owner, max_owner, &
    margin_branch, level_dot, support_min_dot, support_max_dot, &
    activation_margin_dot, status)
call evaluate_enrichment_support_activation_vjp( &
    level_values, support_offsets, support_nodes, min_owner, max_owner, &
    margin_branch, support_min_bar, support_max_bar, activation_margin_bar, &
    level_bar, status)
```

The JVP and VJP are valid only with fixed unique extrema and a fixed nonzero
margin branch.  They reject an owner change, an extremum tie, an invalid owner,
or a margin tie instead of silently differentiating through a topology event.
The VJP accumulates all support contributions into the shared global
level-set samples, so overlapping supports and componentwise XIGA activation
retain a normal real dot-product adjoint contract.

`test_enrichment_support_activation` provides the independent CSR sign/min/max
oracle, central-difference JVP checks, real dot-product VJP check, and
invalid-map/topology/tie rejection checks.  Cut-cell geometry and
Piola-aware conditioning remain separate contracts: this primitive only
defines support activation and its fixed-topology derivatives.
