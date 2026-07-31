---
title: Shifted Heaviside enrichment
---

# Shifted Heaviside enrichment

`fortfem_heaviside_enrichment` provides the first XFEM/GFEM enrichment
primitive. For a level set `phi` and an anchor value `phi_i`, it returns

```
Psi_i = H(phi) - H(phi_i)
H(phi) = 1 for phi > 0, and 0 for phi < 0.
```

The shifted form preserves a zero value at the anchor and can therefore be
multiplied by a partition-of-unity basis without changing the un-enriched
interpolation point. The routine is deliberately geometry-agnostic; cut-cell
classification, support activation, and blending correction remain separate
layers.

```fortran
call evaluate_shifted_heaviside_enrichment( &
    level_values, anchor_values, values, status)
call evaluate_shifted_heaviside_enrichment_jvp( &
    level_values, anchor_values, level_dot, anchor_dot, values_dot, status)
call evaluate_shifted_heaviside_enrichment_vjp( &
    level_values, anchor_values, values_bar, level_bar, anchor_bar, status)
```

The activation sign is a fixed-topology discrete choice. Away from `phi=0`
the JVP and VJP are exactly zero. A zero level or anchor value is reported as
a topology event instead of being assigned an arbitrary derivative. This
makes the piecewise-smooth differentiation contract explicit before cut
geometry and enrichment activation are composed.

`test_heaviside_enrichment` checks the independent sign oracle, fixed-sign
zero derivative identities, and topology-event rejection.
