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
call evaluate_shifted_enriched_basis( &
    base_values, level_values, anchor_values, enriched_values, status)
call evaluate_shifted_enriched_basis_jvp( &
    base_values, level_values, anchor_values, base_dot, level_dot, anchor_dot, &
    enriched_dot, status)
call evaluate_shifted_enriched_basis_vjp( &
    base_values, level_values, anchor_values, enriched_bar, base_bar, &
    level_bar, anchor_bar, status)
call evaluate_shifted_vector_enriched_basis( &
    base_values, level_values, anchor_values, enriched_values, status)
call evaluate_shifted_vector_enriched_basis_jvp( &
    base_values, level_values, anchor_values, base_dot, level_dot, anchor_dot, &
    enriched_dot, status)
call evaluate_shifted_vector_enriched_basis_vjp( &
    base_values, level_values, anchor_values, enriched_bar, base_bar, &
    level_bar, anchor_bar, status)
```

The activation sign is a fixed-topology discrete choice. Away from `phi=0`
the JVP and VJP are exactly zero. A zero level or anchor value is reported as
a topology event instead of being assigned an arbitrary derivative. This
makes the piecewise-smooth differentiation contract explicit before cut
geometry and enrichment activation are composed.

`evaluate_shifted_enriched_basis` applies the product rule
`N_i*(H(phi)-H(phi_i))` to a caller-supplied base value. Its JVP propagates
base-value increments and its VJP returns the corresponding base cotangent;
the level-set cotangents remain zero on a fixed sign pattern. This keeps
partition-of-unity composition explicit without pretending that enrichment
activation or cut connectivity is differentiable across a topology event.

`evaluate_shifted_vector_enriched_basis` applies the same scalar activation
componentwise to a rank-two array of vector basis values. It is deliberately
agnostic about the vector representation, so the caller can supply covariant
or contravariant Piola values from an H(curl) or H(div) space. The JVP uses
the componentwise product rule and the VJP returns the base-value cotangent;
the scalar level-set VJP remains zero on a fixed topology.

`test_heaviside_enrichment` checks the independent sign oracle, fixed-sign
zero derivative identities, and topology-event rejection. The companion
`test_shifted_enriched_basis` checks the scalar product oracle and real
adjoint identity. `test_shifted_vector_enriched_basis` repeats those checks
for a two-component vector base and rejects a zero-level topology event.
