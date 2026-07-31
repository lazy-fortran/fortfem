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
call evaluate_vector_enrichment_differential_3d( &
    base_values, base_gradient, activation, activation_gradient, &
    enriched_values, curl_values, divergence, status)
call evaluate_vector_enrichment_differential_3d_jvp( &
    base_values, base_gradient, activation, activation_gradient, &
    base_values_dot, base_gradient_dot, activation_dot, activation_gradient_dot, &
    enriched_values_dot, curl_values_dot, divergence_dot, status)
call evaluate_vector_enrichment_differential_3d_vjp( &
    base_values, base_gradient, activation, activation_gradient, &
    enriched_values_bar, curl_values_bar, divergence_bar, base_values_bar, &
    base_gradient_bar, activation_bar, activation_gradient_bar, status)
call evaluate_blending_corrected_enrichment( &
    base_values, enriched_mask, enrichment_values, corrected_values, status)
call evaluate_blending_corrected_enrichment_jvp( &
    base_values, enriched_mask, enrichment_values, base_dot, enrichment_dot, &
    corrected_dot, status)
call evaluate_blending_corrected_enrichment_vjp( &
    base_values, enriched_mask, enrichment_values, corrected_bar, base_bar, &
    enrichment_bar, status)
call evaluate_vector_blending_corrected_enrichment( &
    base_values, enriched_mask, enrichment_values, corrected_values, status)
call evaluate_vector_blending_corrected_enrichment_jvp( &
    base_values, enriched_mask, enrichment_values, base_dot, enrichment_dot, &
    corrected_dot, status)
call evaluate_vector_blending_corrected_enrichment_vjp( &
    base_values, enriched_mask, enrichment_values, corrected_bar, base_bar, &
    enrichment_bar, status)
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

For a three-dimensional vector basis, the differential companion makes the
de Rham product terms explicit:

```
curl(psi*b) = psi*curl(b) + grad(psi) x b
 div(psi*b) = psi*div(b)  + grad(psi) . b
```

`evaluate_vector_enrichment_differential_3d` accepts caller-owned base values,
gradients, activation, and activation gradients. Its JVP differentiates all
four inputs and its VJP returns all four cotangents. This is a diagnostic and
assembly primitive for covariant/contravariant Piola values, H(curl), H(div),
and IGA callers; it does not claim that an arbitrary enrichment preserves an
exact sequence. The explicit correction terms are the independent quantity
that a commuting-space implementation must reproduce or intentionally report
as a jump contribution.

`evaluate_blending_corrected_enrichment` implements the corrected-XFEM ramp
without choosing an element or geometry representation. For a fixed logical
node mask (a_i), it forms

```
r(x) = sum_i a_i N_i(x),       Psi_corr(x) = r(x) Psi(x).
```

Thus a standard element has zero enrichment, a fully enriched element
reproduces `Psi`, and a blending element transitions through the partition of
unity. Its JVP/VJP differentiate the two products while treating the mask as
discrete topology. This is the composition described by [Fries' corrected
XFEM construction](https://doi.org/10.1002/nme.2259); cut-cell activation and
rank/conditioning diagnostics stay outside this primitive.

The vector companion applies the same ramp to every component. Its reverse
action contracts the vector cotangent with the vector enrichment before
accumulating the base-function cotangent, so it is suitable for covariant or
contravariant Piola values and for IGA coefficient arrays without changing
the fixed-mask topology contract.

`test_heaviside_enrichment` checks the independent sign oracle, fixed-sign
zero derivative identities, and topology-event rejection. The companion
`test_shifted_enriched_basis` checks the scalar product oracle and real
adjoint identity. `test_shifted_vector_enriched_basis` repeats those checks
for a two-component vector base and rejects a zero-level topology event.
`test_vector_enrichment_differential_3d` checks the independent curl/divergence
product formulas, central-difference JVP, real adjoint identity, and output
shape validation.
`test_xfem_blending_correction` checks the ramp, full-enrichment limit, and
its real adjoint identity. `test_xfem_vector_blending_correction` checks the
componentwise ramp and vector reverse contraction.
