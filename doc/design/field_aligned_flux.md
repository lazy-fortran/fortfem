# Field-aligned anisotropic flux

`evaluate_field_aligned_flux` is the generic pointwise constitutive block for
strongly anisotropic transport and wave operators.  For a unit magnetic
direction (b) and gradient (g), it evaluates

\[
  F = k_\perp g + (k_\parallel-k_\perp)b(b\cdot g).
\]

The non-negative coefficients allow the isotropic limit
(k_\parallel=k_\perp), the perpendicular-off limit, and strongly separated
ratios without changing the API.  The wrapper validates vector layout,
coefficient signs, and the unit-direction contract; the value, JVP, and VJP
products are emitted by the pinned FortSym generator
`gen_field_aligned_flux_products`.

The focused test uses an independent tensor contraction, a central-difference
JVP on a tangent direction, and a real dot-product VJP identity.  This is a
pointwise constitutive ingredient, not a complete assembled diffusion or
Braginskii operator; quadrature, field-line maps, and anisotropy-aware
preconditioners remain higher-level concerns.

## Provenance

The decomposition follows the standard parallel/perpendicular transport tensor
used in anisotropic plasma-fluid and MHD models, and is recorded in the
FortFEM roadmap alongside the GRILLIX/BOUT++/PARALLAX ingredient lineage.  The
implementation is independently generated and contains no downstream code.
