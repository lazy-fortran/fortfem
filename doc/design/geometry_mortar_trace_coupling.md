# Geometry-aware mortar trace coupling

`assemble_geometry_mortar_trace_coupling` is the small physical-measure layer
between a geometry provider and an independently discretized mortar pair. It
accepts test and trial trace values on the same fixed reference quadrature,
positive reference quadrature weights, and the positive surface Jacobian at
each quadrature row. It forms

\[
    w_q^{\mathrm{phys}} = w_q^{\mathrm{ref}} J_\Gamma(q),
    \qquad
    M_{ij}=\sum_q w_q^{\mathrm{phys}} T_i(q)S_j(q).
\]

The output `physical_weights` is intentionally exposed. It is an independent
surface-measure ledger that lets callers compare a NURBS, Fourier, cut-surface,
or panel geometry map against an analytical area or overlap measure before
assembling a large FEM/BEM/DtN block.

The module does not evaluate coordinates, normals, knots, or physical laws.
Those remain caller-owned geometry contracts. A NURBS sampler can use
`evaluate_nurbs_surface_geometry` to obtain the two tangent columns and pass
their norm as `surface_jacobian`; a panel or cut-cell provider can pass its own
metric using the same API. This keeps geometry representation independent of
the mortar basis and works for scalar, tangential-vector, normal-flux, and
Fourier/IGA traces.

## Differentiation

The JVP differentiates trace values, reference weights, and the surface metric
with a fixed quadrature topology:

\[
 \dot w_q = \dot w_q^{\mathrm{ref}}J_\Gamma(q)
            +w_q^{\mathrm{ref}}\dot J_\Gamma(q).
\]

The VJP first accumulates the cotangent of the physical weights from the cross
mass and then applies the product reverse:

\[
 \bar w_q^{\mathrm{ref}}=\bar w_qJ_\Gamma(q),\qquad
 \bar J_\Gamma(q)=\bar w_qw_q^{\mathrm{ref}}.
\]

Changing the owner of a quadrature row, adding/removing a cut row, or changing
the surface orientation is a discrete topology event and requires rebuilding
the fixed contract. Nonpositive metrics and reference weights are rejected.

## Verification

`test_geometry_mortar_trace_coupling` builds the cross-mass with an independent
triple loop, checks the physical-measure output, compares the complete JVP to
central reassembly, and checks the real VJP dot-product identity. It also
rejects a zero surface metric. This test contains no mesh reader or plasma
model and therefore serves as a compact geometry-to-mortar oracle for later
arbitrary multipatch IGA composition.
