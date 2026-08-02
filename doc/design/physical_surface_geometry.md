# Physical surface geometry sampler

`sample_physical_surface_geometry` is the neutral geometry contract between a
surface-map provider and FEM, BEM, DtN, or mortar assembly.  For every fixed
reference quadrature coordinate (\xi_q=(\xi_q^1,\xi_q^2)), a provider supplies
the evaluated physical map (x_q\in\mathbb R^3) and its two tangent columns

\[
    t_{1,q}=\partial_{\xi^1}x(\xi_q),\qquad
    t_{2,q}=\partial_{\xi^2}x(\xi_q).
\]

The sampler returns the physical coordinates, the positive measure

\[
    J_\Gamma(q)=\|t_{1,q}\times t_{2,q}\|>0,
\]

and the oriented unit normal

\[
    n_q=s\frac{t_{1,q}\times t_{2,q}}{J_\Gamma(q)},
    \qquad s\in\{-1,+1\}.
\]

The map evaluation is deliberately provider-owned.  Thus the same contract
can be fed by NURBS/IGA, Fourier-toroidal, panel, or cut-cell geometry without
embedding a particular representation in the mortar or boundary operators.
The reference coordinates are retained as a fixed row identity and are
validated, but are not differentiated: coordinate motion must be included by
the upstream map provider in `map_points_dot` and `map_tangents_dot`.

## Differentiation

The JVP uses

\[
    \dot J_\Gamma=
    \frac{(t_1\times t_2)\cdot
      (\dot t_1\times t_2+t_1\times\dot t_2)}{J_\Gamma},
\]

and differentiates the normalized normal analytically.  The VJP reverses the
normalization and cross product, returning map-point and both tangent-column
cotangents in one sweep.  The orientation sign and quadrature topology are
discrete inputs and have no derivative.

Non-finite data, incompatible ranks, orientation values other than ±1, and
degenerate tangent pairs are rejected before normalization.  This prevents a
zero or NaN surface measure from silently entering a physical cross-mass.

## Verification

`test_physical_surface_geometry` uses an independent cross-product oracle,
checks the normal orientation and positive measure at three samples, compares
all JVP outputs with central re-evaluation, verifies the real VJP dot-product
identity, and rejects both degenerate tangents and invalid orientation signs.
