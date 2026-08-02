# Fixed-topology weighted flux-surface average

`fortfem_flux_surface_average` is a geometry-neutral reduction for values that
have already been sampled on a fixed set of surface quadrature points.  It is a
building block for nested-surface, Fourier--FEM, FEM--BEM, and external
equilibrium clients; it does not label flux surfaces, construct coordinates,
or impose a plasma closure.

For a scalar or vector sample (s_{cq}), quadrature weight (w_q), and an
optional positive geometry/volume factor (g_q), the contract is

\[
D = \sum_q w_q g_q,
\qquad
\bar{s}_c = \frac{1}{D}\sum_q w_q g_q s_{cq}.
\]

The sample array is component-major, `(component, sample)`.  A scalar is
represented by one component, so scalar and vector diagnostics use the same
interface.  If `measure_factors` is absent, (g_q=1).  The denominator (D)
is returned explicitly so callers can report the surface measure and detect a
bad geometry reduction instead of silently normalizing it away.

## Derivative contract

The JVP accepts sample, quadrature-weight, and (when supplied) geometry-factor
directions.  With (e_q=w_qg_q),

\[
\dot D = \sum_q \dot e_q,
\qquad
\dot{\bar{s}}_c =
  \frac{\sum_q e_q\dot{s}_{cq} + \dot e_q s_{cq}
        - \bar{s}_c\dot D}{D},
\]

where (dot e_q=\dot w_qg_q+w_q\dot g_q).  The VJP accepts cotangents for
the average and denominator and returns cotangents for samples, quadrature
weights, and optionally geometry factors.  It uses the real Euclidean
dot-product convention; complex clients must apply their declared real-part
convention at their complex boundary.

The implementation rejects empty sample sets, shape mismatches, non-finite
values, non-positive quadrature weights or geometric factors, and a non-finite
or non-positive denominator.  The topology and quadrature points are fixed
while differentiating.  Changing the number or ordering of points is a
discrete rebuild, not a differentiable operation.

## Provenance and scope

The normalized weighted reduction is the common numerical operation behind
surface averages in nested-surface equilibrium workflows.  DESC documents
its nested-surface coordinate and objective/constraint workflow at
[`desc-docs.readthedocs.io`](https://desc-docs.readthedocs.io/en/stable/), while
the VMEC variational context is described by the
[VMEC documentation](https://w3.pppl.gov/ntcc/VMEC/).  FortFEM deliberately
implements only the reusable value/JVP/VJP reduction.  Surface labels,
Fourier/radial bases, Boozer transforms, profiles, force balance, and
free-boundary physics remain external client responsibilities.

The independent test uses manufactured scalar and vector samples, central
differences for the JVP, and the real dot-product identity for the VJP.  It
also checks scalar one-component use, the optional factor path, and rejection
of empty and non-positive measures.
