# Source-conforming quadratic-level triangle quadrature

`fortfem_triangle_level_quadrature` integrates a vector callback on the reference
triangle `xi >= 0`, `eta >= 0`, `xi + eta <= 1`. It splits the domain at prescribed
levels of a P1/P2 scalar. This resolves piecewise source laws whose breakpoints
cross a finite element; increasing a fixed unsplit rule can miss those branches.
The caller retains ownership of the physical source, geometry determinant and
component error budget.

## Interface

```fortran
type(triangle_level_workspace_t) :: workspace
call initialize_triangle_level_workspace(workspace, nvalues, max_levels, status)
call triangle_level_coefficients(order, nodal_values, coefficients, status)
call integrate_triangle_levels(coefficients, levels, callback, epsabs, epsrel, &
    workspace, integral, error, status, message)
```

- Coefficients are ordered `[constant, xi, eta, xi**2, xi*eta, eta**2]`.
- Order one uses three vertex values. Order two uses six values, with nodes four,
  five and six on edges `(1,2)`, `(2,3)` and `(3,1)`. The nodal helper centers on
  node one before generating derivative coefficients and restores its constant.
- `levels` is finite and sorted nondecreasingly; repeated levels are accepted.
- The callback has arguments `point(2)`, `values(:)` and integer `status`, with
  respective intents `in`, `out` and `out`. `values` has `nvalues` components.
  A callback must give repeatable values for the same point throughout a call.
- `epsabs`, `integral` and `error` each have `nvalues` components. The stopping
  condition is `error <= epsabs + epsrel*abs(integral)` component by component.
  Use an absolute budget for components that can cancel or have zero integrals.
- A successful call exposes `workspace%points(:,1:npoints)` and positive
  `workspace%weights(1:npoints)`. Their measure is the reference triangle;
  physical geometry factors belong in the callback. These arrays can be read
  to evaluate additional components on the accepted trace. Those extra
  components do not inherit an error qualification automatically.

Initialize the workspace once and reuse it. Integration allocates no arrays.
Optional initialization limits default to `max_panels=256`,
`max_evaluations=200000` and `max_points=65536`; `max_levels` limits the input
knots. `nevaluations` counts actual callback invocations, including rejected
rules. The accepted trace is materialized without callback reevaluation.

Nonzero status means failure: invalid input/workspace, failed or nonfinite
callback values, or an exhausted evaluation/panel/trace/root-resolution budget.
Failure clears the integral and `npoints`; a partial trace must not be used.
An available error estimate describes the last completed trial, not a failed
partial trial. Workspace control fields and allocated arrays must not be
modified between initialization and integration.

## Algorithm and limits

The quadratic restrictions to `eta=0` and `eta=1-xi`, and the discriminant of
the inner quadratic, determine outer topology events. Each outer point splits
its inner interval at all level roots. FortNum provides bracketed root finding
and Gauss nodes. Native FortSym generates the scalar restrictions, derivatives,
nodal coefficient maps and interval measures; checked-in runtime outputs need
no CAS and introduce no external dependency.

Ordinary intervals use outer Gauss 8/16 and inner Gauss 4/8. Intervals adjoining
discriminant events use a generated sine-squared endpoint map with outer Gauss
16/32 to resolve square-root tangencies. Difficult smooth inner branches first
upgrade to Gauss 8/16, then subdivide; their outer pair is also upgraded. The
component estimate sums inner pair differences and outer pair differences
before global refinement. Accepted high-rule weights remain positive.

The estimate is an embedded-rule numerical diagnostic, **not a rigorous error
enclosure for an arbitrary callback**. Integrands with additional unresolved
jumps or rapid variation need their own qualified knot set and independent
refinement/oracle. Extremely ill-conditioned level roots and floating-point
limits may exhaust the explicit budget. Common gauges must be representable:
centering cannot recover nodal information already rounded away by the caller.

The module supplies primal integration and a retained trace, not a moving-level
JVP/VJP. For a continuous piecewise source, cancellation of moving-interface
terms is a separate mathematical argument. A caller must independently qualify
its continuous-integral derivative; it must not claim differentiation of the
finite-node adaptive algorithm or apply that argument to discontinuous sources.

## Independent checks

`test_triangle_level_quadrature` checks exact simplex moments, arbitrary P1/P2
nodal quadratics and representable large gauges, vector cancellation, horizontal
quadratic cuts, affine wedges, disk/quarter-disk/annulus polar moments, continuous
kinks and discontinuous slopes, logarithmic rational integrals, trace replay
and explicit failures. In particular, `eta**2 > 1/4` has area `1/8`, and
`max(eta**2-1/4,0)` integrates to `5/192`.

`test_triangle_level_conditioning` checks coefficient scaling from `1e-20` to
`1e20`, a representable `2**40` gauge, seven shrinking circles, component
estimates against independent analytical values and a floating-point floor,
and positive trace replay. At absolute `1e-13` and relative `1e-11`, these
circles each require 2,880 callbacks and three panels. An earlier affine-endpoint
prototype required 56,736 callbacks on the radius-one-eighth circle at absolute
`1e-12` and relative `1e-11`; the generated tangent map requires 2,880 at the
same budget. This measures quadrature work, not equilibrium solver speed.

FPM runs both tests directly. A standalone CMake/CTest source consumer is in
`test/fixtures/triangle_level_cmake`; configure it with an existing pinned
`FORTNUM_DIR`. FortFEM's existing C API and triangle compatibility CMake targets
do not expose this Fortran interface. A consumer enumerating its sources must
include the reference products, triangle-level generated geometry and the new
quadrature module alongside its existing FortNum root/quadrature/status modules.
