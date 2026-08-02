# Arbitrary multipatch signed DOF graphs

`build_multipatch_signed_dof_map` is the topology-only numbering layer for
fitted multipatch FEEC and IGA discretizations. It accepts a packed list of
patch-local degrees of freedom and interface relations of the form

\[
    x_{\mathrm{left}} = s\,x_{\mathrm{right}}, \qquad s\in\{-1,+1\}.
\]

The result is a signed local-to-global map compatible with the existing glued
FEEC dense and CSC assemblers. A negative entry means that the local basis
orientation is opposite to the chosen global representative. Interface
cycles are checked for consistency; an inconsistent cycle returns status 2,
while malformed offsets, endpoints, or signs return status 1.

The map is intentionally independent of geometry, spline knots, trace
quadrature, material laws, and interface physics. Callers construct those
objects separately and use this map to compose arbitrary patch graphs,
including periodic and cyclic graphs. The numbering is a fixed-topology
operation: shape derivatives do not differentiate this integer map. If a
topology event changes the graph, rebuild it and treat that event as a
piecewise-smooth boundary of the differentiable problem.

The implementation uses a signed disjoint-set forest, so it does not create a
dense identification matrix. The independent test checks the equivalence
classes and every declared orientation relation, and also exercises rejection
of an inconsistent cycle.
