---
title: Complex boundary trace residual
---

# Complex boundary trace residual

`assemble_complex_boundary_trace_residual` provides the semantic boundary-port
contract used by free-boundary and interface clients:

\[
  r_n=w_n(Nu-g_n),\qquad r_t=w_t(Tu-g_t).
\]

`N` and `T` are caller-owned normal and tangential trace maps. Their targets
may be prescribed traces, jumps, surface-current data, or total-pressure data;
the module does not select a constitutive law, equilibrium normalization, or
geometry representation. FEM, BEM, DtN, PML, DG, IGA, and wall blocks can all
populate the maps.

Value, JVP, and VJP actions include trace maps, positive work weights, state,
and both supplied targets. VJPs use the real-part complex inner product. The
implementation is matrix-free and only uses caller-sized residual/cotangent
arrays. The independent test checks separate normal and tangential matrix
oracles, a complete product-rule difference, and the complex adjoint identity.
