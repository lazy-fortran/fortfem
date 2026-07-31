---
title: Surface-delta weak loads
---

# Surface-delta weak loads

`assemble_surface_delta_load` assembles an explicit interface contribution

\[
  \ell_i=\int_\Gamma g\,\operatorname{tr}(v_i)\,dS
  \approx\sum_q T_{qi}\,w_q\,g_q.
\]

The trace basis, surface quadrature, and source stay separate from volume
assembly.  This is the algebraic primitive for a surface current, pressure
jump, or other distribution-valued source on a fitted interface; it does not
replace a cut-cell delta approximation.  Positive finite surface weights and
compatible arrays are required.

The focused test compares against an independent trace-transpose oracle and
checks incompatible quadrature sizes.  Vector/tangential current pairings and
`assemble_surface_vector_delta_load` provides the corresponding work-conjugate
tangential surface-current pairing. Both scalar and vector primitives expose
fixed-topology analytic JVP and real VJP actions, so a current-sheet residual
can be differentiated without finite-differencing quadrature or basis data:

```fortran
call assemble_surface_delta_load_jvp( &
    trace_basis, surface_weights, surface_source, trace_basis_dot, &
    surface_weights_dot, surface_source_dot, load_dot, status)
call assemble_surface_delta_load_vjp( &
    trace_basis, surface_weights, surface_source, load_bar, &
    trace_basis_bar, surface_weights_bar, surface_source_bar, status)
call assemble_surface_vector_delta_load_jvp( &
    tangential_trace_basis, surface_weights, surface_current, &
    tangential_trace_basis_dot, surface_weights_dot, surface_current_dot, &
    load_dot, status)
call assemble_surface_vector_delta_load_vjp( &
    tangential_trace_basis, surface_weights, surface_current, load_bar, &
    tangential_trace_basis_bar, surface_weights_bar, surface_current_bar, &
    status)
```

The JVPs apply the full product rule to basis, geometry weights, and source or
current values. The VJPs are the real transpose actions and are suitable for
adjoint or implicit-differentiation compositions. Topology, basis activation,
and quadrature connectivity remain fixed during these derivatives; a cut or
fitted geometry rebuild is a separate event.

`assemble_surface_triangle_measures_3d`
supplies oriented triangle areas and normals, including geometry JVP/VJP
products, for fitted surface quadrature.
