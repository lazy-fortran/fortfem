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
`assemble_surface_vector_delta_load` now provides the corresponding
work-conjugate tangential surface-current pairing.  Parameter derivatives and
surface measures will build on these contracts in later slices.
