# 3-D CGL pressure tensor gallery

This neutral manufactured gallery is the physical-first 3-D fixture for a
tensor-valued pressure/stress contract.  At every sample it evaluates the
gyrotropic CGL-shaped tensor

\[
  P = p_\perp I + (p_\parallel-p_\perp)\,b b^T,
  \qquad |b|=1,
\]

with a spatially varying, oblique magnetic direction.  The first figure shows
the pressure trace as 3-D coloured samples, blue oblique-`B` arrows, orange
surface-traction arrows `P n`, and the cube edges.  A later plot contains
component diagnostics only; the solution/traction visualization is intentionally
first.

The executable uses the public `fortfem_api` CGL tensor, traction, and
pressure-work operators.  Its in-program oracle independently reconstructs the
projector contraction and product rules, and checks tensor/traction/work JVP and
VJP adjoint identities.  `solution.csv` contains the complete 3-D sample field,
including tensor components, normals, tractions, and the velocity-gradient
contraction.  `diagnostics.csv` records all oracle errors and the small kernel
benchmark.

This is a constitutive and differentiation foundation, not an MHD equilibrium
solver: it contains no EQDSK reader, closure, force-balance iteration, or
plasma-specific input.  The CGL form follows the standard gyrotropic pressure
decomposition used in collisionless/anisotropic-fluid models; the numerical
fixture itself is independently manufactured.

Outputs:

- `solution.png` and `solution_3d.png`: pressure trace, oblique direction, and
  boundary traction in 3-D;
- `pressure_components.png`: diagnostic component plot;
- `solution.csv`: reproducible 3-D field samples;
- `diagnostics.csv`: contraction, JVP, VJP, positivity, and timing checks;
- `benchmark.txt` and `gallery_sequence.txt`.
