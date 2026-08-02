# Direct force-balance campaign on an analytic torus

This small, closure-neutral campaign exercises FortFEM's public direct
force-balance objective on a manufactured state sampled on an exact torus.
The torus and its scalar state are analytical; the vector residual and
positive surface weights are caller-owned samples.  No equilibrium closure,
plasma profile, reader, or optimizer is selected here.

The first image is deliberately a physical view: a parametric torus carrying
the manufactured state and short Cartesian force-residual vectors.  A
parameter-space field map and objective derivative diagnostics follow it.

Generated, uncommitted artifacts are:

- `direct_force_torus_solution_3d.png`: physical torus and force vectors;
- `direct_force_state_2d.png`: scalar manufactured state on `(theta, phi)`;
- `direct_force_objective_diagnostics_1d.png`: objective/JVP and independent
  finite-difference diagnostics;
- `direct_force_surface.csv` and `direct_force_diagnostics.csv`: numerical
  plot and oracle inputs;
- `benchmark.json`: objective, derivative errors, sample counts, and timing
  provenance.

The campaign is a tiny operator contract rather than a performance claim for
a production equilibrium code.
