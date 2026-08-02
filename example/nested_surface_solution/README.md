# Nested toroidal surface solution

This physical-first gallery fixture evaluates FortFEM's neutral nested-surface
map for a circular torus and colors the physical surface with a smooth
manufactured scalar solution. Dark parameter lines make both periodic
directions visible. The geometry is produced by the same Fourier/radial map
available to equilibrium-code adapters; the example does not select an
equilibrium model, coordinate convention, or file format.

The first plot is the physical three-dimensional solution. The later plots
show its periodic parameter-space representation and the axis-regular radial
profile used to modulate it. The program checks the evaluated Cartesian
coordinates against the closed-form circular-torus embedding before plotting.

Generated artifacts:

- `nested_surface_solution_3d.png`: scalar solution on the physical torus;
- `nested_surface_solution_parameter_2d.png`: periodic surface field in
  `(theta,zeta)` coordinates;
- `nested_surface_radial_profile_1d.png`: axis-regular radial amplitude;
- `benchmark.txt`: sample counts, oracle error, timings, and provenance.

Generated images and data remain under `output/example/` and are not committed.
