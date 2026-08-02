# Three-dimensional mixed-wave structure

This neutral gallery fixture advances three manufactured Cartesian oscillator
components with the common first-order mixed-wave midpoint block.  The first
plot is the physical (q_x,q_y,q_z) trajectory and overlays the independent
closed-form solution.  Component and energy plots follow as diagnostics.

The example is a time-integration foundation only: it contains no plasma,
material, damping, or boundary physics.  A caller can provide compatible
pressure/velocity, displacement/momentum, Maxwell, elasticity, or tensor-
pressure blocks to the same structure-preserving API.
