# Mixed elasticity wave

This example is a manufactured one-dimensional elastic-bar modal problem. It
uses the structure-preserving mixed midpoint step for displacement/velocity,
reconstructs the physical displacement and stress fields, and checks the
analytical modal solution, reversibility, energy invariant, and the neutral
mixed elasticity residual. The first plot is the physical bar solution; the
energy and error diagnostics follow it.

It is a numerical foundation example, not a plasma or production elasticity
solver. Element spaces, weak symmetry maps, and material laws remain caller
owned and are exposed through the public contracts documented in
`mixed_elasticity_residual` and `elasticity_symmetry_constraint`.
