# Mixed elasticity wave

This example is a manufactured mixed elastic-wave modal problem. It uses the
structure-preserving mixed midpoint step for displacement/velocity,
reconstructs the physical displacement and stress fields in the original
one-dimensional bar view and in a small two-dimensional contour/quiver view,
and checks the analytical modal solution, reversibility, energy invariant, and
the neutral mixed elasticity residual. The physical plots come first; the
energy and error diagnostics follow them. The 2-D sample values are also
written to `mixed_elasticity_solution_2d.csv` for an independent array oracle.

It is a numerical foundation example, not a plasma or production elasticity
solver. Element spaces, weak symmetry maps, and material laws remain caller
owned and are exposed through the public contracts documented in
`mixed_elasticity_residual` and `elasticity_symmetry_constraint`.
