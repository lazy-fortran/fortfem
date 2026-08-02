# Neutral nonlinear resistive-MHD island/wall gallery

This bounded manufactured fixture is the physical-first gallery seed for
MHD-14.  It uses the public eight-block nonlinear composition with a small
caller-owned state:

\[
  u=(a_{\rm island},i_{\rm wall},p') .
\]

The first plot shows a smooth manufactured island-flux profile and its
localized wall-current response.  The later plots show forward/reverse
continuation branches and the residual/input/dissipation ledger.  The
machine-readable CSV and JSON files retain the branch multiplicity,
hysteresis, path metric, residual norm, input power, dissipation, and timing.

This is deliberately not an equilibrium or resistive-MHD solver.  The
constitutive callbacks, state selection, geometry, closure, and continuation
policy remain external.  The fixture only verifies that a neutral client can
compose those callbacks and expose a physical state before diagnostics.

Outputs (generated and ignored by git):

- `island_wall_solution_1d.png`
- `island_wall_branches_1d.png`
- `island_wall_ledger_1d.png`
- `island_wall_solution.csv`
- `island_wall_continuation.csv`
- `benchmark.json`
