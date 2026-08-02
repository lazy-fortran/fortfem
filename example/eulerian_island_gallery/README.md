# Eulerian island/separatrix foundation

This small gallery samples the analytic slab-island flux

\[
\psi(x,y)=(x^2-a)^2+y^2,\qquad
\mathbf B=(\partial_y\psi,-\partial_x\psi),
\]

on a fixed Eulerian domain.  The two closed O-point islands and the white
`psi = psi_X` figure-eight expose the separatrix in the first plot;
normalized tangent-field arrows are overlaid for orientation.  A second plot
shows the caller-owned manufactured force residual magnitude.

The example calls `assemble_eulerian_nonnested_residual` and its JVP with
caller-owned force/divergence samples.  It deliberately selects no plasma
closure.  Outputs are written below
`output/example/eulerian_island_gallery/`:

- `island_flux_solution_2d.png` (physical solution and field arrows)
- `island_flux_diagnostics_2d.png` (residual diagnostic)
- `island_flux.csv` (point samples)
- `provenance.json` (grid, event, and residual metadata)
