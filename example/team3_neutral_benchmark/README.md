# TEAM-3-shaped neutral magnetostatic benchmark

This example is a small, solution-first foundation fixture shaped like the
classic TEAM-3 C-core/air-gap magnetostatic problem.  It intentionally uses
only supplied manufactured arrays: a smooth vector potential, a positive
relative-permeability surrogate, a compact source envelope, and a polygonal
core/coil layout.  The field is

\[
 B=(\partial_y A_z,-\partial_x A_z),\qquad
 J_z=-\Delta A_z,
\]

so the divergence-free and Ampere identities are available as independent
oracles without embedding a TEAM reader, nonlinear B-H curve, restricted
geometry, or reference measurement data.  This is **not** an exact TEAM-3
reproduction; exact benchmark input belongs in the separately licensed
benchmark-data repository.

The first plot, `solution.png`, shows the manufactured magnetic-field
magnitude, vector directions, C-core and coil layout.  `solution_3d.png` is a
3-D height surface of the same field with lifted vector arrows.  `probe.png`
shows a cut through the air gap.  `solution.csv`, `material.csv`,
`geometry.csv`, `probe.csv`, and `diagnostics.csv` contain reproducible arrays
and diagnostics under `output/example/team3_neutral_benchmark`.

Provenance: the shape is motivated by the public TEAM benchmark catalogue and
workshop report ([TEAM catalogue](https://www.osti.gov/biblio/7179128),
[workshop report](https://www.osti.gov/servlets/purl/7179128)).  No source
arrays or measurement data are copied from those works.

Run with:

```text
fo run --example team3_neutral_benchmark
```
