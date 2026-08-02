# TEAM-13-shaped neutral electromagnetic benchmark

This is a small, solution-first foundation fixture shaped like the TEAM 13
coil/channel/plate arrangement. It is not a TEAM reader or a claim to solve
the nonlinear benchmark exactly: the geometry, source, and material response
are supplied manufactured arrays so that the FortFEM gallery remains
license-safe and fast. Exact TEAM geometry, B-H data, probe data, and reference
curves belong in the sister benchmark-data repository when redistribution is
allowed.

The first output is `solution.png`, showing the manufactured magnetic-field
magnitude, vector field, and benchmark-shaped bodies. `probe.png` shows a
probe curve. CSV and diagnostic files are generated under
`output/example/team13_neutral_benchmark`; no images are committed.

The provenance target is the [TEAM workshop report](https://www.osti.gov/servlets/purl/7179128)
and the public [TEAM 13 reproduction description](https://docs.feelpp.org/toolboxes/latest/maxwell/Tws/index.html).

Run with:

```text
fo run --example team13_neutral_benchmark
```
