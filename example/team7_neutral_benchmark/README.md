# TEAM-7-shaped neutral eddy-current benchmark

This solution-first fixture is a small, license-safe 2-D analogue of the
coil/plate eddy-current layout used in the TEAM benchmark family.  A complex
manufactured stream function defines a divergence-free magnetic field,
`B = (dA/dy,-dA/dx)`, and the corresponding complex Ampere current
`J_z = -laplacian(A)`.  The first plot, `solution.png`, shows the field
magnitude, vector directions, coil outline, conducting plate, and slot.  The
additional `current.png` and `probe.png` plots expose the induced-current
response and a probe cut.  `solution.csv`, `probe.csv`, and `diagnostics.csv`
are generated beside the plots.

The geometry, source amplitudes, and material response are intentionally
manufactured.  This example does **not** embed a TEAM reader, nonlinear B-H
curve, proprietary geometry, or external reference data; an exact benchmark
comparison belongs in the separate benchmark-data repository when its source
terms and license permit redistribution.

The provenance target is the [TEAM workshop report](https://www.osti.gov/servlets/purl/7179128)
and the public [TEAM problem catalogue](https://www.osti.gov/biblio/7179128).  The
fixture is intended to exercise FortFEM's 2-D vector-field plotting, analytic
manufactured forcing, and diagnostic hooks before a specialized solver is
added.

Run it with:

```text
fo run --example team7_neutral_benchmark
```
