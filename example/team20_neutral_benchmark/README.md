# TEAM-20-shaped neutral magnetostatic benchmark

This is a small, solution-first foundation fixture shaped like TEAM Problem
20: a three-dimensional solenoid with a central pole and an enclosing yoke.
The field is generated from a smooth manufactured vector potential, so the
curl construction supplies an analytical divergence-free solution and a fast
plotting/diagnostics path.

The fixture is deliberately **not** an exact TEAM-20 reproduction. It does
not contain a TEAM reader, the nonlinear steel B--H curve, the workshop mesh,
force measurements, or redistributed reference data. Those inputs belong in
the sister benchmark-data repository after provenance and license review. The
supplied outlines are only a visual geometry surrogate. `solution.png` is a
2-D pole/yoke cut with field vectors; `solution_3d.png` is a 3-D solenoid
surface with sparse vector segments; `probe.png` is a one-dimensional cut.
CSV and diagnostic files are generated under
`output/example/team20_neutral_benchmark`; images are not committed.

The provenance targets are the [TEAM workshop report](https://www.osti.gov/servlets/purl/7179128),
the [TEAM-20 static-force description](https://www.simscale.com/docs/validation-cases/team-20-magnetostatics/),
and the [Bíró/Ostergaard TEAM-20 proceedings reference](https://ansyshelp.ansys.com/public/Views/Secured/corp/v252/en/ans_vm/Hlp_V_VM233.html).

Run it with:

```text
fo run --example team20_neutral_benchmark
```
