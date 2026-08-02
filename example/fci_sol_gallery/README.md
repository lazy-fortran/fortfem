---
title: fci_sol_gallery Example
---

This neutral FCI/SOL fixture traces three helical field lines on a toroidal
surface and applies the same fixed-topology maps to FortFEM's conservative
parallel gradient and diffusion support operators. The first image is the
physical three-dimensional field-line solution; the second is a two-dimensional
field-coordinate-plane solution. The diagnostic plot and JSON/CSV ledgers come
after those solution views.

The example deliberately does not choose a plasma species, sheath closure,
transport model, equilibrium reader, or edge-physics code. It provides a small
geometry and operator contract that later SOL applications can consume. The
support pairing follows the PARALLAX field-coordinate-independent construction
described by [Stegmeir et al.](https://doi.org/10.1016/j.cpc.2015.09.016).

Outputs include:

- `fci_sol_field_lines_3d.png`: toroidal surface and traced field lines (the
  physical primary plot);
- `fci_sol_solution_2d.png`: manufactured scalar on the FCI plane grid;
- `fci_sol_gradient_diagnostic_1d.png`: support-gradient diagnostic against the
  analytical line derivative;
- `fci_sol_solution.csv` and `fci_sol_field_lines.csv`: reproducible values;
- `benchmark.json`: timing, conservation, closure, and provenance metadata.

Run it with:

```bash
fpm run --example fci_sol_gallery
```
