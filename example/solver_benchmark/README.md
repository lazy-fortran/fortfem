---
title: Solver benchmark
---

This example solves the unit-source Poisson system, shows the physical
solution first, and reports direct/PCG solve timings together with
memory-scalable row-oriented ILUT and ICHOL construction timings. The
ICHOL comparison uses a shifted symmetric copy of the assembled matrix so
the SPD requirement is explicit.

The first plot, `poisson_solution_2d.png`, is sampled from the finest PCG
solution (not a copied source field); `poisson_solution.csv` records it
against an independent unit-source square-Poisson series oracle. Timing and
factorization plots follow it.
