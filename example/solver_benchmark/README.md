---
title: Solver benchmark
---

This example solves a manufactured Poisson system, shows the physical
solution first, and reports direct/PCG solve timings together with
memory-scalable row-oriented ILUT and ICHOL construction timings. The
ICHOL comparison uses a shifted symmetric copy of the assembled matrix so
the SPD requirement is explicit.
