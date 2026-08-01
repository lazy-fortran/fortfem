---
title: Long-time mixed-wave invariants
---

# Long-time mixed-wave invariants

`test_mixed_wave_long_time` is a small, deterministic structure-preservation
campaign for the ideal mixed-wave Strang step.  It advances a non-identity,
non-commuting two-mode block for 2000 steps and evaluates the quadratic
Hamiltonian independently from the caller-owned mass matrices.  It then
applies the same number of signed backward steps and checks recovery of the
initial state.

This is a regression oracle for long-time energy drift and reversibility.  It
does not claim that a nonlinear or dissipative model is symplectic: damping,
resistivity, viscosity, and PML factors remain separate declared substeps.
Discrete-gradient/average-vector-field methods and broader problem-size
campaigns remain future extensions of the same port-Hamiltonian contract.
