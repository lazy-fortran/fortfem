---
title: Linear perturbation block composition
---

# Linear perturbation block composition

`fortfem_linear_perturbation_composition` is the closure-neutral algebraic
boundary for a linear ideal or resistive perturbation client. The caller
supplies seven independently inspectable square blocks,

\[
  A(\omega)=L+P+V+W+S-\omega^2M+i\omega R,
\]

where `L` is the Lorentz block, `P` is the pressure/stress block, `M` is
inertia, `V` is the vacuum response, `W` is the wall response, `R` is the
resistive block, and `S` is a supplied singular-layer response. The frequency
may be complex. FortFEM does not derive any of these blocks, choose a pressure
closure, locate a rational surface, or read a plasma-code format.

The existing `assemble_linear_response_residual` and
`assemble_generalized_eigen_residual` consume the composed matrix for forced
and eigenvalue problems. `assemble_singular_layer_matching` remains the
separate rectangular inner/outer trace constraint; a caller may condense that
constraint into `S` or retain it in a larger complex block graph.

The composition exposes exact JVP and VJP actions for all seven matrices and
the complex frequency. The VJP uses the real part of the complex inner product,

\[
  \langle x,y\rangle_{\mathbb R}=\operatorname{Re}(x^H y),
\]

matching the complex residual and FEM/BEM/DtN/PML contracts elsewhere in the
library. Invalid shapes or nonfinite data produce a nonzero status and a zero
output rather than a partially assembled operator.

`linear_perturbation_metadata_t` records a positive normalization scale, a
normalization label, provenance, and explicit declarations for transpose
reciprocity and Hermitian nonnegative passivity. `UNDECLARED` values are valid:
the record describes a caller's convention and never asserts that an arbitrary
matrix satisfies it. Numerical diagnostics remain available through
`evaluate_linear_response_diagnostics`.

`test_linear_perturbation_composition` compares the value with an independent
elementwise manufactured oracle, differentiates every block and frequency by
a central difference, checks the real-complex adjoint identity, and verifies
metadata and invalid-shape rejection.
