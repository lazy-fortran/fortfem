# Finite single-layer residual energy

`fortfem_bem_residual_energy` provides a diagnostic for an SPD single-layer
Galerkin matrix V on a fixed test space Y_h. It is independent of the kernel,
so the GS and toroidal-mode boundary operators can use the same routine.
Pass the explicitly symmetric matrix used by the solve: the existing GS
`v_solve` uses `(ops%v+transpose(ops%v))/2`.

Factor V once with `bem_residual_factor`. For the same immutable matrix and
its factor, call `bem_residual_energy(V, factor, b, lambda, work, eta, status)`.
The caller supplies work of length dim(Y_h); it returns the reference-space
correction. No allocation occurs in the repeated action. Matrix-vector
evaluation is OpenMP-parallel; LAPACK triangular solves remain serial with a
single-thread BLAS. Setup costs O(dim(Y_h)^3), and each action costs
O(dim(Y_h)^2). These are dense BEM costs, not an O(N) coupled error bound.

## Exact identity and claim limit

In exact arithmetic let V=L L^T, r=b-V lambda, and V lambda_h=b. The routine
computes

    eta = ||L^-1 r||_2 = ||lambda_h-lambda||_V,
    work = V^-1 r.

This identity is a finite-space statement. With exact continuous V and b,
lambda in Y_h, and e=lambda_exact-lambda, Galerkin orthogonality gives

    ||e||_V^2 = eta^2 + ||lambda_exact-lambda_h||_V^2.

Consequently eta is a **lower** bound on the continuous flux error. For
nested flux spaces on the same geometry and the same trace, it also equals
the square root of the increase in the symmetric Steklov quadratic form
when the flux space is enriched, provided the W term and all pairings are
unchanged. Geometry changes, trace interpolation and numerical quadrature
do not automatically satisfy these assumptions.

A counterexample is supplied by the unit-circle Laplace single layer in
normalized cosine modes: V_k=1/(2k), k>=1. A residual entirely above the
highest tested mode has zero measured norm and positive energy. The test
also uses lambda_k=rho^k and the independent infinite-tail identity
`sum rho^(2k)/(2k)=-log(1-rho^2)/2` to show increasing finite norms approach
the error from below. The reported quantity must stay a diagnostic in
ordinary floating arithmetic; no rounding enclosure is implemented.

## What remains for FEM-BEM certification

A continuous coupled upper bound must control both the Neumann residual
against arbitrary H1 traces and the single/double-layer residual outside
Y_h. A measured reference norm or a nonlinear stability lower bound alone
does not control those complements. No empirical saturation factor is
provided. A functional lifting of the full boundary residual is a possible
route; [Kurz et al. (2021)](https://doi.org/10.1007/s00211-021-01188-6)
give guaranteed functional estimates for Laplace BEM. Applying that route to
the weighted axisymmetric operator requires its own lifting, geometry,
quadrature and arithmetic bounds.

`fpm test test_bem_residual_energy --flag '-fopenmp'` checks exact recovery,
the Fourier tail, an analytically invertible non-diagonal SPD matrix and
input rejection. `fpm run --example bem_residual_study --profile release
--flag '-fopenmp'` measures setup, solve and repeated diagnostic cost on a
synthetic dense SPD fixture, with 1/2/4 threads and three repetitions.
It is a kernel cost study, not a GS error-efficiency measurement.

Measured 2026-10-15 on a shared Ryzen 9 5950X (32 logical CPUs), gfortran
16.2.1, FPM release `-fopenmp`, OpenBLAS limited to one thread. Each action
time is the median of three runs of 20 actions; the final column gives the
range over the three runs. Setup/factorization is excluded from the action
and included in the denominator. If a factor must be built solely for the
diagnostic, its setup cost must be added to the diagnostic cost as well.

| Flux DOF | Action 1 thread (s) | 2 threads (s) | 4 threads (s) | 4-thread action / (factor + solve) |
|---:|---:|---:|---:|---:|
| 512 | 1.77e-4 | 1.22e-4 | 1.01e-4 | 1.6--4.0% |
| 1024 | 2.77e-3 | 4.29e-4 | 3.59e-4 | 1.6--4.8% |
| 1536 | 9.74e-3 | 7.35e-3 | 4.37e-3 | 4.3--8.4% |
| 2048 | 2.71e-2 | 1.70e-2 | 1.43e-2 | 8.3--11.8% |

At 2048 DOF the paired 1-to-4-thread speedup is 1.48--2.39. The 1024-DOF
timings are particularly noisy (paired speedup 0.73--11.84); they do not
support a scaling claim. Cache transitions, host contention and the serial
triangular solves prevent interpreting these wall times as an asymptotic
complexity study. No full coupled-estimator cost or continuous-bound
efficiency was measured. The Fourier diagnostic/error ratio is
0.98634--0.999997 as the retained coarse modes increase from 2 to 8; this
approaches one from below, as the identity requires.
