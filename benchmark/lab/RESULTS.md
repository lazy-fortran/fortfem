# Laboratory results

[V] These are initial synthetic MMS measurements, not a general ranking of
finite-element families or a VMEC-class equilibrium comparison. All timings
below were collected on faepmac2, GNU Fortran 16.1.0, release optimisation
`-O3 -funroll-loops -fopenmp`, four requested OpenMP threads and one BLAS thread.
The raw rows retain the exact source revisions and executable SHA256s. The
problem-definition tests, PCG scale regression and existing advanced-solver
checks passed natively on that host. SymPy independently verified the exact
Poisson, curl-curl, closed-line tangency and coupled-model identities.

[V] The L1 tet P3 ladder reaches relative L2 error 5.3895e-6 at 103860 DOFs
in 328.19 s. Its measured 1e-3 target is reached by the smaller P3 point in
0.23618 s at 2940 DOFs. No assembled run reaches 1e-7. The spline ladder's
best measured L1 error is 4.6736e-5 (p=4, 2752 DOFs, 43.310 s).
All are single measurements, including setup and error integration.

[V] On the 32-by-32 P1 Guenter grid, increasing the ratio from 1 to 1e6
raises relative L2 error from 0.0027009 to 0.99913 and the amplitude-equivalent
numerical perpendicular conductivity from 0.00080345 to 1256.3. This is a
closed-line pollution failure, not a successful anisotropy-uniform baseline.
For the 3D sheared test the signed amplitude proxy is negative, and its
nonzero trace on the faceted boundary is included. That proxy is not a
physical negative conductivity or evidence of anisotropy-uniform convergence.

[V] The old `ichol` option was traced by a process sample to its dense reference
setup in `sparse_incomplete_cholesky.f90:124`. The measured ladder instead uses
`ichol_row`, which exposes the existing bounded row primitive with a 64-entry
row cap. The old alias remains available and unchanged. PCG's old absolute
p.Ap cutoff was independently reproduced under RHS scaling and repaired.
The lab's symmetric cell matrices are stored canonically for LDLT's exact
input-symmetry contract. No solver acceptance check was weakened.

[H] Norms are over each approximate physical mesh, with exact MMS trace data.
Tet geometry is piecewise affine; spline geometry is the Greville control-net
approximation documented in PROBLEMS.md. The rows do not bound the error of
the exact curved boundary or establish a BEM/DtN coupling. L4 is a constant
coupled elliptic model, not a nonlinear equilibrium solve. L5's assembled
axis diagnostic is a tube, not a ball or a polar-chart comparison.

## Reproduction and pilot

[V] Raw ledgers: [baselines-20261022.csv](baselines-20261022.csv) and
[additional-20261022.csv](additional-20261022.csv). Commands are in every row.
Build and execution recipes are in [README.md](README.md); regenerate the
following tables with:

```sh
python3 benchmark/lab/summarize.py benchmark/lab/baselines-20261022.csv \
  benchmark/lab/additional-20261022.csv --output assembled-tables.md
```

[V] One sCluster pilot was submitted as job **3355332**, four cores, one hour,
at ad40b17a2ceab358a41e4491b49704faa8abe44e. It was not waited on. Its worktree is
`/home/ert/wt/fortfem/ben1-20261022`, log `ben1-pilot-3355332.log`, and outputs
`/home/ert/ben1-runs/3355332`. It predates the lab's canonical-symmetry fix,
so it is not verification of the final head. The template guards 47.2 GiB,
which is the enforced four-core memory share. Strong/weak scaling templates
are prepared but no scaling campaign was launched before pilot qualification.

[V] KIN6D's N1 cases and adapter are reused by L6.py. Its oracle work is in
[KIN6D PR 16](https://github.com/itpplasma/kin6d/pull/16), already merged.
No external equilibrium accuracy or speed row is claimed here. The matched
B, p, iota and flux contract, including zero-pressure normalization and input
conventions, is in PROBLEMS.md.

## MF core

[H] The matrix-free path exists on `origin/lab/mf1` (PR 110). Its lab adapter
and rows are maintained on `lab/ben1-mf`, stacked on that branch. This section
is separate from the assembled measurements. An aligned periodic sheared-mesh
adapter still needs a nonmatching seam constraint; it is rejected rather than
relabelled as an aligned run.

## Assembled baselines

[V] Process wall time includes assembly, constraints, preconditioner/factor setup,
solve and error integration. It excludes build and Python mesh preparation.
RSS is the child high-water mark from wait4/getrusage, in bytes in the CSV.
Every measured row includes its exact command, source revision, executable SHA256,
compiler and flags, host/job ID, tetra mesh checksum and linear residual.

| Problem | Criterion | Family | preconditioner | p | 1/h | Ratio | DOFs | rel L2 | axis rel L2 | chi_perp,num | iterations | wall (s) | RSS (MiB) |
|---|---:|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| L1 | 1 | tet | ichol_row | 1 | 1 | 1.0 | 180 | 7.2896e-02 | 6.6446e-02 | — | 3 | 0.01531 | 3.97 |
| L1 | 1 | tet | ichol_row | 1 | 2 | 1.0 | 1320 | 2.2431e-02 | 1.9535e-02 | — | 5 | 0.12013 | 7.56 |
| L1 | 1 | tet | ichol_row | 2 | 1 | 1.0 | 1000 | 5.6289e-03 | 8.3081e-03 | — | 4 | 0.033948 | 5.94 |
| L1 | 1 | tet | ichol_row | 2 | 2 | 1.0 | 9040 | 7.6971e-04 | 1.3425e-03 | — | 11 | 2.3341 | 35.02 |
| L1 | 1 | tet | ichol_row | 3 | 1 | 1.0 | 2940 | 2.9291e-04 | 4.6786e-04 | — | 8 | 0.23618 | 13.16 |
| L1 | 1 | tet | ichol_row | 3 | 2 | 1.0 | 28920 | 2.6146e-05 | 1.8945e-05 | — | 19 | 24.564 | 110.31 |
| L3a | 3 | tet | ldlt | 1 | 8 | 1.0 | 81 | 4.2266e-02 | — | 1.2917e-02 | 0 | 0.0072027 | 3.59 |
| L3a | 3 | tet | ldlt | 1 | 8 | 100.0 | 81 | 6.5259e-01 | — | 2.0309e+00 | 0 | 0.0071706 | 3.59 |
| L3a | 3 | tet | ldlt | 1 | 8 | 10000.0 | 81 | 9.9460e-01 | — | 1.9701e+02 | 0 | 0.007021 | 3.53 |
| L3a | 3 | tet | ldlt | 1 | 8 | 1000000.0 | 81 | 9.9995e-01 | — | 1.9694e+04 | 0 | 0.006871 | 3.59 |
| L3a | 3 | tet | ldlt | 1 | 16 | 1.0 | 289 | 1.0755e-02 | — | 3.2169e-03 | 0 | 0.0093093 | 3.91 |
| L3a | 3 | tet | ldlt | 1 | 16 | 100.0 | 289 | 3.2146e-01 | — | 5.5204e-01 | 0 | 0.0096692 | 3.89 |
| L3a | 3 | tet | ldlt | 1 | 16 | 10000.0 | 289 | 9.7865e-01 | — | 5.0117e+01 | 0 | 0.009629 | 3.91 |
| L3a | 3 | tet | ldlt | 1 | 16 | 1000000.0 | 289 | 9.9978e-01 | — | 5.0034e+03 | 0 | 0.0091691 | 3.80 |
| L3a | 3 | tet | ldlt | 1 | 32 | 1.0 | 1089 | 2.7009e-03 | — | 8.0345e-04 | 0 | 0.026738 | 4.98 |
| L3a | 3 | tet | ldlt | 1 | 32 | 100.0 | 1089 | 1.0762e-01 | — | 1.5468e-01 | 0 | 0.027685 | 5.02 |
| L3a | 3 | tet | ldlt | 1 | 32 | 10000.0 | 1089 | 9.1988e-01 | — | 1.2644e+01 | 0 | 0.028047 | 5.00 |
| L3a | 3 | tet | ldlt | 1 | 32 | 1000000.0 | 1089 | 9.9913e-01 | — | 1.2563e+03 | 0 | 0.027436 | 4.97 |
| L3a | 3 | tet | ldlt | 1 | 64 | 1.0 | 4225 | 6.7598e-04 | — | 2.0081e-04 | 0 | 0.30954 | 10.61 |
| L3a | 3 | tet | ldlt | 1 | 64 | 100.0 | 4225 | 2.9644e-02 | — | 4.4180e-02 | 0 | 0.30953 | 10.50 |
| L3a | 3 | tet | ldlt | 1 | 64 | 10000.0 | 4225 | 7.4261e-01 | — | 3.2203e+00 | 0 | 0.30403 | 10.48 |
| L3a | 3 | tet | ldlt | 1 | 64 | 1000000.0 | 4225 | 9.9653e-01 | — | 3.1458e+02 | 0 | 0.30369 | 10.55 |
| L3 | 3 | tet | ldlt | 2 | 1 | 1.0 | 1000 | 9.7175e-03 | 6.2204e-03 | -2.4033e-04 | 0 | 0.041314 | 6.12 |
| L3 | 3 | tet | ldlt | 2 | 1 | 100.0 | 1000 | 1.7129e-02 | 2.4104e-02 | -7.8616e-03 | 0 | 0.042595 | 6.11 |
| L3 | 3 | tet | ldlt | 2 | 1 | 10000.0 | 1000 | 4.2370e-02 | 7.7100e-02 | -3.0951e-02 | 0 | 0.04226 | 6.11 |
| L3 | 3 | tet | ldlt | 2 | 1 | 1000000.0 | 1000 | 4.3320e-02 | 7.8915e-02 | -3.1768e-02 | 0 | 0.043184 | 6.11 |
| L1 | 1 | tet | ichol_row | 3 | 3 | 1.0 | 103860 | 5.3895e-06 | 3.4217e-06 | — | 29 | 328.19 | 355.55 |
| L1 | 1 | bspline | ichol_row | 2 | 1 | 1.0 | 336 | 9.1669e-03 | 8.9444e-03 | — | 3 | 0.42605 | 6.50 |
| L1 | 1 | bspline | ichol_row | 2 | 2 | 1.0 | 1476 | 1.0783e-03 | 1.7824e-03 | — | 6 | 1.7076 | 20.62 |
| L1 | 1 | bspline | ichol_row | 3 | 1 | 1.0 | 550 | 2.8723e-03 | 2.3846e-03 | — | 4 | 0.9158 | 15.66 |
| L1 | 1 | bspline | ichol_row | 3 | 2 | 1.0 | 2058 | 1.7377e-04 | 1.4281e-04 | — | 9 | 8.9964 | 90.94 |
| L1 | 1 | bspline | ichol_row | 4 | 1 | 1.0 | 828 | 2.8115e-03 | 2.1392e-03 | — | 9 | 5.1828 | 49.67 |
| L1 | 1 | bspline | ichol_row | 4 | 2 | 1.0 | 2752 | 4.6736e-05 | 3.2027e-05 | — | 492 | 43.31 | 298.12 |
| L2 | 1 | tet | ichol_row | 1 | 1 | 1.0 | 820 | 3.1298e-01 | — | — | 4 | 0.32969 | 4.91 |
| L4 | 4 | tet | ichol_row | 2 | 1 | 1.0 | 3000 | 5.6044e-03 | 8.2743e-03 | — | 8 | 0.072784 | 5.92 |
| L4 | 4 | tet | ichol_row | 2 | 2 | 1.0 | 27120 | 7.6652e-04 | 1.3336e-03 | — | 21 | 4.5597 | 35.03 |
| L5 | 2 | tet | ichol_row | 2 | 1 | 1.0 | 1000 | 5.6289e-03 | 8.3081e-03 | — | 4 | 0.033394 | 5.98 |
| L5 | 2 | tet | ichol_row | 2 | 2 | 1.0 | 9040 | 7.6971e-04 | 1.3425e-03 | — | 11 | 2.3325 | 35.11 |

### Time to measured accuracy

[V] Minimum measured process wall time among qualified rows reaching the target.
A dash means no run reached it; no fitted or extrapolated timing is substituted.

| Problem | Family | Ratio | target rel L2 | DOFs | wall (s) |
|---|---|---:|---:|---:|---:|
| L1 | bspline | 1.0 | 1e-03 | 2058 | 8.9964 |
| L1 | bspline | 1.0 | 1e-05 | — | — |
| L1 | bspline | 1.0 | 1e-07 | — | — |
| L1 | tet | 1.0 | 1e-03 | 2940 | 0.23618 |
| L1 | tet | 1.0 | 1e-05 | 103860 | 328.19 |
| L1 | tet | 1.0 | 1e-07 | — | — |
| L2 | tet | 1.0 | 1e-03 | — | — |
| L2 | tet | 1.0 | 1e-05 | — | — |
| L2 | tet | 1.0 | 1e-07 | — | — |
| L3 | tet | 1.0 | 1e-03 | — | — |
| L3 | tet | 1.0 | 1e-05 | — | — |
| L3 | tet | 1.0 | 1e-07 | — | — |
| L3 | tet | 100.0 | 1e-03 | — | — |
| L3 | tet | 100.0 | 1e-05 | — | — |
| L3 | tet | 100.0 | 1e-07 | — | — |
| L3 | tet | 10000.0 | 1e-03 | — | — |
| L3 | tet | 10000.0 | 1e-05 | — | — |
| L3 | tet | 10000.0 | 1e-07 | — | — |
| L3 | tet | 1000000.0 | 1e-03 | — | — |
| L3 | tet | 1000000.0 | 1e-05 | — | — |
| L3 | tet | 1000000.0 | 1e-07 | — | — |
| L3a | tet | 1.0 | 1e-03 | 4225 | 0.30954 |
| L3a | tet | 1.0 | 1e-05 | — | — |
| L3a | tet | 1.0 | 1e-07 | — | — |
| L3a | tet | 100.0 | 1e-03 | — | — |
| L3a | tet | 100.0 | 1e-05 | — | — |
| L3a | tet | 100.0 | 1e-07 | — | — |
| L3a | tet | 10000.0 | 1e-03 | — | — |
| L3a | tet | 10000.0 | 1e-05 | — | — |
| L3a | tet | 10000.0 | 1e-07 | — | — |
| L3a | tet | 1000000.0 | 1e-03 | — | — |
| L3a | tet | 1000000.0 | 1e-05 | — | — |
| L3a | tet | 1000000.0 | 1e-07 | — | — |
| L4 | tet | 1.0 | 1e-03 | 27120 | 4.5597 |
| L4 | tet | 1.0 | 1e-05 | — | — |
| L4 | tet | 1.0 | 1e-07 | — | — |
| L5 | tet | 1.0 | 1e-03 | 9040 | 2.3325 |
| L5 | tet | 1.0 | 1e-05 | — | — |
| L5 | tet | 1.0 | 1e-07 | — | — |

### Axis refinement diagnostics

[V] Rates use the declared 1/n mesh parameter and physical quadrature norms.
These through-axis rates do not establish a rate loss relative to a polar chart.

| Problem | Family | p | refinement | global L2 rate | axis-tube L2 rate | global minus axis |
|---|---|---:|---|---:|---:|---:|
| L1 | bspline | 2 | 1 to 2 | 3.088 | 2.327 | 0.760 |
| L1 | bspline | 3 | 1 to 2 | 4.047 | 4.062 | -0.015 |
| L1 | bspline | 4 | 1 to 2 | 5.911 | 6.062 | -0.151 |
| L1 | tet | 1 | 1 to 2 | 1.700 | 1.766 | -0.066 |
| L1 | tet | 2 | 1 to 2 | 2.870 | 2.630 | 0.241 |
| L1 | tet | 3 | 1 to 2 | 3.486 | 4.626 | -1.140 |
| L1 | tet | 3 | 2 to 3 | 3.895 | 4.221 | -0.326 |
| L4 | tet | 2 | 1 to 2 | 2.870 | 2.633 | 0.237 |
| L5 | tet | 2 | 1 to 2 | 2.870 | 2.630 | 0.241 |
