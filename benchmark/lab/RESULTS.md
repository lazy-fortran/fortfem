# Laboratory results

[V] These are initial synthetic MMS measurements, not a general ranking of
finite-element families or a VMEC-class equilibrium comparison. The assembled timings
were collected on faepmac2, GNU Fortran 16.1.0, release optimisation
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

## BEN1 curved torus MF consumer

[V] These runs use the MF1 core on faepmac1, GNU Fortran 16.2.0,
`-O3 -funroll-loops -fopenmp`, four OpenMP threads and one BLAS thread.
They are separate from the assembled faepmac2/GNU 16.1.0 measurements, so
compiler, geometry, host and algorithm differ. No isolated speedup is inferred.
The p=10, n=1 L1 point reaches 3.6082e-8 relative L2 error at 104200 DOFs in
1.2029 s. At ratio 1e6, unaligned L3 p=8 reaches 4.6640e-8 at 53920 DOFs in
9.1464 s, with 816 iterations and signed amplitude proxy 2.0131e-8.
High order resolves this case without requiring a field-aligned mesh, but the
iteration count still shows the cost of anisotropy.

[V] The regular O-grid is warped by the L1 map and includes the axis inside
Cartesian patch cells. The polar comparison instead identifies all axis angular
nodes at each phi and uses the assembled oracle with row ICHOL. On n=2 to n=4,
global/axis-ball rates are 2.934/3.055 for the O-grid and 2.993/3.444 for the
polar chart. There is no positive axis-ball rate loss on that refinement pair.
The coarsest pair has transient rates, so these short sequences do not prove
asymptotic chart independence. Polar axis constraints are not the structural
approach adopted by the laboratory.

[V] The independent quad-precision Cartesian finite-difference source check
agrees with generated L3 forcing within 256 double-precision epsilons times
its scale. The native problem test passed in 0.27 s. MF errors agree between
one and four threads to 1e-10, and all six FortSym products regenerate
byte-for-byte. All raw rows have their commands and binary fingerprints in
`mf-*-20261022.csv`. The initial polar rows have an unreported NNZ placeholder
of zero, not zero matrix storage; their RSS is measured normally.

[H] A genuine aligned periodic comparison still needs a nonmatching seam
constraint for the radially varying field-line shift. The CLI rejects that
request, because a periodic sinusoidal twist or an artificial Dirichlet cut
would change the benchmark. No aligned timing row is claimed.

[V] Process wall time includes assembly, constraints, preconditioner/factor setup,
solve and error integration. It excludes build and Python mesh preparation.
RSS is the child high-water mark from wait4/getrusage, in bytes in the CSV.
Every measured row includes its exact command, source revision, executable SHA256,
compiler and flags, host/job ID, tetra mesh checksum and linear residual.

| Problem | Criterion | Family | preconditioner | p | 1/h | Ratio | DOFs | rel L2 | axis rel L2 | chi_perp,num | iterations | wall (s) | RSS (MiB) |
|---|---:|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| L1 | 1 | gll-mf | pmg | 2 | 1 | 1.0 | 1000 | 2.5496e-02 | 8.9678e-03 | — | 6 | 0.046058 | 10.38 |
| L1 | 1 | gll-mf | pmg | 3 | 1 | 1.0 | 3120 | 2.7118e-03 | 4.6617e-04 | — | 6 | 0.047797 | 17.39 |
| L1 | 1 | gll-mf | pmg | 4 | 1 | 1.0 | 7120 | 9.7542e-04 | 1.2644e-04 | — | 7 | 0.073974 | 27.30 |
| L1 | 1 | gll-mf | pmg | 4 | 2 | 1.0 | 53920 | 2.1158e-05 | 1.4437e-06 | — | 9 | 0.43957 | 192.98 |
| L3 | 3 | gll-mf | pmg | 2 | 1 | 1.0 | 1000 | 5.5963e-03 | 1.9407e-04 | -2.5300e-04 | 7 | 0.027928 | 10.33 |
| L3 | 3 | gll-mf | pmg | 2 | 1 | 100.0 | 1000 | 6.3222e-03 | 2.2112e-03 | 2.5136e-03 | 24 | 0.050086 | 10.36 |
| L3 | 3 | gll-mf | pmg | 2 | 1 | 10000.0 | 1000 | 9.4873e-02 | 1.2899e-01 | 9.8073e-02 | 33 | 0.062388 | 10.31 |
| L3 | 3 | gll-mf | pmg | 2 | 1 | 1000000.0 | 1000 | 3.8342e-01 | 5.7130e-01 | 5.4549e-01 | 34 | 0.063866 | 10.39 |
| L5 | 2 | gll-mf | pmg | 2 | 1 | 1.0 | 1000 | 2.5496e-02 | 8.9678e-03 | — | 6 | 0.02641 | 10.34 |
| L5 | 2 | gll-polar | ichol_row | 2 | 1 | 1.0 | 1320 | 6.5497e-03 | 5.5125e-03 | — | 5 | 0.027761 | 14.38 |
| L5 | 2 | gll-mf | pmg | 2 | 2 | 1.0 | 7120 | 2.4427e-03 | 6.1290e-04 | — | 7 | 0.12856 | 57.86 |
| L5 | 2 | gll-polar | ichol_row | 2 | 2 | 1.0 | 10320 | 7.8235e-04 | 1.3926e-03 | — | 10 | 1.74 | 96.17 |
| L1 | 1 | gll-mf | pmg | 6 | 1 | 1.0 | 23160 | 3.6784e-05 | 1.3077e-06 | — | 9 | 0.22847 | 63.48 |
| L1 | 1 | gll-mf | pmg | 8 | 1 | 1.0 | 53920 | 1.1977e-06 | 1.4156e-08 | — | 11 | 0.53605 | 117.75 |
| L1 | 1 | gll-mf | pmg | 4 | 3 | 1.0 | 178800 | 3.6129e-06 | 1.8398e-07 | — | 10 | 1.9734 | 642.80 |
| L5 | 2 | gll-mf | pmg | 2 | 4 | 1.0 | 53920 | 3.1962e-04 | 7.7896e-05 | — | 8 | 6.3087 | 447.72 |
| L5 | 2 | gll-polar | ichol_row | 2 | 4 | 1.0 | 82080 | 9.8302e-05 | 1.2557e-04 | — | 20 | 144.13 | 755.84 |
| L1 | 1 | gll-mf | pmg | 10 | 1 | 1.0 | 104200 | 3.6082e-08 | 1.4539e-10 | — | 13 | 1.2029 | 208.36 |
| L3 | 3 | gll-mf | pmg | 6 | 1 | 1000000.0 | 23160 | 1.6577e-05 | 1.9204e-05 | 1.5136e-05 | 390 | 2.921 | 63.06 |
| L3 | 3 | gll-mf | pmg | 8 | 1 | 1000000.0 | 53920 | 4.6640e-08 | 1.8158e-08 | 2.0131e-08 | 816 | 9.1464 | 120.44 |

### Time to measured accuracy

[V] Minimum measured process wall time among qualified rows reaching the target.
A dash means no run reached it; no fitted or extrapolated timing is substituted.

| Problem | Family | Ratio | target rel L2 | DOFs | wall (s) |
|---|---|---:|---:|---:|---:|
| L1 | gll-mf | 1.0 | 1e-03 | 7120 | 0.073974 |
| L1 | gll-mf | 1.0 | 1e-05 | 53920 | 0.53605 |
| L1 | gll-mf | 1.0 | 1e-07 | 104200 | 1.2029 |
| L3 | gll-mf | 1.0 | 1e-03 | — | — |
| L3 | gll-mf | 1.0 | 1e-05 | — | — |
| L3 | gll-mf | 1.0 | 1e-07 | — | — |
| L3 | gll-mf | 100.0 | 1e-03 | — | — |
| L3 | gll-mf | 100.0 | 1e-05 | — | — |
| L3 | gll-mf | 100.0 | 1e-07 | — | — |
| L3 | gll-mf | 10000.0 | 1e-03 | — | — |
| L3 | gll-mf | 10000.0 | 1e-05 | — | — |
| L3 | gll-mf | 10000.0 | 1e-07 | — | — |
| L3 | gll-mf | 1000000.0 | 1e-03 | 23160 | 2.921 |
| L3 | gll-mf | 1000000.0 | 1e-05 | 53920 | 9.1464 |
| L3 | gll-mf | 1000000.0 | 1e-07 | 53920 | 9.1464 |
| L5 | gll-mf | 1.0 | 1e-03 | 53920 | 6.3087 |
| L5 | gll-mf | 1.0 | 1e-05 | — | — |
| L5 | gll-mf | 1.0 | 1e-07 | — | — |
| L5 | gll-polar | 1.0 | 1e-03 | 10320 | 1.74 |
| L5 | gll-polar | 1.0 | 1e-05 | — | — |
| L5 | gll-polar | 1.0 | 1e-07 | — | — |

### Axis refinement diagnostics

[V] Rates use the declared 1/n mesh parameter and physical quadrature norms.
These through-axis rates do not establish a rate loss relative to a polar chart.

| Problem | Family | p | refinement | global L2 rate | axis-tube L2 rate | global minus axis |
|---|---|---:|---|---:|---:|---:|
| L1 | gll-mf | 4 | 1 to 2 | 5.527 | 6.453 | -0.926 |
| L1 | gll-mf | 4 | 2 to 3 | 4.359 | 5.081 | -0.722 |
| L5 | gll-mf | 2 | 1 to 2 | 3.384 | 3.871 | -0.487 |
| L5 | gll-mf | 2 | 2 to 4 | 2.934 | 2.976 | -0.042 |
| L5 | gll-polar | 2 | 1 to 2 | 3.066 | 1.985 | 1.081 |
| L5 | gll-polar | 2 | 2 to 4 | 2.993 | 3.471 | -0.479 |

### Physical axis ball

[V] Radius 0.2 about X_axis(phi=0), with physical volume weights.
The reported norms use indicator quadrature, not certified integrals.

| Family | p | 1/h | global rel L2 | ball rel L2 |
|---|---:|---:|---:|---:|
| gll-mf | 2 | 1 | 2.54964e-02 | 7.78964e-04 |
| gll-mf | 2 | 2 | 2.44268e-03 | 4.48426e-04 |
| gll-mf | 2 | 4 | 3.19618e-04 | 5.39560e-05 |
| gll-polar | 2 | 1 | 6.54969e-03 | 3.81327e-03 |
| gll-polar | 2 | 2 | 7.82348e-04 | 1.23283e-03 |
| gll-polar | 2 | 4 | 9.83016e-05 | 1.13272e-04 |

| Family | p | refinement | global rate | ball rate | global minus ball |
|---|---:|---|---:|---:|---:|
| gll-mf | 2 | 1 to 2 | 3.384 | 0.797 | 2.587 |
| gll-mf | 2 | 2 to 4 | 2.934 | 3.055 | -0.121 |
| gll-polar | 2 | 1 to 2 | 3.066 | 1.629 | 1.436 |
| gll-polar | 2 | 2 to 4 | 2.993 | 3.444 | -0.452 |

## L6: reused N1 external reference measurements

[V] `l6-n1-reference-20261022.csv` imports thirty existing VMEC++/GVEC records
from N1's completed sCluster job 3355100 on node7. The adapter snapshot is
checked against `sweep-binaries.sha256`; no equilibrium cases or solver
harness are copied. Eighteen VMEC++ processes and four GVEC QA processes
completed; eight GVEC integer-family processes failed and retain blank
accuracy fields. The table below selects one-core points, not a best-code
ranking. Full requested CPU-affinity counts, input hashes, package versions,
DOF conventions and exact commands are retained in the CSV.

[V] Native solve time, adapter solve time and whole-process time have different
scopes. External wheel build compilers/flags were not exposed by N1 and remain
explicitly unavailable. N1 pressure errors use its own relative-pressure norm,
flux errors are maximum label defects, and sample positions belong to each
candidate. They therefore do not establish the common-physical-sample L6
contract defined in PROBLEMS.md. The QA case has no exact oracle, so its
accuracy is not filled with a force stopping tolerance. No lab-family versus
VMEC-class matched-accuracy speed claim is made.

| Code | Case | m | status / exit | DOFs | native solve (s) | adapter solve (s) | process (s) | peak RSS (MiB) | B sample error | p sample error (N1) | max flux-label defect |
|---|---|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|
| vmecpp | landreman-integer | 8 | sampled / 0 | 14355 | 0.2519 | 1.3704 | 7.29 | 469.62 | 0.0003222 | 0.0053373 | 0.0079158 |
| vmecpp | issan-integer | 8 | sampled / 0 | 28710 | 0.50646 | 1.6327 | 9.22 | 478.97 | 0.000687 | 0.0053031 | 0.010207 |
| gvec | landreman-integer | 4 | incomplete / 1 | 2178 | — | — | 15.79 | 240.93 | — | — | — |
| gvec | issan-integer | 4 | incomplete / 1 | 2178 | — | — | 10.83 | 243.54 | — | — | — |
| gvec | qa | 4 | sampled / 0 | 1089 | 17.087 | 18.882 | 19.44 | 243.41 | — | — | — |
