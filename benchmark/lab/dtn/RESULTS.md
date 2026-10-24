# Dense exterior DtN build measurements

The exterior DtN build at fdf48aa7 was measured by sCluster job 3356679 using synthetic Fibonacci-sphere operators. The new implementation uses BLAS DGEMM for the Schur product and OpenMP symmetrisation. The old implementation is evaluated against the same operators in the same job.

| Trace size | Build, 1 thread (s) | Build, 16 threads (s) | Speedup |
|---|---:|---:|---:|
| 1536 | 0.2140 | 0.0390 | 5.49 |
| 3456 | 1.9700 | 0.2410 | 8.17 |
| 6144 | 10.6560 | 1.0950 | 9.73 |

Across the rows that also run the legacy oracle, the maximum relative matrix difference is 8.8743e-16 and the symmetry defect is zero. The six FEM+BEM consumer rows retain the previously recorded L2 errors: 3.379e-3, 2.454e-4, 1.127e-4 and 6.879e-5 for n=1,2,3,4 respectively. The n=2 and n=3 cases were run at both one and sixteen threads. Their DtN times fall from 0.2143/1.962 seconds at one thread to 0.04127/0.2603 seconds at sixteen threads.

Consumer provenance is FB1 commit 7d625ef7 with the DtN source file from fdf48aa7 staged into that build; the collected source hash matches the PR source exactly. This is mixed-tree consumer evidence, not a claim that the consumer was built from a single PR #119 commit. All original inputs, command lines and checksums remain outside Git. No new numerical run is needed to record this completed sweep.

The linked BLAS/LAPACK and OpenMP thread counts were both set to the tabulated count. Times are wall-clock build measurements from the same synthetic operators, rather than end-to-end discretisation timings. The synthetic regression also uses rectangular operators (37 flux unknowns and 23 trace unknowns), checks the flux equation against the original operators, and compares the resulting Schur matrix with an explicit-loop oracle.

The retained scaling log has SHA256 78d89c667089b3ef2df2cf22041a4c7b9d2623eb62f3bbda14e8b5b777aa4eaa. The staged consumer source has SHA256 2df4e79ea8936ae627b1da5134c7fd3890b61dac04c01dbb7449b3ce3adcb98e. Raw logs remain outside the repository. These historical measurements retain their original source identities and do not certify a subsequently restacked head.
