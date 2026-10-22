# Discretisation laboratory results

<!-- MF1 BEGIN -->
## MF1: CPU GLL matrix-free diffusion and p-multigrid

[V] Synthetic data only. Correctness/convergence pilot: sCluster `node8`,
Slurm **3355289**, revision 2f12cc61657ff9846a39f9cbc5425157e44799eb.
GNU Fortran 12.2.0, CMake Release, `-O3 -Wall -Wextra -fcheck=bounds`,
OpenMP, and the repository's GCC <13 `-finline-matmul-limit=0` workaround.
`OPENBLAS_NUM_THREADS=1`; tests use four OpenMP threads. All five named
`test_lab_*` tests passed. Peak test-process RSS was 377,440 KiB.

[V] The updated orientation/seam test also passes on `faepmac1`, GNU Fortran
16.2.0, `-O2 -fcheck=all -Wall -Wextra`, four OpenMP threads. The resident
Gremlin gate covers all five named tests. Focused `fo lint` has no findings.
The pinned FortSym reproduction check `tools/codegen/check_generated.sh`
passes on `faepmac1`. Existing dependency/compiler warnings are separate:
there were no compiler diagnostics from the new laboratory sources.

### Operator and geometry checks

[V] GLL polynomial moments through degree 2p-1 satisfy an absolute 1e-15
threshold for p=1..16. Gauss/GLL interpolation, derivatives, physical line
scaling and periodic/nonperiodic dof maps pass their named test.

[V] Matrix-free apply agrees with independent, unfactored CSR quadrature to
1e-12 relative for p=1..6 on curved 2D and 3D patches, periodic in one direction.
The geometry degree is p+1. Exact diagonals and precompute/recompute agreement
pass. A separate physical-energy identity checks the generated inverse-metric
pullback, and singular geometry is rejected before division. Conforming-map
tests reverse the reference orientation on alternating cells and check shared
face dofs. The five-block disk and periodic torus have the expected seam dof
counts. These results inform criteria 2 and 4.

### MMS convergence (criteria 2 and 4)

[V] Each h level doubles both poloidal subdivisions and, in 3D, phi cells.
The annulus starts with 1 radial by 8 angular cells. The disk starts with one
cell per block, and the torus with five blocks by 8 phi cells. Geometry degree
is p+2; quadrature for error/load integration has p+3 points per direction.
The exact solution is exp(0.1x+0.2y) in 2D and exp(0.1x+0.2y+0.3z) in 3D,
with analytic forcing and nonzero Dirichlet lifting. The 3D boundary has nfp=5.

| Domain | p | Finest DOFs | L2 error | Final h rate | PCG iterations |
|---|---:|---:|---:|---:|---:|
| Annulus | 2 | 576 | 9.190e-6 | 3.00 | 9 |
| Annulus | 3 | 1,248 | 1.218e-7 | 4.00 | 9 |
| Full disk, O-grid | 2 | 337 | 6.262e-5 | 2.99 | 7 |
| Full disk, O-grid | 3 | 745 | 1.651e-6 | 4.00 | 7 |
| Full 3D torus, O-grid | 2 | 21,568 | 8.625e-4 | 2.97 | 14 |
| Full 3D torus, O-grid | 3 | 71,520 | 3.931e-5 | 3.93 | 14 |

[V] The final near-axis nodal maximum-error rates are 3.56/5.18 for p=2/3
on the disk and 3.93/4.42 on the torus. The sampled region is radius <0.2
from the axis. This is a nodal-error check, separate from the global L2 norm.
The axis is an ordinary mesh interior, with no polar boundary condition.

| p at fixed h=2 | Annulus L2 | Disk L2 | 3D torus L2 |
|---:|---:|---:|---:|
| 2 | 7.3307e-5 | 4.9924e-4 | 6.7795e-3 |
| 3 | 1.9469e-6 | 2.6378e-5 | 6.0006e-4 |
| 4 | 5.3366e-8 | 1.4344e-6 | 6.2585e-5 |
| 5 | 1.5475e-9 | 8.1904e-8 | 5.3923e-6 |
| 6 | 4.4429e-11 | 4.6348e-9 | 5.4738e-7 |

### p-multigrid and anisotropy (criteria 3 and 4)

[V] PCG stops at a freshly recomputed relative residual <=1e-10. The coarse
p=1 CSR operator is solved through FortFEM's native FortSparse LDLT hook.
The interpolation/restriction adjoint identity and V-cycle symmetry/positive
energy pass independent vector tests. The following iterations are for 2D
Poisson with a smooth sine-product MMS on h by h cells.

| p | h=2 iterations | h=4 iterations |
|---:|---:|---:|
| 2 | 3 | 3 |
| 3 | 5 | 5 |
| 4 | 6 | 7 |
| 5 | 6 | 6 |
| 6 | 7 | 7 |
| 7 | 7 | 7 |
| 8 | 7 | 7 |

[V] p=4, 4 by 4 cells, chi_perp=1. Alignment uses b=(1,0); the unaligned
case uses b=(1,1)/sqrt(2), including the mixed derivative in the MMS forcing.
The unaligned degradation is retained as evidence.

| chi_parallel/chi_perp | Alignment | Iterations | True residual | L2 error |
|---:|---|---:|---:|---:|
| 1 | Aligned | 7 | 1.0388e-12 | 2.6139e-6 |
| 1e3 | Aligned | 10 | 1.2115e-11 | 2.3486e-6 |
| 1e6 | Aligned | 10 | 1.4789e-12 | 2.3373e-6 |
| 1e6 | 45 degrees | 18 | 8.4180e-12 | 2.7027e-6 |

### Timing (criterion 4)

[V] Complete sweep: sCluster `node4`, Slurm **3355311**, revision
25ec09e24e4375579e67a9376f3c03b458b94de3, GNU Fortran 12.2.0,
the same Release flags as the pilot. All 128 apply cases and 16 solve
cases completed. There are 50 repeated applies per method and case.
The complete 1–16-thread series is in the external log; endpoints follow.

Each method cell reports **microseconds / million DOF/s / GFLOP/s**.
GFLOP/s is the contraction-work estimate; recompute has additional uncounted
metric arithmetic. The precompute/recompute MF error is zero in this sweep.

| p | DOFs | Threads | Precompute | Recompute | CSR | MF/CSR relative error | RSS MiB |
|---:|---:|---:|---|---|---|---:|---:|
| 1 | 27 | 1 | 8.93 / 3.022 / 0.831 | 16.62 / 1.624 / 0.447 | 0.58 / 46.614 / 1.184 | 1.67e-16 | 6.2 |
| 1 | 27 | 16 | 1538.13 / 0.018 / 0.005 | 21.46 / 1.258 / 0.346 | 4.96 / 5.443 / 0.138 | 1.67e-16 | 9.4 |
| 2 | 125 | 1 | 17.13 / 7.296 / 2.068 | 36.86 / 3.391 / 0.961 | 3.43 / 36.466 / 2.866 | 4.22e-16 | 6.5 |
| 2 | 125 | 16 | 6679.10 / 0.019 / 0.005 | 87.21 / 1.433 / 0.406 | 5.40 / 23.166 / 1.821 | 4.22e-16 | 12.9 |
| 3 | 343 | 1 | 41.12 / 8.340 / 2.639 | 76.17 / 4.503 / 1.425 | 18.75 / 18.295 / 3.178 | 2.11e-15 | 7.6 |
| 3 | 343 | 16 | 14519.60 / 0.024 / 0.007 | 17642.01 / 0.019 / 0.006 | 8.60 / 39.904 / 6.932 | 2.11e-15 | 49.8 |
| 4 | 729 | 1 | 89.76 / 8.122 / 2.897 | 162.59 / 4.484 / 1.599 | 89.62 / 8.134 / 2.625 | 6.72e-15 | 9.9 |
| 4 | 729 | 16 | 306.96 / 2.375 / 0.847 | 387.17 / 1.883 / 0.672 | 16.32 / 44.671 / 14.418 | 6.72e-15 | 32.5 |
| 5 | 1,331 | 1 | 193.83 / 6.867 / 2.746 | 318.69 / 4.176 / 1.670 | 252.64 / 5.268 / 2.833 | 1.38e-14 | 16.2 |
| 5 | 1,331 | 16 | 2460.09 / 0.541 / 0.216 | 5907.04 / 0.225 / 0.090 | 128.54 / 10.355 / 5.569 | 1.38e-14 | 44.1 |
| 6 | 2,197 | 1 | 358.80 / 6.123 / 2.723 | 560.53 / 3.919 / 1.743 | 575.93 / 3.815 / 3.169 | 2.37e-14 | 29.3 |
| 6 | 2,197 | 16 | 11474.55 / 0.191 / 0.085 | 893.71 / 2.458 / 1.093 | 70.86 / 31.004 / 25.760 | 2.37e-14 | 65.5 |
| 7 | 3,375 | 1 | 546.54 / 6.175 / 3.028 | 843.93 / 3.999 / 1.961 | 1309.53 / 2.577 / 3.128 | 4.17e-14 | 54.9 |
| 7 | 3,375 | 16 | 6523.73 / 0.517 / 0.254 | 1502.10 / 2.247 / 1.102 | 133.93 / 25.199 / 30.588 | 4.17e-14 | 102.6 |
| 8 | 4,913 | 1 | 884.06 / 5.557 / 2.982 | 1331.32 / 3.690 / 1.980 | 3079.26 / 1.596 / 2.711 | 6.83e-14 | 101.3 |
| 8 | 4,913 | 16 | 4101.94 / 1.198 / 0.643 | 1783.87 / 2.754 / 1.478 | 266.15 / 18.459 / 31.360 | 6.83e-14 | 148.8 |

[V] Full time to solution, including hierarchy setup and coarse LDLT:
3D p=4, 12 by 12 by 12 cells, **117,649 DOFs**, tolerance 1e-10.

| Threads | Iterations | Wall seconds | True relative residual | RSS MiB |
|---:|---:|---:|---:|---:|
| 1 | 7 | 3.536351 | 7.00e-12 | 185.8 |
| 2 | 7 | 2.618415 | 7.00e-12 | 187.3 |
| 3 | 7 | 2.225605 | 7.00e-12 | 188.4 |
| 4 | 7 | 2.043090 | 7.00e-12 | 189.8 |
| 5 | 7 | 1.947137 | 7.00e-12 | 191.1 |
| 6 | 7 | 2.048363 | 7.00e-12 | 192.4 |
| 7 | 7 | 1.845934 | 7.00e-12 | 193.7 |
| 8 | 7 | 1.786545 | 7.00e-12 | 194.8 |
| 9 | 7 | 1.842010 | 7.00e-12 | 195.9 |
| 10 | 7 | 1.889822 | 7.00e-12 | 197.7 |
| 11 | 7 | 1.877476 | 7.00e-12 | 198.9 |
| 12 | 7 | 1.750091 | 7.00e-12 | 200.2 |
| 13 | 7 | 1.976438 | 7.00e-12 | 201.4 |
| 14 | 7 | 1.741875 | 7.00e-12 | 203.2 |
| 15 | 7 | 2.299921 | 7.00e-12 | 203.1 |
| 16 | 7 | 2.051032 | 7.00e-12 | 205.1 |

### Reproduction and limits

Use the measured revision when reproducing numbers. On the remote host:

```bash
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 MPLBACKEND=Agg
unset DISPLAY
cmake -S . -B build/lab -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_Fortran_FLAGS='-O3 -Wall -Wextra -fcheck=bounds' \
  -DFORTFEM_BUILD_TESTING=ON -DFORTFEM_BUILD_BENCHMARKS=ON
cmake --build build/lab -j4 --target test_lab_gll_line test_lab_tensor_mesh \
  test_lab_tensor_apply test_lab_tensor_solver test_lab_tensor_mms fortfem_lab_mf1
ctest --test-dir build/lab -R '^test_lab_' --output-on-failure -V
for p in $(seq 1 8); do
  for t in $(seq 1 16); do
    OMP_NUM_THREADS=$t /usr/bin/time -v build/lab/fortfem_lab_mf1 "$p" 2 50 0
  done
done
for t in $(seq 1 16); do
  OMP_NUM_THREADS=$t /usr/bin/time -v build/lab/fortfem_lab_mf1 4 12 1 1
done
```

[V] Actual batch wrappers and full logs are outside Git at
`scluster:~/runs/fortfem-mf1-20261022/`: `mf1-pilot.sh`, `mf1-sweep.sh`,
`pilot-3355289.log` and `sweep-3355311.log`. Each wrapper records host,
Slurm job ID and revision. The memory guard reserves 16 cores and enforces
188.8 GiB; short sweep processes use a 0.1-second polling interval. RSS in the
tables is the full benchmark-process peak, including both operator variants
and CSR, rather than an isolated metric-storage measurement.

[V] The Mac kernel-generation environment additionally sets
`FORTFEM_CODEGEN_RUNNER=fpm`, `FPM_CXXFLAGS=-std=c++17`, `FPM_LDFLAGS=-lc++`,
Homebrew include/library paths, and SymEngine include/library paths from the
existing FortSym CMake dependency build. The FortSym source worktree is pinned
to 1725c8e16c0af9b7c6650cb82ec05aa6aec99aa4. The exact shell environment and
logs are retained on `faepmac1` under `~/wt/fortfem/`.

The apply sweep uses eight cells, so cell colouring offers little parallel
work at high thread counts. The large solve is the more representative CPU
parallelism check. The recompute variant rebuilds G from sampled Jacobians and
coefficients; it does not reinterpolate nodal geometry. The reported flop rate
counts contraction work and excludes the additional metric recomputation.
The solve error column is a true discrete residual, with rhs formed by applying
the discrete operator to a smooth nodal field. Continuum errors are the separate
MMS tables above. These measurements establish neither a VMEC/GVEC comparison
nor the FEM/BEM free-boundary coupling in criterion 1.

Mixed B-spline/Fourier operator adapters remain separate work. General automatic
coarsening assumes a regular cell complex, and propagation of explicit
periodic multigraph entities to coarse levels needs further validation.
The exact-diagonal validation covers cells without dof self-identifications;
one-cell periodic self-identifications are outside the verified scope.
<!-- MF1 END -->
