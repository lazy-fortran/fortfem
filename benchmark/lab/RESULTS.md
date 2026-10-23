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

Mixed B-spline/Fourier operator adapters are in the IGA-AX section. General automatic
coarsening assumes a regular cell complex, and propagation of explicit
periodic multigraph entities to coarse levels needs further validation.
The exact-diagonal validation covers cells without dof self-identifications;
one-cell periodic self-identifications are outside the verified scope.
<!-- MF1 END -->

<!-- IGAAX BEGIN -->
## IGA-AX: GLL, B-spline and Fourier lines through one MF operator

All rows use the same driver, operator, generated metric kernel, MG-PCG and
coarse Cholesky. Only the line tables and the transfer rule change between
families. Solves use relative tolerance 1e-10 from a zero initial guess, with
Chebyshev-Jacobi smoothing of degree 3. Wall times include hierarchy setup and
coarse factorisation but not error evaluation. RSS is the process VmHWM.

[V] Evidence run: sCluster `node7` (exclusive), Slurm **3356269**, revision
3bf52b951e475b9793c480d3de35aa19d929c1ec (pre-rebase; its source tree equals the
source of the commit adding this section, which only adds text), GNU Fortran 12.2.0, fpm 0.13.0
`--profile release` (`-O3 -funroll-loops -fcoarray=single -fPIC -fopenmp`),
`OMP_PROC_BIND=close OMP_PLACES=cores`. Driver SHA256 388f41ff…, libc
6b4a4535…. The CSVs `disk.csv`, `iters.csv`, `scan.csv` and `threads.csv` (with
SHA256SUMS) are in `scluster:~/runs/fortfem-igaax-20261023/sweep-3356269/`.
The pilot 3356244 ran on a shared node, so its timings are not used.

### Operator checks (criterion 4)

[V] `test_lab_line_operator` checks the spline MF apply against the existing
assembled B-spline H1 operators for 2D/3D and p=1..4 (worst relative error
2.9e-14), mixed spline x GLL x Fourier against the CSR oracle (3.2e-16) and
spline x Fourier exactness on trigonometric polynomials (L2 7.6e-15). It also
checks Galerkin nestedness P^T A_f P = A_c across O-grid interfaces and
Fourier mode truncation (2.1e-14). `test_lab_line_hierarchy` checks rates,
bounded iterations and the GLL p-prolongation.

### Disk MMS through the axis vs C1 polar splines (criterion 2)

u = exp(0.1x + 0.2y) on the unit disk. The O-grid has five C0-glued patches with
n spans or cells per patch direction, and the C1 polar comparison uses n radial
spans, 4n angular spans and a dense direct solve. "Axis" columns restrict to
r < 0.2. The rows show n=16 at 4 threads, and the rates are log2 of the error
ratio between n=8 and n=16.

| Family | p | DOFs | Its | Global rel L2 | Axis max | Axis L2 | Global rate |
|---|---:|---:|---:|---:|---:|---:|---:|
| spline O-grid | 2 | 1,480 | 6 | 1.31e-7 | 4.25e-9 | 8.83e-10 | 3.12 |
| GLL O-grid | 2 | 5,185 | 7 | 1.27e-7 | 4.28e-9 | 8.83e-10 | 2.98 |
| C1 polar | 2 | 1,027 | direct | 1.42e-7 | 2.94e-8 | 3.80e-9 | 3.07 |
| spline O-grid | 3 | 1,657 | 14 | 7.64e-9 | 1.01e-10 | 2.98e-11 | 4.02 |
| GLL O-grid | 3 | 11,617 | 7 | 1.46e-9 | 1.09e-11 | 3.13e-12 | 3.90 |
| C1 polar | 3 | 1,091 | direct | 6.15e-9 | 8.36e-10 | 8.43e-11 | 4.07 |
| spline O-grid | 4 | 1,844 | 27 | 2.31e-10 | 9.57e-11 | 2.12e-11 | 5.03 |
| GLL O-grid | 4 | 20,609 | 9 | 1.66e-11 | 8.27e-12 | 2.80e-12 | 4.81 |
| C1 polar | 4 | 1,155 | direct | 2.11e-10 | 2.66e-11 | 2.70e-12 | 5.11 |

[V] All three discretisations converge globally at rate p+1. The O-grid has no
axis degradation. At p=2 its near-axis maximum error is 7x below the polar
spline at the same n. At p>=3 the O-grid axis columns stop at about 1e-10,
which is the PCG tolerance. A local rerun with `tol=1e-14` (revision
3bf52b95, 4 threads) gives an axis maximum of 6.4e-12 at n=16 and 3.5e-13 at
n=32 for spline p=3, 130x below polar at n=16, and 9.7e-14 for GLL p=3 at
n=32. Per DOF, the polar spline space is the most compact, but it needs the
extraction operator and a C1 pole constraint, and it has no multigrid here.

### Multigrid iterations vs p, h and anisotropy (criterion 4)

Disk MMS, 4 threads, PCG iterations to 1e-10 for n = 4, 8, 16, 32, 64
(spline h-MG by span halving) or n = 2, 4, 8, 16 (GLL p-MG on fixed cells):

| Hierarchy | p=2 | p=3 | p=4 | p=5 | p=6 | p=8 |
|---|---|---|---|---|---|---|
| spline h-MG, max continuity | 6,6,6,6,6 | 15,14,14,14,14 | 27 (all n) | 61,49,54,54,53 | | |
| GLL p-MG | 6,7,7,7 | | 7,7,8,9 | | 8,8,10,10 | 9,10,11,12 |

[V] Spline h-multigrid is h-independent but not p-robust, while GLL
p-multigrid is close to robust in both p and h. [H] The cause is the
point-Jacobi Chebyshev smoother acting on maximally smooth splines. Hofreither
and Takacs report p-robust smoothers built on stable splittings of the spline
space (SIAM J. Numer. Anal. 55(4):2004-2024, 2017; arXiv:1607.05035), which
are not implemented here.

Aligned anisotropy uses u = sin(pi x) sin(pi y) on one spline box with
K = diag(chi, 1), at n = 16 and n = 64 spans:

| chi | Coarsening | p=2 | p=3 | p=4 |
|---|---|---|---|---|
| 1e3 | full | 17, 49 | 24, 64 | 29, 89 |
| 1e3 | semi (strong direction only) | 6, 6 | 11, 9 | 16, 17 |
| 1e6 | full | 17, 95 | 25, 121 | 31, 174 |
| 1e6 | semi | 6, 6 | 10, 8 | 16, 15 |

[V] Semi-coarsening along the aligned strong direction makes the iteration
counts independent of chi and h, while full coarsening degrades with both.

### L1 torus at matched error (criterion 4)

L1 is the shaped nfp=5 torus from `benchmark/lab/PROBLEMS.md`, u = sin x +
cos y + sin z + (x+y+z)/10, solved on the five-patch O-grid through the axis
and extended in phi. Phi is discretised with periodic lines of the same family
(12n cells or spans) or with M Fourier modes. A 16-thread scan covered 69 runs
(`scan.csv`). For each family the cheapest run at rel L2 <= 1e-6 and <= 1e-8
was then timed at 1, 4 and 16 threads:

| Target | Family | p | n | phi | DOFs | Its | rel L2 | Wall 1 / 4 / 16 thr (s) | RSS 16 thr MiB |
|---|---|---:|---:|---|---:|---:|---:|---|---:|
| 1e-6 | GLL | 8 | 1 | 12 cells | 32,352 | 14 | 6.74e-7 | 0.49 / 0.23 / **0.23** | 262 |
| 1e-6 | GLL x Fourier | 8 | 1 | M=24 | 16,513 | 20 | 5.72e-7 | 0.64 / **0.52** / 0.68 | 485 |
| 1e-6 | spline x Fourier | 4 | 8 | M=24 | 30,772 | 31 | 2.92e-7 | 11.7 / 4.13 / 2.63 | 1,390 |
| 1e-6 | spline | 4 | 8 | 96 spans | 60,288 | 40 | 4.02e-7 | 67.2 / 19.5 / 7.75 | 1,776 |
| 1e-8 | GLL | 8 | 2 | 24 cells | 252,096 | 15 | 1.84e-9 | 4.10 / 1.81 / **1.26** | 972 |
| 1e-8 | GLL x Fourier | 6 | 3 | M=32 | 107,705 | 31 | 9.06e-9 | 6.66 / 2.94 / 2.51 | 5,505 |
| 1e-8 | spline x Fourier | 4 | 16 | M=32 | 119,860 | 30 | 7.82e-9 | 71.8 / 22.8 / 10.6 | 6,296 |
| 1e-8 | spline (closest, 16 thr) | 4 | 16 | 192 spans | 354,048 | 32 | 1.06e-8 | – / – / 49.1 | 11,889 |

[V] At matched error the tensor GLL discretisation is the fastest by a wide
margin. At 1e-8 and 16 threads it is 8x faster than spline x Fourier, 2x
faster than GLL x Fourier and 39x faster than the pure spline closest to the
target. Fourier in phi gives the fewest DOFs per accuracy, but its global phi
cell multiplies the quadrature and geometry storage, and its mode-truncation
multigrid needs 20-31 iterations instead of 14-15. Pure splines need no fewer
DOFs than GLL p=8 here (60,288 against 32,352 near 5e-7), and their cost per DOF
is far higher. With p+1 Gauss points per span and direction, each span has
(p+1)^3 quadrature points for about one new DOF. The spline h-MG iteration counts also grow with
p. [H] Weighted quadrature (Sangalli and Tani, arXiv:1712.08565) and p-robust
spline smoothers are the known remedies, and neither is implemented here.

[V] Thread scaling is poor for GLL. The 1e-8 run has a speedup of 3.2 at 16
threads, because setup (0.88 s of 1.26 s) does not scale. The compute-bound
pure spline reaches 8.7 and spline x Fourier 4.5 to 6.8. The small GLL x
Fourier run slows down beyond 4 threads.

### Verdict and limits

[V] For criterion 4 on these problems, use high-order isoparametric GLL with
p-multigrid as the default discretisation. Splines and Fourier remain
available through the same operator. The axis-regular O-grid satisfies
criterion 2 for both GLL and splines without any polar treatment, and the C1
polar spline remains a comparison only. [H] A spline route would need
weighted quadrature and p-robust smoothing to compete. Fourier in phi pays off
only if its setup and storage are reduced, for example with a phi rule and
geometry per level or with FFT-based phi application. Neither has been tried
here.

Not measured here: VMEC++/GVEC comparisons, FEM/BEM coupling (criterion 1),
non-nested/island cases, G1 patch interfaces (the O-grid is C0 across patches)
and NURBS weights. All runs were single-node, with pointwise kernels from the
FortSym-generated metric module.

### Reproduction

```bash
# sCluster, exclusive node; script and CSVs in ~/runs/fortfem-igaax-20261023/
# source tree of 3bf52b95 = the commit that adds this section
OMP_NUM_THREADS=16 fpm build --profile release
B=$(find build -path '*app/fortfem_lab_igaax' | head -1)
export OPENBLAS_NUM_THREADS=1 OMP_PROC_BIND=close OMP_PLACES=cores
OMP_NUM_THREADS=4 $B problem=disk family=spline p=3 n=16 header=1
OMP_NUM_THREADS=4 $B problem=disk family=polar p=3 n=16 nphi=64
OMP_NUM_THREADS=4 $B problem=aniso family=spline p=3 n=64 chi=1e6 coarsen=semi
OMP_NUM_THREADS=16 $B problem=L1 family=gll p=8 n=2 nphi=24
OMP_NUM_THREADS=16 $B problem=L1 family=spline-fourier p=4 n=16 m=32
OMP_NUM_THREADS=16 $B problem=L1 family=gll-fourier p=6 n=3 m=32
```

The full grid is in `igaax-sweep.sh` next to the CSVs. Each row in the CSVs
records threads, DOFs, levels, iterations, the final residual, setup, solve
and total time, global and near-axis errors and RSS.
<!-- IGAAX END -->

<!-- UNS1 BEGIN -->
## UNS1: curved unstructured quad/hex input and island-adapted meshes

[V] Synthetic data only. Gmsh fixtures in `test/data/gmsh/` are regenerated by
`make_fixtures.py` (Gmsh 4.15.2). `test_lab_unstructured_mesh` and the
`lab_uns1_island_*` CTest cases pass on `faepmac1` (GNU Fortran 16.2.0,
`-O2 -fcheck=all -Wall -Wextra`, four threads; `fpm test
test_lab_unstructured_mesh` also passes) and in sCluster pilot **3356132**
(`node10`, revision 3c2fe688, GNU Fortran 12.2.0, Release
`-O3 -Wall -Wextra -fcheck=bounds`). Focused `fo lint` reports no reliance on
short-circuit evaluation and no unused imports.

### Mesh checks (criteria 2 and 4)

[V] Matrix-free apply agrees with the assembled CSR to 1e-12 on every mesh of
the test: Gmsh disks of geometry order 1..3, an order-2 hex cylinder, order-2/3
periodic hex boxes, a mirrored (reoriented) disk, the island-chain quad mesh and
its periodic quad x phi extrusion. Each mesh has positive Jacobians on a
(q+2)-point GLL lattice and reproduces its shared-node coordinates. The
reader flags the three cells that Gmsh's own minimum scaled Jacobian marks as
inverted in the unoptimised order-3 disk. The island mesh has Euler
characteristic 1, every X-point is a vertex of valence 8, and the separatrix
vertices lie on the level set to 1e-12.

### Island chain, anisotropic diffusion (criteria 3 and 4)

Problem: single-helicity 4/3 island of Hölzl et al., Phys. Plasmas 14, 052501
(2007), doi:10.1063/1.2722980, Sec. II.A, in the helical-frame plane.
[V] The q profile, perturbation form, r_s = 0.502a and w = 0.068a were checked
against the paper body. [H] R0 = 3.3a is inferred from the paper's w/w_c
values, as the paper does not state it. The manufactured T is constant on
field lines plus a perturbation (`island_diffusion check` verifies its
gradients). `ogrid` is the unaligned MF1 five-block O-grid through the axis,
and `island` is the island-adapted mesh. `rel_l2` is the relative L2 error,
`interp` that of the nodal interpolant (best-approximation proxy) and `spur`
is (E_total - E_perp)/E_perp of T_h, the spurious parallel energy.

[V] faepmac1, four threads, p = 4, p-MG with LDLT coarse (tol 1e-10, 500
iterations at most). The chi = 1e9 rows use `direct` (LDLT of the same assembled
fine operator), since p-MG does not converge there. Where p-MG converges,
direct and p-MG agree in every printed digit.

| mesh | DOFs | interp | rel_l2 1e3 | rel_l2 1e6 | rel_l2 1e9 | spur 1e6 | spur 1e9 | its 1e3/1e6 |
|---|---|---|---|---|---|---|---|---|
| ogrid L8 | 5185 | 8.03e-3 | 7.16e-3 | 9.24e-3 | 3.83e-2 | 1.41e-2 | 1.79e-1 | 8/181 |
| ogrid L16 | 20609 | 1.83e-3 | 1.75e-3 | 2.37e-3 | 1.74e-2 | 4.30e-3 | 2.46e-2 | 9/248 |
| island L2 | 4641 | 6.76e-3 | 9.10e-3 | 8.32e-3 | 1.25e-2 | 1.57e-3 | 1.62e-2 | 27/113 |
| island L4 | 16449 | 1.98e-3 | 1.96e-3 | 1.93e-3 | 4.82e-3 | 6.91e-4 | 1.32e-2 | 34/186 |

[V] At chi = 1e3 both meshes reach their interpolation error at similar DOFs.
At 1e9 the unaligned O-grid error is 9.5 times its interpolation error
(L16), against 2.4 times for the island mesh (L4), whose error is 3.6 times
lower at 0.8 times the DOFs. Neither mesh is pollution-free at 1e9.
[V] p-MG with Jacobi/Chebyshev smoothing is not robust in chi: 9 -> 248 -> >500
iterations (ogrid L16) and 34 -> 186 -> >500 (island L4) for chi = 1e3, 1e6,
1e9. The SA-AMG coarse solve (FortSparse #29) does not help: 150 vs 113 (island
L2) and 237 vs 181 (ogrid L8) iterations at 1e6 in pilot 3356132.

[V] sCluster pilot 3356132, `node10`, 16 cores (not exclusive), p = 4,
chi = 1e6, p-MG/LDLT solve wall time in seconds (setup excluded):

| threads | 1 | 2 | 4 | 8 | 16 |
|---|---|---|---|---|---|
| ogrid L16, 20609 DOFs, 248 its | 5.16 | 2.78 | 1.67 | 1.11 | 1.09 |
| island L4, 16449 DOFs, 186 its | 3.03 | 1.69 | 1.01 | 0.69 | 0.70 |

Peak RSS at 16 threads is 302,508 KiB (ogrid) and 231,540 KiB (island). These
2D problems are too small to scale beyond eight threads.

[V] Full sweep, Slurm array **3356137** (tasks p = 2/4/6 on `node6`, `node13`,
`node4`, exclusive, 16 threads, revision 3c2fe688, GNU Fortran 12.2.0, `-O3
-Wall -Wextra -fcheck=bounds`): meshes island L1..4 and ogrid L4..16, chi =
1e3, 1e5, 1e7, 1e9, solvers ldlt, amg and direct, 288 rows without failure.
The p = 4 rows reproduce the faepmac1 numbers above digit for digit. The
largest mesh of each family (island L4 against ogrid L16, 0.8 times the DOFs)
gives the direct-solve rel_l2 and the p-MG/LDLT iterations (* = not converged
in 500):

| p | DOFs island/ogrid | interp island/ogrid | chi | rel_l2 island | rel_l2 ogrid | its island/ogrid |
|---|---|---|---|---|---|---|
| 2 | 4129/5185 | 8.65e-3/7.01e-3 | 1e3 | 8.79e-3 | 6.61e-3 | 21/8 |
| 2 | | | 1e5 | 8.68e-3 | 7.87e-3 | 36/53 |
| 2 | | | 1e7 | 1.22e-2 | 3.49e-2 | 167/166 |
| 2 | | | 1e9 | 8.07e-2 | 2.27e-1 | 500*/500* |
| 4 | 16449/20609 | 1.98e-3/1.83e-3 | 1e3 | 1.96e-3 | 1.75e-3 | 34/9 |
| 4 | | | 1e5 | 1.91e-3 | 1.80e-3 | 59/81 |
| 4 | | | 1e7 | 2.17e-3 | 4.14e-3 | 500*/500* |
| 4 | | | 1e9 | 4.82e-3 | 1.74e-2 | 500*/500* |
| 6 | 36961/46273 | 7.82e-4/7.30e-4 | 1e3 | 1.05e-3 | 6.62e-4 | 43/12 |
| 6 | | | 1e5 | 1.01e-3 | 6.69e-4 | 72/99 |
| 6 | | | 1e7 | 1.10e-3 | 1.24e-3 | 500*/500* |
| 6 | | | 1e9 | 1.30e-3 | 3.72e-3 | 500*/500* |

[V] Up to chi = 1e5 the unaligned O-grid is at its interpolation error and
matches or beats the island mesh at matched DOFs (criterion 3: alignment is
not needed there). At 1e9 the O-grid error exceeds its interpolation error
by 32, 9.5 and 5.1 times for p = 2, 4 and 6, and the island mesh by 9.3, 2.4
and 1.7 times. High order therefore reduces the pollution on both meshes, and
the island-adapted mesh is 2.8 to 3.6 times more accurate at 1e9. [V] No
solver here gives a fast stiff-anisotropy solve: p-MG stalls from 1e7 (p = 4,
6), and the direct LDLT reference takes 125 s (ogrid L16) and 95 s (island L4)
at p = 6. Wall times for every row are in the sweep logs.

### Reproduction and limits

```bash
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=1 MPLBACKEND=Agg
cmake -S . -B build/lab -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_Fortran_FLAGS='-O3 -Wall -Wextra -fcheck=bounds' \
  -DFORTFEM_BUILD_TESTING=ON -DFORTFEM_BUILD_BENCHMARKS=ON
cmake --build build/lab -j16 --target test_lab_unstructured_mesh \
  fortfem_lab_uns1_island
ctest --test-dir build/lab -R 'test_lab_unstructured|lab_uns1' --output-on-failure
# KIND LEVEL P RATIO COARSE: one CSV row, header in the sweep wrapper
build/lab/fortfem_lab_uns1_island island 4 4 1e9 direct
build/lab/fortfem_lab_uns1_island ogrid 16 4 1e6 ldlt
```

[V] Wrappers and logs live outside Git at `scluster:~/runs/fortfem-uns1-20261023/`
(`uns1-pilot.sh`, `uns1-sweep.sh`, `uns1_table.py`, `pilot-3356132.log`,
`sweep-3356137_{0,1,2}.log`, `SHA256SUMS`). `python3 -I uns1_table.py
sweep-3356137_*.log` prints every row grouped by p, chi and mesh.
The faepmac1 rows were produced by the same binary arguments with
`OMP_NUM_THREADS=4`.
The island mesh uses a polar chart about the disk centre only away from the
centre (the core is a transfinite square), and it assumes one ridge of the
level function along every ray. The comparison is 2D (single helicity). A 3D
multi-helicity or stochastic case, a chi-robust preconditioner and
time-to-solution against VMEC-class codes are not covered.
<!-- UNS1 END -->

<!-- FB1 BEGIN -->
## FB1: tensor FEM interior with 3D BEM free boundary (criterion 1)

[V] Synthetic data only (point charges, quadratic interior field, prescribed
surface current). Correctness rows from revision 1cb4b23e plus the
allocatable-work-array edit to `parametric_bem_tensor.f90`, workstation,
GNU Fortran 16.2.1, `fpm test` default debug profile (no OpenMP, loaded
machine), so the time columns below are orders of magnitude only. Timing and
thread scaling are cluster rows in the section below.

### Transmission oracle (`test_lab_fem_bem_transmission`)

Exact transmission data on the discrete O-grid torus (5-block O-grid x phi,
`n` radial cells per block, `nphi` toroidal cells), interior and exterior fields
from point charges. Columns: relative volume L2 and H1 error, relative coupled
Neumann error (`relq`) and DtN-only Neumann error (`relq_dtn`), PCG iterations
for relative tolerance 1e-11 with the block preconditioner.

| n | nphi | p | pf | dofs | trace | its | relL2 | relH1 | relq | relq_dtn |
|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| 1 | 6 | 2 | 1 | 300 | 96 | 18 | 3.04e-2 | 1.28e-1 | 3.78e-1 | 3.65e-1 |
| 1 | 6 | 3 | 2 | 936 | 216 | 18 | 9.23e-3 | 5.38e-2 | 1.32e-1 | 1.56e-1 |
| 1 | 6 | 4 | 3 | 2136 | 384 | 22 | 3.38e-3 | 2.71e-2 | 1.16e-1 | 1.22e-1 |
| 2 | 12 | 2 | 1 | 2136 | 384 | 20 | 6.99e-3 | 4.86e-2 | 9.08e-2 | 1.07e-1 |

[V] All errors decrease under p- and h-refinement. The data are too coarse and
the BEM quadrature error too large to confirm an asymptotic rate of p+1 in h
or exponential convergence in p. Neumann errors stagnate near 1e-1 for p=3..4,
so the BEM quadrature limit, not the FEM space, dominates there [H].

[V] Schur scale sigma of the preconditioner (n=1, nphi=6, p=3, pf=2):

| sigma | 1 | 2 | 3 | 5 | 8 |
|--:|--:|--:|--:|--:|--:|
| PCG its | 29 | 20 | 18 | 24 | 30 |

### Vacuum-field current sheet (`test_lab_fem_bem_current_sheet`)

Prescribed divergence-free surface current on the torus boundary, exterior B from
the representation formula (`ext_rel_max`) and interior B on the axis core cells
(`core_rel_l2`) against an independent Biot-Savart surface integral.

| n | nphi | p | dofs | trace | its | ext_rel_max | core_rel_l2 |
|--:|--:|--:|--:|--:|--:|--:|--:|
| 1 | 6 | 2 | 300 | 96 | 17 | 1.57e-2 | 2.46e-2 |
| 1 | 6 | 4 | 2136 | 384 | 20 | 6.69e-4 | 3.84e-3 |
| 2 | 12 | 2 | 2136 | 384 | 18 | 8.51e-4 | 4.79e-3 |

### Thread scaling and size scan (sCluster job 3356163)

[V] Revision 46fefb79, node6 (exclusive), GNU Fortran 12.2.0, `fpm test --target
test_lab_fem_bem_transmission --profile release --flag "-fopenmp
-finline-matmul-limit=0" --link-flag -fopenmp`, run as
`OMP_NUM_THREADS=T OPENBLAS_NUM_THREADS=T <bin> n nphi 4 3 0 0 3` (p=4, pf=3,
default BEM quadrature, sigma=3), so the dense LU/Cholesky of fortnum/OpenBLAS is
inside the scaling. Times are the per-stage wall times printed by the test, in seconds,
`wall` is the whole process (`/usr/bin/time`). Rows: `fb1_sweep_3356163.csv`, log sha256
8e5599000c9d6b4bbb643866bcf08755666afecaa058ea3f38a68e0b7864db2b.

| n | nphi | trace | T | t_bem | t_dtn | t_mg | wall | peak RSS MiB | its |
|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| 2 | 12 | 1536 | 1 | 12.48 | 1.66 | 0.083 | 15.1 | 223 | 28 |
| 2 | 12 | 1536 | 2 | 8.48 | 1.59 | 0.061 | 10.7 | 232 | 28 |
| 2 | 12 | 1536 | 4 | 5.36 | 1.57 | 0.051 | 7.4 | 252 | 28 |
| 2 | 12 | 1536 | 8 | 2.70 | 1.55 | 0.046 | 4.6 | 308 | 28 |
| 2 | 12 | 1536 | 16 | 1.00 | 1.56 | 0.043 | 2.9 | 408 | 28 |
| 3 | 18 | 3456 | 1 | 28.61 | 22.88 | 0.497 | 55.7 | 1009 | 31 |
| 3 | 18 | 3456 | 2 | 20.00 | 22.44 | 0.423 | 45.1 | 1011 | 31 |
| 3 | 18 | 3456 | 4 | 12.25 | 22.06 | 0.383 | 36.2 | 1019 | 31 |
| 3 | 18 | 3456 | 8 | 6.61 | 21.97 | 0.368 | 30.1 | 1066 | 31 |
| 3 | 18 | 3456 | 16 | 2.59 | 21.86 | 0.385 | 25.9 | 1171 | 31 |

[V] `t_bem` scales 12.5x (n=2) and 11.1x (n=3) from 1 to 16 threads. `t_dtn`
does not scale at all (1.66 s to 1.56 s and 22.9 s to 21.9 s), and `t_mg` gains
only 1.9x (n=2) and 1.3x (n=3). The wall time therefore saturates at 5.1x
(n=2) and 2.2x (n=3) because the serial dense DtN build is Amdahl-limited.

| n | nphi | dofs | trace N | its | t_bem | t_dtn | t_mg | wall | peak RSS MiB |
|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| 1 | 6 | 2136 | 384 | 22 | 0.49 | 0.021 | 0.010 | 0.7 | 163 |
| 2 | 12 | 16176 | 1536 | 28 | 1.01 | 1.57 | 0.044 | 3.0 | 380 |
| 3 | 18 | 53640 | 3456 | 31 | 2.61 | 21.78 | 0.386 | 25.8 | 1170 |
| 4 | 24 | 126048 | 6144 | 33 | 5.74 | 127.5 | 2.32 | 138.0 | 3109 |

Size scan at 16 threads (same job, same flags). [V] `t_dtn` grows like N^3
(0.021, 1.57, 21.8, 127.5 s at N=384, 1536, 3456, 6144: ratios 75, 14, 5.9
against N^3 ratios 64, 11.4, 5.6), while `t_bem` grows sub-quadratically at 16
threads (0.49, 1.01, 2.61, 5.74 s, local exponents 0.5 to 1.2) and `t_mg` reaches
2.3 s at N=6144. Dense DtN construction dominates the setup above
N~1200 (interpolated between N=384 and N=1536) and is 92% of the wall time at
N=6144. Below N~1200 and in serial (1 thread: `t_bem` 12.5 s against `t_dtn` 1.7 s
at N=1536) the dense BEM assembly dominates. [H] The DtN build, i.e. the dense solve
and the transposed-jump product in `build_exterior_dtn`, runs on one thread
(traced in the follow-up PR); ACA/FMM would not remove it.

### Reproduction

```bash
fpm test --target test_lab_fem_bem_transmission
fpm test --target test_lab_fem_bem_current_sheet
# scan row: n nphi p [flux_degree target_points fan_points schur_scale]
build/gfortran_*/test/test_lab_fem_bem_transmission 1 6 3 2 0 0 3
```

[V] Not run: `tools/codegen/check_generated.sh` (needs the pinned FortSym
checkout, absent here); the generated file is unchanged by the edit.
[V] Not registered in CMake: the BEM sources are built by fpm only.
<!-- FB1 END -->
