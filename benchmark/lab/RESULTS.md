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

[V] Actual batch wrappers and full logs are outside Git under
`<run root>/fortfem-mf1-20261022/`: `mf1-pilot.sh`, `mf1-sweep.sh`,
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
logs are retained outside Git next to the build worktree.

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

<!-- BEN1 BEGIN -->
## BEN1: problems, harness and curved-torus consumers

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

### Reproduction and pilot

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
a dedicated clean worktree, log `ben1-pilot-3355332.log`, and outputs under
`<run root>/3355332`. It predates the lab's canonical-symmetry fix,
so it is not verification of the final head. The template guards 47.2 GiB,
which is the enforced four-core memory share. Strong/weak scaling templates
are prepared but no scaling campaign was launched before pilot qualification.

[V] KIN6D's N1 cases and adapter are reused by L6.py. Its oracle work is in
[KIN6D PR 16](https://github.com/itpplasma/kin6d/pull/16), already merged.
No matched-accuracy equilibrium speed comparison is claimed here. The matched
B, p, iota and flux contract, including zero-pressure normalization and input
conventions, is in PROBLEMS.md.

### Assembled baselines

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

#### Time to measured accuracy

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

#### Axis refinement diagnostics

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

### BEN1 curved torus MF consumer

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

[V] `test_lab_tensor_mms` exceeded the 10-second Gremlin CPU budget because
its 3D torus ladder (levels 0..2, p = 2..6) cost about 10 of its 11 CPU seconds
on faepmac1 (h = 4 solves 2.4 s at p = 2 and 3.5 s at p = 3; the p = 5 and 6 solves
take 1.0 and 1.4 s). The standard case now runs every 2D ladder and the 3D
levels 0..1 with p = 2..4 in 2.5 CPU seconds, and the full ladder is the
registered `test_lab_tensor_mms_slow`. No asserted rate or factor changed.
The Gremlin gate of the eight selected targets is green.

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

#### Time to measured accuracy

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

#### Axis refinement diagnostics

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

#### Physical axis ball

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

### L6: reused N1 external reference measurements

[V] `l6-n1-reference-20261022.csv` imports thirty existing VMEC++/GVEC records
from N1's completed sCluster job 3355100 on node7. The adapter snapshot is
checked against `sweep-binaries.sha256`; no equilibrium cases or solver
harness are copied. Eighteen VMEC++ processes and four GVEC QA processes
completed; eight GVEC integer-family processes failed and retain blank
accuracy fields. The table below selects one-core points, not a best-code
ranking. Full requested CPU-affinity counts, input hashes, package versions,
DOF conventions and exact commands are retained in the CSV. Failed GVEC DOFs
refer to the last saved native stage, not to a converged state.

[V] Native solve time, adapter solve time and whole-process time have different
scopes. External wheel build compilers/flags were not exposed by N1 and remain
explicitly unavailable. N1 pressure errors use its own relative-pressure norm,
flux errors are maximum label defects, and sample positions belong to each
candidate. The common-physical-sample vector of PROBLEMS.md is evaluated in
the next subsection. The QA case has no exact oracle, so its accuracy is not
filled with a force stopping tolerance.

| Code | Case | m | status / exit | DOFs | native solve (s) | adapter solve (s) | process (s) | peak RSS (MiB) | B sample error | p sample error (N1) | max flux-label defect |
|---|---|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|
| vmecpp | landreman-integer | 8 | sampled / 0 | 14355 | 0.2519 | 1.3704 | 7.29 | 469.62 | 0.0003222 | 0.0053373 | 0.0079158 |
| vmecpp | issan-integer | 8 | sampled / 0 | 28710 | 0.50646 | 1.6327 | 9.22 | 478.97 | 0.000687 | 0.0053031 | 0.010207 |
| gvec | landreman-integer | 4 | incomplete / 1 | 2178 | — | — | 15.79 | 240.93 | — | — | — |
| gvec | issan-integer | 4 | incomplete / 1 | 2178 | — | — | 10.83 | 243.54 | — | — | — |
| gvec | qa | 4 | sampled / 0 | 1089 | 17.087 | 18.882 | 19.44 | 243.41 | — | — | — |

### L6: matched-accuracy admission of the N1 records (criterion 4)

[V] `l6-n1-reference-20261022.csv` now carries the four-component vector of
PROBLEMS.md for every N1 record, from `import_n1.py --admit` run as sCluster job
3355869 (node8, 2:40, peak 606 MiB, revision b892a4a2, `slurm/l6_admission.sbatch`;
the first run 3355859 at 99108211 differed only in the replicate labels). It
reads each record's retained `samples.txt` (1728 equal-weight interior samples,
three surfaces) and inverts the reference coordinates at the sample positions
as N1's `compare_qa_final.py` does, with forward position checks below 3e-12 m.
The reference is the finest VMEC++ run of the case: `m12-t1` for both integer
families and N1's `m16/n12/ns201` run (8 cores, 42 s) for QA. `e_B` is the
Cartesian B relative RMS, `e_p` uses the denominator max(RMS p_ref, B_ref,rms^2/2)
in mu0-scaled pressure, `e_iota` uses max(RMS iota_ref, 1) on the three sample
surfaces in input-sign convention, and `e_flux` is the RMS label defect divided
by the unit flux range. A record is admitted at epsilon only if it converged and
all four components are at most epsilon. Records that are the reference or its
thread replicate are labelled and not admitted.

[V] The integer-family `e_B` against the m12 reference (landreman m4/m8
1.10e-3/3.80e-4, issan m4/m8 6.07e-4/2.53e-4) has the magnitude of N1's
exact-oracle B error (1.12e-3/3.22e-4, 6.93e-4/6.87e-4), except for issan m8,
where the exact-oracle error is 2.7 times larger. This confirms that the units
and signs are aligned. The
m12 records themselves differ from the exact oracle by 4.2e-4 (landreman) and
5.2e-4 (issan), so the m12 reference is not exact either, and the vector is a
difference from a finite-resolution reference, not an error bound.
[V] e_iota is zero to roundoff for the integer families, whose iota is a
constraint. QA has p = 0, so e_p is zero by construction.

[V] Outcome at epsilon = 1e-2 and 1e-3: all eight non-reference VMEC++
integer-family records are admitted at 1e-2 (landreman m4 and m8, issan m4);
issan m8 is also admitted at 1e-3, and landreman m4/m8 are not, because e_flux
(4.3e-3, 1.4e-3) exceeds it. All six VMEC++ QA records are not admitted at 1e-2:
e_flux saturates at 1.2e-2 for m8 and m12 (1.7e-2 at m4) while e_B and e_iota
are at most 1.5e-3 and 1.8e-3 there. The four GVEC QA records are admitted at 1e-3 (maximum
5.3e-4 at m4 and 3.4e-4 at m8). The eight GVEC integer-family records failed and
are not admitted. [H] The VMEC++ QA flux plateau suggests a radial-resolution
floor at the default radial grid relative to the ns201 reference, and the GVEC
agreement with this VMEC++ reference is a cross-code comparison rather than a
refinement study. There is no second finer reference, so no independent
refinement difference exists for QA.

| Code | Case | m | threads | e_B | e_p | e_iota | e_flux | max | admitted |
|---|---|---:|---:|---:|---:|---:|---:|---:|---|
| vmecpp | landreman-integer | 4 | 1 | 1.10e-03 | 2.95e-04 | 2.64e-16 | 4.29e-03 | 4.29e-03 | yes, 1e-2 |
| vmecpp | landreman-integer | 4 | 8 | 1.10e-03 | 2.95e-04 | 2.64e-16 | 4.29e-03 | 4.29e-03 | yes, 1e-2 |
| vmecpp | landreman-integer | 8 | 1 | 3.80e-04 | 9.82e-05 | 1.43e-16 | 1.43e-03 | 1.43e-03 | yes, 1e-2 |
| vmecpp | landreman-integer | 8 | 8 | 3.80e-04 | 9.82e-05 | 1.43e-16 | 1.43e-03 | 1.43e-03 | yes, 1e-2 |
| vmecpp | landreman-integer | 12 | 1 | 1.12e-15 | 2.08e-16 | 0.00e+00 | 3.02e-15 | 3.02e-15 | reference |
| vmecpp | landreman-integer | 12 | 8 | 2.13e-13 | 6.06e-15 | 0.00e+00 | 8.82e-14 | 2.13e-13 | reference |
| vmecpp | issan-integer | 4 | 1 | 6.07e-04 | 1.36e-04 | 2.64e-16 | 1.70e-03 | 1.70e-03 | yes, 1e-2 |
| vmecpp | issan-integer | 4 | 8 | 6.07e-04 | 1.36e-04 | 2.64e-16 | 1.70e-03 | 1.70e-03 | yes, 1e-2 |
| vmecpp | issan-integer | 8 | 1 | 2.53e-04 | 6.78e-05 | 1.43e-16 | 8.51e-04 | 8.51e-04 | yes, 1e-3 |
| vmecpp | issan-integer | 8 | 8 | 2.53e-04 | 6.78e-05 | 1.43e-16 | 8.51e-04 | 8.51e-04 | yes, 1e-3 |
| vmecpp | issan-integer | 12 | 1 | 1.09e-15 | 2.24e-16 | 0.00e+00 | 2.81e-15 | 2.81e-15 | reference |
| vmecpp | issan-integer | 12 | 8 | 4.12e-13 | 5.92e-14 | 0.00e+00 | 7.43e-13 | 7.43e-13 | reference |
| vmecpp | qa | 4 | 1 | 2.20e-03 | 0.00e+00 | 3.39e-03 | 1.71e-02 | 1.71e-02 | no |
| vmecpp | qa | 4 | 8 | 2.20e-03 | 0.00e+00 | 3.39e-03 | 1.71e-02 | 1.71e-02 | no |
| vmecpp | qa | 8 | 1 | 1.46e-03 | 0.00e+00 | 1.67e-03 | 1.18e-02 | 1.18e-02 | no |
| vmecpp | qa | 8 | 8 | 1.46e-03 | 0.00e+00 | 1.67e-03 | 1.18e-02 | 1.18e-02 | no |
| vmecpp | qa | 12 | 1 | 1.49e-03 | 0.00e+00 | 1.82e-03 | 1.20e-02 | 1.20e-02 | no |
| vmecpp | qa | 12 | 8 | 1.49e-03 | 0.00e+00 | 1.82e-03 | 1.20e-02 | 1.20e-02 | no |
| gvec | landreman-integer | 4 | 1 | — | — | — | — | — | no (failed) |
| gvec | landreman-integer | 4 | 8 | — | — | — | — | — | no (failed) |
| gvec | landreman-integer | 8 | 1 | — | — | — | — | — | no (failed) |
| gvec | landreman-integer | 8 | 8 | — | — | — | — | — | no (failed) |
| gvec | issan-integer | 4 | 1 | — | — | — | — | — | no (failed) |
| gvec | issan-integer | 4 | 8 | — | — | — | — | — | no (failed) |
| gvec | issan-integer | 8 | 1 | — | — | — | — | — | no (failed) |
| gvec | issan-integer | 8 | 8 | — | — | — | — | — | no (failed) |
| gvec | qa | 4 | 1 | 1.72e-04 | 0.00e+00 | 2.23e-04 | 5.25e-04 | 5.25e-04 | yes, 1e-3 |
| gvec | qa | 4 | 8 | 1.72e-04 | 0.00e+00 | 2.23e-04 | 5.25e-04 | 5.25e-04 | yes, 1e-3 |
| gvec | qa | 8 | 1 | 2.71e-05 | 0.00e+00 | 1.86e-05 | 3.37e-04 | 3.37e-04 | yes, 1e-3 |
| gvec | qa | 8 | 8 | 2.71e-05 | 0.00e+00 | 1.86e-05 | 3.37e-04 | 3.37e-04 | yes, 1e-3 |

Admitted records only. Cores are the requested CPU-affinity count of the N1
sweep; wall times and peak RSS (process maximum) are N1's, and the three times
have different scopes as stated above. External build compilers and flags
remain unavailable.

| Code | Case | m | cores | native solve (s) | adapter solve (s) | process (s) | peak RSS (MiB) | admitted |
|---|---|---:|---:|---:|---:|---:|---:|---|
| vmecpp | landreman-integer | 4 | 1 | 0.03069 | 1.171 | 5.97 | 434.1 | yes, 1e-2 |
| vmecpp | landreman-integer | 4 | 8 | 0.05762 | 1.22 | 6.02 | 443.0 | yes, 1e-2 |
| vmecpp | landreman-integer | 8 | 1 | 0.2519 | 1.37 | 7.29 | 469.6 | yes, 1e-2 |
| vmecpp | landreman-integer | 8 | 8 | 0.1537 | 1.263 | 7.22 | 482.1 | yes, 1e-2 |
| vmecpp | issan-integer | 4 | 1 | 0.06571 | 1.196 | 6.37 | 436.5 | yes, 1e-2 |
| vmecpp | issan-integer | 4 | 8 | 0.06781 | 1.193 | 6.18 | 445.8 | yes, 1e-2 |
| vmecpp | issan-integer | 8 | 1 | 0.5065 | 1.633 | 9.22 | 479.0 | yes, 1e-3 |
| vmecpp | issan-integer | 8 | 8 | 0.2135 | 1.338 | 9.14 | 489.3 | yes, 1e-3 |
| gvec | qa | 4 | 1 | 17.09 | 18.88 | 19.44 | 243.4 | yes, 1e-3 |
| gvec | qa | 4 | 8 | 14.91 | 16.68 | 17.26 | 251.1 | yes, 1e-3 |
| gvec | qa | 8 | 1 | 58.42 | 60.26 | 60.83 | 471.7 | yes, 1e-3 |
| gvec | qa | 8 | 8 | 24.12 | 25.85 | 26.47 | 482.6 | yes, 1e-3 |

[V] Admitted VMEC++ integer-family records take 0.03-0.51 s native, 5.97-9.22 s
whole-process and 434-489 MiB, and the admitted GVEC QA records take 14.9-58.4 s
native and 243-483 MiB, with the whole-process times of the table. They belong
to different cases, so the table is no code ranking. No lab-family solve is
timed against these records and no VMEC-class speed claim is made.
<!-- BEN1 END -->
