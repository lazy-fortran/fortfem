# Discretisation laboratory

The discretisation-independent physical definitions and limitations are in
[PROBLEMS.md](PROBLEMS.md). All inputs are synthetic or published. Every result
is associated with a decision criterion: 1 FEM+BEM interior, 2 magnetic-axis
regularity, 3 optional field alignment, or 4 time against VMEC-class solvers.
A criterion tag records relevance, not proof that the criterion is satisfied.

## Build and verify

Run long builds on faepmac1/2 or in a Slurm allocation. Use the project's native
Fortran runner and explicit resource limits, for example on faepmac2:

```sh
export FO_JOBS=4 OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 MPLBACKEND=Agg
export PKG_CONFIG_PATH=/opt/homebrew/opt/openblas/lib/pkgconfig
unset DISPLAY
fo test --profile release --flag -fopenmp \
  test_lab_problem_definitions test_pcg_scale_invariance test_advanced_solvers
fo build --profile release --flag -fopenmp
```

On GNU Fortran older than 13, append `-finline-matmul-limit=0`, as required by
the root CMake contract. Never run a bare `fpm test` or a full test suite here.
`check_definitions.py` is an additional exact-algebra independent oracle
(requires SymPy). The native test checks sampled Cartesian differential
identities and through-axis inverse geometry using finite differences.
Finite differences have stated truncation tolerances; the exact symbolic
identities vanish algebraically, not merely at a loose FD tolerance.

The pointwise exact fields, forcing and geometry are emitted by
`tools/codegen/app/gen_lab_problems.f90` against `tools/codegen/fortsym.lock`.
The root generation and comparison scripts include these six products.
Do not edit generated sources. Mac code generation additionally needs the
SymEngine/FLINT/MPFR headers and libraries; the run used Homebrew and an
explicit LLVM libc++ search path, recorded in RESULTS.md.

## Measured rows

Use one Python entry point per problem, after building once:

```sh
python3 benchmark/lab/L1.py --family tet --p 2 --h 0.5 \
  --preconditioner ichol_row --threads 4 --executable /absolute/lab_scalar_run \
  --flags '-O3 -funroll-loops -fopenmp' --output lab.csv
python3 benchmark/lab/L3a.py --p 1 --n 32 --ratio 1000000 \
  --preconditioner ldlt --threads 4 --executable /absolute/lab_scalar_run \
  --flags '-O3 -funroll-loops -fopenmp' --output lab.csv
python3 benchmark/lab/L1.py --family bspline --p 4 --n 2 \
  --preconditioner ichol_row --threads 4 --executable /absolute/lab_bspline_run \
  --flags '-O3 -funroll-loops -fopenmp' --output lab.csv
```

L1, L3, L3a, L4 and L5 use `lab_scalar_run`; L2 uses `lab_hcurl_run` and accepts
Nedelec order one. L1/L5 spline runs use `lab_bspline_run`, with degrees 2–4.
L4 diagonalises its constant component block and combines the modal norms.
L6.py calls the authoritative KIN6D N1 adapter and records an **unqualified**
external timing until all four matched-accuracy metrics have been supplied.
It neither copies the cases nor invents equilibrium accuracy.

`ichol` retains the old dense reference setup. `ichol_row` explicitly selects
the existing row ICHOL primitive, drop tolerance 1e-4 and maximum 64 retained
strict-lower entries per row. The old aliases are unchanged. A new PCG
regression prevents its old absolute curvature cutoff from stopping a small
scaled RHS before convergence. ILU-PCG is available for comparison, but a
nonsymmetric incomplete factor does not provide the SPD theory of ICC-PCG.

CSV fields are defined in `run.py:FIELDS`. `h=1/n` is the cross-section mesh
parameter; toroidal cell counts scale with n as well. `size` is the element
count. DOFs include prescribed boundary DOFs. The geometry column distinguishes
planar-faced tetrahedral and spline approximations from the exact domain.
The spline patch uses a regular Cartesian square-to-disk control map, including
the axis, and identifies the toroidal endpoint coefficient planes for C0
periodicity. The exact trace is interpolated in the spline basis; Greville
values are not used as coefficients. Tetrahedral mesh files have a SHA256;
spline runs leave that field empty and are reproduced from n, p and the source.

Each process supplies true wall time, per-child peak RSS (wait4/getrusage),
DOFs, L2 and H1 seminorm or H(curl) error, iterations and an independently
recomputed linear residual. Timing includes assembly, setup, solve and error
integration, and excludes build and Python mesh preparation. Threads are the
requested OpenMP team size. Existing serial setup paths remain serial.
The command and binary SHA256 are retained beside revision, compiler, flags,
host and Slurm job ID. Do not combine different timing scopes without relabeling.

```sh
LAB_EXE=/absolute/lab_scalar_run LAB_SPLINE_EXE=/absolute/lab_bspline_run \
  LAB_FLAGS='-O3 -funroll-loops -fopenmp' OMP_NUM_THREADS=4 \
  bash benchmark/lab/baselines.sh lab.csv
python3 benchmark/lab/summarize.py lab.csv --output assembled-tables.md
```

The summary reports measured DOFs and process time at relative L2 targets
1e-3, 1e-5, 1e-7, with no extrapolation. Axis norms currently use `s<0.04`, a
physical-volume-weighted tube. A polar-chart comparison and an aligned periodic
mesh adapter remain required before drawing their respective conclusions.

## Cluster recipes

`slurm/pilot.sbatch` builds and runs the ladder on one four-core allocation.
`--mem` is documentation on sCluster, so the mandatory external memory guard
uses the reserved-core share, 4*11.8 GiB. Copy the cluster skill's guard into
the run root, set LAB_REPO, LAB_REVISION and MEM_GUARD, and submit once.
The pilot checks a clean immutable revision before building. Archive the job
ID, logs, CSV and checksums. Do not wait on a queued pilot.

`slurm/scaling.sbatch` serves both strong and weak scaling for a selected
family. Set LAB_MODE=strong or weak and LAB_THREADS to 1,2,4,8,16,32,64.
Override `--cpus-per-task` to at least the thread count and
`ceil(1.5*measured_peak_GiB/11.8)`. Set LAB_PILOT_PEAK_GIB from a completed
pilot, not an estimate. The template validates that reservation and guards
its core share. Weak scaling holds p fixed and uses `n proportional to
threads^(1/3)`, reporting actual sizes because integer rounding matters.
No scaling jobs are submitted before the pilot qualifies memory.
