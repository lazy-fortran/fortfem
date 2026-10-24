# Discretisation laboratory

The discretisation-independent physical definitions and limitations are in
[PROBLEMS.md](PROBLEMS.md). All inputs are synthetic or published. Every result
is associated with a decision criterion: 1 FEM+BEM interior, 2 magnetic-axis
regularity, 3 optional field alignment, or 4 time against VMEC-class solvers.
A criterion tag records relevance, not proof that the criterion is satisfied.

## Build and verify

Run long builds on a remote workstation or in a Slurm allocation. Use the project's
native Fortran runner and explicit resource limits, for example:

```sh
export FO_JOBS=4 OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 MPLBACKEND=Agg
export PKG_CONFIG_PATH=/opt/homebrew/opt/openblas/lib/pkgconfig
unset DISPLAY
fo test --profile release --flag -fopenmp \
  test_lab_problem_definitions test_pcg_scale_invariance test_advanced_solvers
fo build --profile release --flag -fopenmp
```

On GNU Fortran older than 13, append `-finline-matmul-limit=0`, as required by
the root CMake contract. On an offload-enabled GCC (the cluster's 12.2), compile
and link OpenMP code with `-foffload=disable`, or the link fails with an
undefined `GOMP_offload_unregister_ver`. Never run a bare `fpm test` or a full test suite here.
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
`import_n1.py --admit` supplies them for the thirty retained N1 records.
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
1e-3, 1e-5, 1e-7, with no extrapolation. The original axis norm uses `s<0.04`,
a physical-volume-weighted tube. `lab_mf_run` also measures the physical ball
of radius 0.2 about the axis point at phi=0. These are indicator-quadrature
norms, not certified integrals. The polar comparison collapses all r=0 angular
nodes at each phi, which is an explicit axis constraint used only in that
comparison. The primary GLL path uses MF1's regular O-grid, warped by the
analytic L1 map, and needs no axis collapse.

The stacked MF consumer accepts any positive GLL degree and uses p-multigrid:

```sh
python3 benchmark/lab/L1.py --family gll-mf --p 10 --n 1 \
  --preconditioner pmg --threads 4 --executable /absolute/lab_mf_run \
  --flags '-O3 -funroll-loops -fopenmp' --output mf.csv
LAB_MF_EXE=/absolute/lab_mf_run LAB_FLAGS='-O3 -funroll-loops -fopenmp' \
  OMP_NUM_THREADS=4 bash benchmark/lab/mf_baselines.sh mf.csv
```

L5's `gll-polar` comparison uses the same operator's assembled oracle with
`ichol_row`, because generic p-coarsening of collapsed polar cells is outside
MF1's topology contract. A polar run has n radial, 8n angular and 20n toroidal
cells. The regular O-grid has five cross-section patches and 20n toroidal cells.
The MF load uses the solver's quadrature, and error integration uses p+3 points
per direction. The initial MF adapter left `nnz=0` unreported for the polar
runs; their memory is measured by RSS. Current runs report polar fine-matrix
NNZ, while matrix-free rows use zero for an absent fine matrix. Coarse storage
is included in RSS. An aligned periodic sheared mesh still needs an interface
constraint and is explicitly rejected by the driver.

The existing N1 sweep can be imported without rerunning it or copying cases:

```sh
# N1_ROOT is the retained N1 sweep tree on the machine that holds it.
python3 -I benchmark/lab/import_n1.py "$N1_ROOT" --output l6-n1-reference.csv
```

The importer checks the retained adapter SHA256 manifest. Its thirty rows
include failed GVEC integer-family runs, and keep native solve, adapter and
whole-process times separate. They retain N1's original error normalizations.
With `--admit` it also evaluates the four-component L6 admission vector of
PROBLEMS.md against the finest VMEC++ reference of each case, reading only the
retained `samples.txt`, `metrics.json` and `wout_reference.nc` files. This needs
numpy, scipy and SIMSOPT, so it runs as `slurm/l6_admission.sbatch` with the N1
Python; `l6_admission_table.py` renders the tables of RESULTS.md from its CSV.
Prebuilt external-package compiler flags remain explicitly unavailable.

## Cluster recipes

`slurm/pilot.sbatch` builds and runs the ladder on one four-core allocation.
`--mem` is documentation on sCluster, so the mandatory external memory guard
uses the reserved-core share, 4*11.8 GiB. Copy the cluster skill's guard into
the run root, set LAB_REPO, LAB_REVISION and MEM_GUARD, and submit once.
The pilot checks a clean immutable revision before building. Archive the job
ID, logs, CSV and checksums. Do not wait on a queued pilot. Both templates
write their Slurm log to `/home/%u/ben1-runs/%x-%j.log`, outside the worktree,
so the cleanliness check cannot see it. Run `mkdir -p ~/ben1-runs` before the
first submission, because Slurm does not create the log directory.

`slurm/scaling.sbatch` serves both strong and weak scaling for a selected
family. Set LAB_MODE=strong or weak and LAB_THREADS to 1,2,4,8,16,32,64.
Override `--cpus-per-task` to at least the thread count and
`ceil(1.5*measured_peak_GiB/11.8)`. Set LAB_PILOT_PEAK_GIB from a completed
pilot, not an estimate. The template validates that reservation and guards
its core share. Weak scaling holds p fixed and uses `n proportional to
threads^(1/3)`, reporting actual sizes because integer rounding matters.
No scaling jobs are submitted before the pilot qualifies memory.
