# Native sparse migration: unchanged FortFEM tests

Measured 2026-10-16 on an AMD Ryzen 9 5950X, GNU Fortran 16.2.1
(20260810), one thread, CPU affinity 0. Source checkpoints:

- Baseline FortFEM `72a2c44bac62d61f08541359c79590b1dfcf99ab`.
- Native FortFEM `4ed292ffac67d6d8dc7357e93942493b9a01d6c3`.
- Both use FortSparse `8dfac488ba5384f7a577beeac2d6f7c23ba6ac16` to
  isolate the migration. The older FortSparse revisions originally declared
  by baseline CMake/fpm were not measured.
- Both compile the exact original baseline test sources. The original adapter
  selects AUTO, which resolves to system SuperLU 7; the new adapter selects
  NATIVE. All other optional sparse backends and AMD/METIS are disabled.
- Both select system OpenBLAS. Prebuilt SuperLU additionally depends on system
  `libblas.so.3`; its loaded dependency cost is included in the baseline.

Release flags are `-O2 -g -fcheck=all -fbacktrace` for Fortran and
`-O3 -DNDEBUG` for C. `OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=MKL_NUM_THREADS=1`.
The wrapper builds the ordinary FortFEM shared C API and five unchanged test
programs. Both cohorts first pass their five registered CTests.

Nine measured processes per variant and case, plus one excluded warmup.
Variant order alternates within each pair. All 90 measured invocations and
10 warmups pass the unchanged independent test assertions and take <1 s each.
The additional native-only pivot/lifetime/failure assertions added by the PR
have no original equivalent and are excluded from this matched comparison.

The table gives medians. CPU is user+system from `wait4`; peak RSS is each
individual process's Linux `ru_maxrss` in KiB. A small C observer uses
`posix_spawn` followed by `wait4`, avoiding the Python fork/high-water floor.
Wall covers spawn through completion. These are complete tiny test-program
measurements, including dynamic loading, test oracles and output. They are
not factor-only timings, allocation-payload estimates, concurrent aggregate
RSS, compilation time or a Pareto-front admission. Background unrelated
services remained present; other owned CPU builds/tests were quiescent during
the short measurement lease. Min/max, separate user/system CPU and every trial
are retained in the numeric receipts.

| Original test driver (matrix n) | CPU baseline→native ms | Δ CPU ms (%) | Peak RSS baseline→native KiB | Δ RSS KiB (%) | Wall baseline→native ms |
|---|---:|---:|---:|---:|---:|
| `sparse_direct_backend` (2) | 2.242→1.830 | -0.412 (-18.4%) | 9612→6476 | -3136 (-32.6%) | 2.345→1.897 |
| `sparse_direct_derivatives` (3) | 2.159→1.804 | -0.355 (-16.4%) | 9124→6312 | -2812 (-30.8%) | 2.275→1.879 |
| `sparse_direct_complex_derivatives` (3) | 2.201→1.850 | -0.351 (-15.9%) | 9228→6428 | -2800 (-30.3%) | 2.262→1.956 |
| `sparse_constrained_derivatives` (4) | 2.156→1.800 | -0.356 (-16.5%) | 9024→6324 | -2700 (-29.9%) | 2.245→1.898 |
| `capi_solver` (2) | 2.143→1.773 | -0.370 (-17.3%) | 8992→6344 | -2648 (-29.4%) | 2.234→1.858 |

## Reproduction

Create detached checkouts at the three revisions above. Configure this directory
twice with `FORTFEM_SOURCE` pointing to the baseline/native checkout,
`DRIVER_SOURCE` always pointing to the original baseline, and
`FETCHCONTENT_SOURCE_DIR_FORTSPARSE` pointing to the common sparse checkout:

```sh
cmake -S <this-directory> -B <build> \
  -DFORTFEM_SOURCE=<baseline-or-native> -DDRIVER_SOURCE=<baseline> \
  -DFETCHCONTENT_SOURCE_DIR_FORTSPARSE=<sparse-8dfac48> \
  -DCMAKE_BUILD_TYPE=Release \
  '-DCMAKE_Fortran_FLAGS_RELEASE=-O2 -g -fcheck=all -fbacktrace' \
  -DFORTSPARSE_USE_SYSTEM_DEPS=ON -DFORTSPARSE_ENABLE_SUPERLU=<ON-or-OFF> \
  -DFORTSPARSE_ENABLE_MUMPS=OFF -DFORTSPARSE_ENABLE_UMFPACK_IPC=OFF \
  -DFORTSPARSE_ENABLE_AMD=OFF -DFORTSPARSE_ENABLE_METIS=OFF
cmake --build <build> -j2
ctest --test-dir <build> --output-on-failure -j2
cc -O2 -Wall -Wextra measure.c -o <measure>
python3 run.py --baseline-build <baseline-build> --native-build <native-build> \
  --measure <measure> --output <results> --repeats 9
```

The native checkout enables its native-only default; baseline SuperLU must be
ON, native OFF. For this receipt, binaries and shared libraries were copied to
immutable cohort directories and executable RPATHs set to `$ORIGIN` before
measurement. Binary hashes are in `manifest.json`. Frozen ELFs, logs and loaded
library listings remain outside Git under
`/home/ert/code/sparse-dev/benchmarks/2026-10-16-fortfem-native-migration/`.
