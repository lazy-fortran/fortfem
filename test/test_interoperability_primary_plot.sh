#!/usr/bin/env bash
set -euo pipefail

repository_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
fixture_dir=$(mktemp -d "${TMPDIR:-/tmp}/fortfem-interop-primary.XXXXXX")
trap 'rm -rf "$fixture_dir"' EXIT
cd "$repository_dir"
unset DISPLAY
export MPLBACKEND=Agg OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1

# No external solver records exist: the physical reference must stand alone.
fpm run --example interoperability_benchmarks --profile release \
    --runner 'timeout --kill-after=10s 180s' -- --primary-only "$fixture_dir/plots"
test -s "$fixture_dir/plots/manufactured_poisson_solution_2d.png"
test "$(find "$fixture_dir/plots" -name '*.png' | wc -l)" -eq 1
test ! -e "$fixture_dir/plots/poisson_accuracy.png"
echo "interoperability physical primary plot needs no external solver records"
