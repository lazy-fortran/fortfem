#!/usr/bin/env bash
set -euo pipefail

repo_root=$(cd "$(dirname "$0")/.." && pwd)
cd "$repo_root"
fo exec solver_benchmark >/tmp/fortfem_solver_benchmark_gallery.log

output_dir=output/example/solver_benchmark
test "$(head -n 1 "$output_dir/gallery_sequence.txt")" = physical_solution
test -s "$output_dir/poisson_solution_2d.png"
test -s "$output_dir/poisson_solution.csv"

python3 - "$output_dir/poisson_solution.csv" <<'PY'
import csv
import math
import pathlib
import sys

rows = list(csv.DictReader(pathlib.Path(sys.argv[1]).open()))
assert len(rows) == 32 * 32
errors = []
for row in rows:
    x, y = float(row["x"]), float(row["y"])
    numerical = float(row["numerical"])
    exact = float(row["exact"])
    assert all(math.isfinite(float(row[k])) for k in ("x", "y", "numerical", "exact"))
    expected = 0.0
    for m in range(1, 80, 2):
        for n in range(1, 80, 2):
            expected += (16.0 / math.pi**4 * math.sin(m * math.pi * x) *
                         math.sin(n * math.pi * y) /
                         (m * n * (m*m + n*n)))
    assert abs(exact - expected) < 2.0e-10
    errors.append(abs(numerical - exact))
assert max(errors) < 2.0e-3
assert max(errors) > 1.0e-8, "plot must contain solved FEM values, not copied exact data"
print("solver benchmark gallery renders the solved Poisson field before timing diagnostics")
PY
