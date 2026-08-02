#!/usr/bin/env bash
set -euo pipefail

repo_root=$(cd "$(dirname "$0")/.." && pwd)
cd "$repo_root"
fo exec mixed_poisson >/tmp/fortfem_mixed_poisson_gallery.log

output_dir=output/example/mixed_poisson
test "$(head -n 1 "$output_dir/gallery_sequence.txt")" = physical_solution
test -s "$output_dir/mixed_poisson_solution_2d.png"
test -s "$output_dir/mixed_poisson_solution.csv"

python3 - "$output_dir/mixed_poisson_solution.csv" <<'PY'
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
    errors.append(abs(numerical - exact))
assert max(errors) < 2.0e-2
assert max(errors) > 1.0e-8, "plot must contain solved DG values, not copied exact data"
print("mixed Poisson gallery renders the solved pressure before diagnostics")
PY
