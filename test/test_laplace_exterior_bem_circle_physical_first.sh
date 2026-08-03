#!/usr/bin/env bash
set -euo pipefail

repository_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
output_directory="$repository_dir/output/example/laplace_exterior_bem_circle"
log_file="${TMPDIR:-/tmp}/fortfem-laplace-exterior-bem.log"

rm -rf -- "$output_directory"
(
    cd "$repository_dir"
    timeout --foreground "${EXAMPLE_TIMEOUT:-10}s" \
        fpm run --example laplace_exterior_bem_circle >"$log_file"
)

test -s "$output_directory/exterior_laplace_solution_2d.png"
test -s "$output_directory/exterior_laplace_trace_1d.png"
test -s "$output_directory/exterior_laplace_density_2d.png"
test -s "$output_directory/solution.csv"
test -s "$output_directory/provenance.json"
test -s "$output_directory/gallery_sequence.txt"
test "$(file -b "$output_directory/exterior_laplace_solution_2d.png" | cut -d' ' -f1)" = PNG
test "$(head -n 1 "$output_directory/gallery_sequence.txt")" = physical_solution
test "$(sed -n '2p' "$output_directory/gallery_sequence.txt")" = diagnostics

primary_name=$(awk '$1 == "laplace_exterior_bem_circle" {print $2}' \
    "$repository_dir/doc/examples/primary_plots.txt")
test "$primary_name" = exterior_laplace_solution_2d.png

python3 - "$output_directory" <<'PY'
import csv
import json
import math
import pathlib
import sys

root = pathlib.Path(sys.argv[1])
metadata = json.loads((root / "provenance.json").read_text(encoding="utf-8"))
assert metadata["schema"] == "fortfem-gallery-provenance-1"
assert metadata["problem"] == "exterior-laplace-unit-circle"
assert metadata["solution_first"] is True
assert metadata["exact_solution"] == "cos(theta)/r"
assert metadata["external_data"] is False
assert metadata["primary_plot"] == "exterior_laplace_solution_2d.png"

with (root / "solution.csv").open(newline="", encoding="utf-8") as stream:
    rows = list(csv.DictReader(stream))
assert len(rows) > 4000
max_error = 0.0
for row in rows[::113]:
    x = float(row["x"])
    y = float(row["y"])
    radius2 = x * x + y * y
    assert radius2 > 1.02 * 1.02
    computed = float(row["computed"])
    exact = x / radius2
    assert abs(float(row["exact"]) - exact) < 1.0e-13
    max_error = max(max_error, abs(computed - exact))
assert max_error < 1.0e-3
PY

echo "exterior Laplace BEM renders the physical field first with an independent solution oracle"
