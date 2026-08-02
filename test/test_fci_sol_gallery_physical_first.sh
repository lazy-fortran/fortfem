#!/usr/bin/env bash
set -euo pipefail

repository_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
output_directory="$repository_dir/output/example/fci_sol_gallery"
log_file="${TMPDIR:-/tmp}/fortfem-fci-sol-gallery.log"

rm -rf -- "$output_directory"
(
    cd "$repository_dir"
    timeout --foreground "${EXAMPLE_TIMEOUT:-10}s" \
        fpm run --example fci_sol_gallery >"$log_file"
)

test -s "$output_directory/fci_sol_field_lines_3d.png"
test -s "$output_directory/fci_sol_solution_2d.png"
test -s "$output_directory/fci_sol_gradient_diagnostic_1d.png"
test -s "$output_directory/fci_sol_solution.csv"
test -s "$output_directory/fci_sol_field_lines.csv"
test -s "$output_directory/benchmark.json"
test "$(file -b "$output_directory/fci_sol_field_lines_3d.png" | cut -d' ' -f1)" = PNG

mapfile -t stages <"$output_directory/gallery_sequence.txt"
test "${#stages[@]}" -eq 2
test "${stages[0]}" = physical_solution
test "${stages[1]}" = diagnostics

primary_name=$(awk '$1 == "fci_sol_gallery" {print $2}' \
    "$repository_dir/doc/examples/primary_plots.txt")
test "$primary_name" = fci_sol_field_lines_3d.png

python3 - "$output_directory" <<'PY'
import csv
import json
import math
import pathlib
import sys

root = pathlib.Path(sys.argv[1])
with (root / "benchmark.json").open(encoding="utf-8") as stream:
    record = json.load(stream)
assert record["provenance"] == "analytic-torus-fci-sol-v1"
assert record["closure"] == "neutral-support-operator"
assert record["primary_plot"] == "fci_sol_field_lines_3d.png"
assert record["field_line_count"] == 3
assert abs(record["mass_rate"]) < 5.0e-10
assert record["dissipation"] < 0.0
assert math.isfinite(record["gradient_action_seconds"])

with (root / "fci_sol_solution.csv").open(newline="", encoding="utf-8") as stream:
    rows = list(csv.DictReader(stream))
assert len(rows) == 24 * 48
for row in rows[::37]:
    theta = float(row["theta"])
    phi = float(row["phi"])
    value = float(row["solution"])
    expected = math.sin(theta - 0.42 * phi) + 0.35 * math.cos(2.0 * phi)
    assert abs(value - expected) < 2.0e-13

with (root / "fci_sol_field_lines.csv").open(newline="", encoding="utf-8") as stream:
    line_rows = list(csv.DictReader(stream))
assert len(line_rows) == 3 * (128 + 1)
for row in line_rows[::29]:
    x, y, z = (float(row[name]) for name in ("x", "y", "z"))
    minor_radius = math.hypot(math.hypot(x, y) - 2.4, z)
    # Each traced line stays on its prescribed toroidal surface; this is an
    # independent geometric oracle, not a copy of the Fortran implementation.
    assert 0.27 < minor_radius < 0.64
PY

echo "FCI/SOL gallery renders physical 3D field lines before diagnostics"
