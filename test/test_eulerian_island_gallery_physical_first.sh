#!/usr/bin/env bash
set -euo pipefail

repository_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
output_directory="$repository_dir/output/example/eulerian_island_gallery"
log_file="${TMPDIR:-/tmp}/fortfem-eulerian-island-gallery.log"

rm -rf -- "$output_directory"
(
    cd "$repository_dir"
    timeout --foreground "${EXAMPLE_TIMEOUT:-30}s" \
        fpm run --example eulerian_island_gallery >"$log_file"
)

test -s "$output_directory/island_flux_solution_2d.png"
test -s "$output_directory/island_flux_diagnostics_2d.png"
test -s "$output_directory/island_flux.csv"
test -s "$output_directory/provenance.json"
test -s "$output_directory/gallery_sequence.txt"
test "$(file -b "$output_directory/island_flux_solution_2d.png" | cut -d' ' -f1)" = PNG

mapfile -t stages <"$output_directory/gallery_sequence.txt"
test "${#stages[@]}" -eq 2
test "${stages[0]}" = physical_solution
test "${stages[1]}" = diagnostics

python3 - "$output_directory" <<'PY'
import csv
import json
import math
import pathlib
import sys

root = pathlib.Path(sys.argv[1])
with (root / "provenance.json").open(encoding="utf-8") as stream:
    record = json.load(stream)
assert record["schema"] == "fortfem-eulerian-island-gallery-v1"
assert record["closure"] == "neutral-caller-owned-force-divergence"
assert record["primary_plot"] == "island_flux_solution_2d.png"
assert record["event_code"] == 1
assert record["event_index"] == 1
assert record["maximum_divergence"] == 0.0
assert record["residual_l2"] > 0.0
assert math.isfinite(record["residual_jvp_l2"])

parameter = record["island_amplitude"]
with (root / "island_flux.csv").open(newline="", encoding="utf-8") as stream:
    rows = list(csv.DictReader(stream))
assert len(rows) == record["grid_x"] * record["grid_y"]
for row in rows[::113]:
    x = float(row["x"])
    y = float(row["y"])
    psi = float(row["psi"])
    bx = float(row["bx"])
    by = float(row["by"])
    divergence = float(row["divergence"])
    assert abs(psi - ((x * x - parameter) ** 2 + y * y)) < 2.0e-12
    assert abs(bx - 2.0 * y) < 2.0e-12
    assert abs(by + 4.0 * x * (x * x - parameter)) < 2.0e-12
    assert abs(divergence) < 1.0e-14

print("Eulerian island gallery renders physical flux before diagnostics")
PY
