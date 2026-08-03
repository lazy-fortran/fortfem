#!/usr/bin/env bash
set -euo pipefail

# Independent behavioral oracle for the physical surface-current gallery.
# The gallery itself calls the canonical parity facade; this gate only checks
# emitted geometry/current data, positivity/orientation, and stage ordering.
repository_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
output_directory="$repository_dir/output/example/sheet_current_surface_gallery"
log_file="${TMPDIR:-/tmp}/fortfem-sheet-current-surface-gallery.log"

(
    cd "$repository_dir"
    timeout --foreground 10s fo exec --no-build --example sheet_current_surface_gallery \
        >"$log_file"
)

test -s "$output_directory/surface_current_solutions_slab_3d.png"
test -s "$output_directory/surface_current_solutions_cylinder_3d.png"
test -s "$output_directory/surface_current_solutions_sphere_3d.png"
test -s "$output_directory/surface_current_solutions_torus_3d.png"
test -s "$output_directory/surface_current.csv"
test -s "$output_directory/geometry_ledger.csv"
test -s "$output_directory/gallery_sequence.txt"
test -s "$output_directory/benchmark.json"

python3 - "$output_directory/surface_current.csv" \
    "$output_directory/geometry_ledger.csv" \
    "$output_directory/gallery_sequence.txt" \
    "$output_directory/benchmark.json" <<'PY'
import csv
import json
import math
import pathlib
import sys

surface_path, ledger_path, sequence_path, metadata_path = map(pathlib.Path, sys.argv[1:])
rows = list(csv.DictReader(surface_path.open(newline="")))
assert len(rows) > 300
names = {row["geometry"] for row in rows}
assert names == {"slab", "cylinder", "sphere", "torus"}
for row in rows:
    values = [float(row[key]) for key in (
        "x", "y", "z", "nx", "ny", "nz", "weight", "kx", "ky", "kz")]
    assert all(math.isfinite(value) for value in values)
    normal = values[3:6]
    current = values[7:10]
    assert values[6] > 0.0
    assert abs(sum(value * value for value in normal) - 1.0) < 2.0e-11
    assert abs(sum(normal[index] * current[index] for index in range(3))) < 2.0e-11

ledger = list(csv.DictReader(ledger_path.open(newline="")))
assert [row["geometry"] for row in ledger] == ["slab", "cylinder", "sphere", "torus"]
for row in ledger:
    assert int(row["points"]) > 0
    assert float(row["area"]) > 0.0
    assert float(row["relative_error"]) < 4.0e-12
    for key in ("fitted_x", "fitted_y", "fitted_z", "resolved_x", "resolved_y", "resolved_z"):
        assert math.isfinite(float(row[key]))

stages = sequence_path.read_text().splitlines()
assert stages == [
    "physical_solution_slab",
    "physical_solution_cylinder",
    "physical_solution_sphere",
    "physical_solution_torus",
    "diagnostics",
]
metadata = json.loads(metadata_path.read_text())
assert metadata["schema"] == "fortfem-sheet-current-surface-gallery-v1"
assert metadata["physical_solution_first"] is True
assert metadata["geometries"] == ["slab", "cylinder", "sphere", "torus"]
PY

echo "sheet-current gallery emits slab/cylinder/sphere/torus physical surfaces before diagnostics"
