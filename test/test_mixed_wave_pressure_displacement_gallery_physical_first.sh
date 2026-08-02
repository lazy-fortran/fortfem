#!/usr/bin/env bash
set -euo pipefail

repository_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
output_directory="$repository_dir/output/example/"\
"mixed_wave_pressure_displacement_gallery"
log_file="${TMPDIR:-/tmp}/fortfem-mixed-wave-pressure-displacement.log"

rm -rf -- "$output_directory"
(
    cd "$repository_dir"
    timeout --foreground 10s fo exec --no-build \
        mixed_wave_pressure_displacement_gallery >"$log_file"
)

for image in solution_1d.png solution_2d.png solution_3d.png \
    diagnostics_energy_error.png diagnostics_structure.png; do
    test -s "$output_directory/$image"
    test "$(file -b "$output_directory/$image" | cut -d' ' -f1)" = PNG
done
test -s "$output_directory/solution_1d.csv"
test -s "$output_directory/solution_2d.csv"
test -s "$output_directory/diagnostics.csv"
test -s "$output_directory/gallery_sequence.txt"

mapfile -t stages <"$output_directory/gallery_sequence.txt"
test "${#stages[@]}" -eq 2
test "${stages[0]}" = physical_solution
test "${stages[1]}" = diagnostics

python3 - "$output_directory" <<'PY'
import csv
import math
import pathlib
import sys

root = pathlib.Path(sys.argv[1])
values = {}
with (root / "diagnostics.csv").open(encoding="utf-8") as stream:
    next(stream)
    for line in stream:
        name, value = line.rstrip().split(",", 1)
        if name != "primary_plot":
            values[name] = float(value)
assert values["maximum_exact_state_error"] < 2.0e-3
assert values["maximum_relative_energy_drift"] < 3.0e-12
assert values["reversal_error"] < 3.0e-12
assert values["symplectic_map_defect"] < 5.0e-11
assert values["midpoint_wall_time_seconds"] < 10.0

with (root / "solution_1d.csv").open(newline="", encoding="utf-8") as stream:
    rows = list(csv.DictReader(stream))
assert len(rows) == 121
assert max(abs(float(row["pressure"])) for row in rows) > 0.1
assert max(abs(float(row["displacement_x"])) for row in rows) > 0.02
assert max(abs(float(row["momentum"])) for row in rows) > 0.01

with (root / "solution_2d.csv").open(newline="", encoding="utf-8") as stream:
    rows = list(csv.DictReader(stream))
assert len(rows) == 13 * 13
assert max(abs(float(row["displacement_x"])) for row in rows) > 0.02
assert max(abs(float(row["displacement_y"])) for row in rows) > 0.01
assert max(abs(float(row["pressure"])) for row in rows) > 0.05
print("mixed pressure/displacement gallery emits physical 1-D/2-D/3-D fields")
PY
