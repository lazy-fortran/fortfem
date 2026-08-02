#!/usr/bin/env bash
set -euo pipefail

repository_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
output_directory="$repository_dir/output/example/mixed_elasticity_wave"
log_file="${TMPDIR:-/tmp}/fortfem-mixed-elasticity-gallery.log"

rm -rf -- "$output_directory"
(
    cd "$repository_dir"
    timeout --foreground "${EXAMPLE_TIMEOUT:-30}s" \
        fpm run --example mixed_elasticity_wave >"$log_file"
)

test -s "$output_directory/mixed_elasticity_solution_1d.png"
test -s "$output_directory/mixed_elasticity_solution_2d.png"
test -s "$output_directory/mixed_elasticity_stress_1d.png"
test -s "$output_directory/mixed_elasticity_energy_1d.png"
test -s "$output_directory/mixed_elasticity_error_1d.png"
test -s "$output_directory/mixed_elasticity_solution_2d.csv"
test -s "$output_directory/benchmark.txt"
test -s "$output_directory/gallery_sequence.txt"
test "$(file -b "$output_directory/mixed_elasticity_solution_2d.png" | \
    cut -d' ' -f1)" = PNG

mapfile -t stages <"$output_directory/gallery_sequence.txt"
test "${#stages[@]}" -eq 2
test "${stages[0]}" = physical_solution
test "${stages[1]}" = diagnostics

primary_name=$(awk '$1 == "mixed_elasticity_wave" {print $2}' \
    "$repository_dir/doc/examples/primary_plots.txt")
test "$primary_name" = mixed_elasticity_solution_1d.png

python3 - "$output_directory" <<'PY'
import csv
import math
import pathlib
import sys

root = pathlib.Path(sys.argv[1])
benchmark = {}
with (root / "benchmark.txt").open(encoding="utf-8") as stream:
    next(stream)
    for line in stream:
        name, value = line.rstrip().split(",", 1)
        benchmark[name] = float(value)

assert abs(benchmark["final_time"] - 4.0) < 1.0e-14
expected_mode_1 = 0.7 * math.cos(math.pi * 4.0) - 0.15 * math.sin(math.pi * 4.0)
expected_mode_2 = -0.25 * math.cos(2.0 * math.pi * 4.0) - 0.35 * math.sin(2.0 * math.pi * 4.0)
assert abs(benchmark["final_displacement_mode_1"] - expected_mode_1) < 5.0e-3
assert abs(benchmark["final_displacement_mode_2"] - expected_mode_2) < 5.0e-3

with (root / "mixed_elasticity_solution_2d.csv").open(newline="", encoding="utf-8") as stream:
    rows = list(csv.DictReader(stream))
assert len(rows) == 9 * 9
d1 = benchmark["final_displacement_mode_1"]
d2 = benchmark["final_displacement_mode_2"]
for row in rows[::7]:
    x = float(row["x"])
    y = float(row["y"])
    ux = float(row["displacement_x"])
    uy = float(row["displacement_y"])
    stress = float(row["stress_magnitude"])
    expected_ux = d1 * math.sin(math.pi * x) * math.sin(math.pi * y)
    expected_uy = d2 * math.sin(2.0 * math.pi * x) * math.sin(2.0 * math.pi * y)
    expected_stress = math.hypot(
        math.pi * d1 * math.cos(math.pi * x) * math.sin(math.pi * y),
        2.0 * math.pi * d2 * math.sin(2.0 * math.pi * x) * math.cos(2.0 * math.pi * y),
    )
    assert abs(ux - expected_ux) < 2.0e-13
    assert abs(uy - expected_uy) < 2.0e-13
    assert abs(stress - expected_stress) < 2.0e-13

print("mixed elasticity gallery renders 2D physical fields before diagnostics")
PY
