#!/usr/bin/env bash
set -euo pipefail

repository_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
output_directory="$repository_dir/output/example/cgl_pressure_tensor"
log_file="${TMPDIR:-/tmp}/fortfem-cgl-pressure-tensor-gallery.log"

rm -rf -- "$output_directory"
(
    cd "$repository_dir"
    timeout --foreground "${EXAMPLE_TIMEOUT:-30}s" \
        fpm run --example cgl_pressure_tensor >"$log_file"
)

test -s "$output_directory/cgl_pressure_components_1d.png"
test -s "$output_directory/cgl_force_divergence_1d.png"
test -s "$output_directory/cgl_tensor_principal_2d.png"
test -s "$output_directory/cgl_profile.csv"
test -s "$output_directory/cgl_tensor_field_2d.csv"
test -s "$output_directory/gallery_sequence.txt"
test "$(file -b "$output_directory/cgl_tensor_principal_2d.png" | cut -d' ' -f1)" = PNG

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
with (root / "cgl_tensor_field_2d.csv").open(newline="", encoding="utf-8") as stream:
    rows = list(csv.DictReader(stream))

assert len(rows) == 25 * 21
for row in rows[::37]:
    x = float(row["x"])
    y = float(row["y"])
    angle = 0.35 * x + 0.22 * y
    p_parallel = 3.0 + 0.2 * x - 0.12 * y
    p_perpendicular = 1.0 + 0.1 * math.cos(math.pi * x) * math.cos(0.5 * math.pi * y)
    bx = math.cos(angle)
    by = math.sin(angle)
    expected = (
        p_perpendicular + (p_parallel - p_perpendicular) * bx * bx,
        p_perpendicular + (p_parallel - p_perpendicular) * by * by,
        (p_parallel - p_perpendicular) * bx * by,
    )
    for name, value in zip(("P11", "P22", "P12"), expected):
        assert abs(float(row[name]) - value) < 3.0e-13
    assert abs(float(row["principal_x"]) - bx) < 3.0e-13
    assert abs(float(row["principal_y"]) - by) < 3.0e-13
    assert abs(float(row["trace_error"])) < 3.0e-13

print("CGL gallery renders the physical 2D tensor field before diagnostics")
PY
