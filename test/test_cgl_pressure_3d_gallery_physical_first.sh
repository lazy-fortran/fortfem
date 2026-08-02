#!/usr/bin/env bash
set -euo pipefail

repository_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
output_directory="$repository_dir/output/example/cgl_pressure_3d_gallery"
log_file="${TMPDIR:-/tmp}/fortfem-cgl-pressure-3d-gallery.log"

rm -rf -- "$output_directory"
(
    cd "$repository_dir"
    timeout --foreground "${EXAMPLE_TIMEOUT:-10}s" \
        fpm run --example cgl_pressure_3d_gallery >"$log_file"
)

test -s "$output_directory/solution.png"
test -s "$output_directory/solution_3d.png"
test -s "$output_directory/pressure_components.png"
test -s "$output_directory/solution.csv"
test -s "$output_directory/diagnostics.csv"
test -s "$output_directory/gallery_sequence.txt"
test "$(file -b "$output_directory/solution.png" | cut -d' ' -f1)" = PNG

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
with (root / "solution.csv").open(newline="", encoding="utf-8") as stream:
    rows = list(csv.DictReader(stream))
assert len(rows) == 4**3

for row in rows:
    x = float(row["x"])
    y = float(row["y"])
    z = float(row["z"])
    raw_b = (
        1.0 + 0.25 * math.sin(math.pi * x),
        2.0 + 0.20 * math.cos(math.pi * y),
        2.0 + 0.15 * math.sin(math.pi * z),
    )
    bnorm = math.sqrt(sum(value * value for value in raw_b))
    b = tuple(value / bnorm for value in raw_b)
    ppar = 3.6 + 0.4 * x - 0.25 * y + 0.18 * z
    pperp = 1.2 + 0.15 * math.cos(math.pi * x) * math.cos(0.5 * math.pi * y) * math.cos(0.75 * math.pi * z)
    delta = ppar - pperp
    expected = [
        pperp + delta * b[0] * b[0],
        delta * b[0] * b[1],
        delta * b[0] * b[2],
        pperp + delta * b[1] * b[1],
        delta * b[1] * b[2],
        pperp + delta * b[2] * b[2],
    ]
    actual = [float(row[name]) for name in ("P11", "P12", "P13", "P22", "P23", "P33")]
    assert max(abs(left - right) for left, right in zip(actual, expected)) < 3.0e-12
    assert max(abs(float(row[name]) - value) for name, value in zip(("bx", "by", "bz"), b)) < 3.0e-12
    normal = [float(row[name]) for name in ("nx", "ny", "nz")]
    traction = [float(row[name]) for name in ("t1", "t2", "t3")]
    expected_traction = [
        expected[0] * normal[0] + expected[1] * normal[1] + expected[2] * normal[2],
        expected[1] * normal[0] + expected[3] * normal[1] + expected[4] * normal[2],
        expected[2] * normal[0] + expected[4] * normal[1] + expected[5] * normal[2],
    ]
    assert max(abs(left - right) for left, right in zip(traction, expected_traction)) < 3.0e-12
    trace = actual[0] + actual[3] + actual[5]
    assert trace > 0.0

with (root / "diagnostics.csv").open(newline="", encoding="utf-8") as stream:
    diagnostics = {row["quantity"]: float(row["value"]) for row in csv.DictReader(stream)}
for name in diagnostics:
    if name.endswith("error"):
        assert diagnostics[name] < 3.0e-10, (name, diagnostics[name])
assert diagnostics["minimum_pressure_trace"] > 0.0
assert diagnostics["kernel_seconds"] < 10.0

print("3-D CGL gallery has independently reconstructed tensor and traction fields")
PY
