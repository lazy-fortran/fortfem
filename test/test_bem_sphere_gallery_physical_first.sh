#!/usr/bin/env bash
set -euo pipefail

repository_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
output_directory="$repository_dir/output/example/bem_sphere_3d"
log_file="/mnt/storage/fortfem-bem-sphere-gallery.log"

(
    cd "$repository_dir"
    timeout --foreground "${EXAMPLE_TIMEOUT:-10}s" \
        fpm run --example bem_sphere_3d >"$log_file"
)

test -s "$output_directory/sphere_exterior_solution_3d.png"
test -s "$output_directory/sphere_exterior_solution_3d.csv"
test -s "$output_directory/gallery_sequence.txt"
test -s "$output_directory/benchmark.txt"
test "$(file -b "$output_directory/sphere_exterior_solution_3d.png" | cut -d' ' -f1)" = PNG
test "$(sed -n '1p' "$output_directory/gallery_sequence.txt")" = physical_solution

primary_name=$(awk '$1 == "bem_sphere_3d" {print $2}' \
    "$repository_dir/doc/examples/primary_plots.txt")
test "$primary_name" = sphere_exterior_solution_3d.png
test "$primary_name" != sphere_capacitance.png

python3 - "$output_directory/sphere_exterior_solution_3d.csv" <<'PY'
import csv
import math
import sys

with open(sys.argv[1], newline="", encoding="utf-8") as stream:
    rows = list(csv.DictReader(stream))

if len(rows) != 13 * 25:
    raise SystemExit(f"unexpected sphere surface sample count: {len(rows)}")
for row in rows:
    values = [float(row[name]) for name in
              ("theta", "phi", "x", "y", "z", "computed", "exact", "absolute_error")]
    if not all(math.isfinite(value) for value in values):
        raise SystemExit("non-finite BEM sphere surface sample")
    theta, phi, x, y, z, computed, exact, error = values
    radius = math.sqrt(x*x + y*y + z*z)
    if radius <= 1.0:
        raise SystemExit("observation surface is not exterior to the sphere")
    if abs(exact - 1.0/radius) > 1.0e-12:
        raise SystemExit("sphere CSV exact column is not the independent 1/r oracle")
    if abs(error - abs(computed - exact)) > 1.0e-12:
        raise SystemExit("sphere CSV absolute error is not value-derived")

# The CSV is a connected tensor-product surface, not an unordered point cloud:
# each azimuth block has monotone theta and the final periodic seam closes.
for offset in range(0, len(rows), 13):
    block = rows[offset:offset + 13]
    theta = [float(row["theta"]) for row in block]
    if any(b <= a for a, b in zip(theta, theta[1:])):
        raise SystemExit("sphere theta samples are not ordered in each surface strip")
first = rows[:13]
last = rows[-13:]
for left, right in zip(first, last):
    if max(abs(float(left[name]) - float(right[name])) for name in ("x", "y", "z")) > 1.0e-12:
        raise SystemExit("sphere surface azimuth seam is not connected")
if max(float(row["absolute_error"]) for row in rows) >= 4.0e-2:
    raise SystemExit("BEM sphere field does not reproduce the exterior monopole")
PY

echo "3-D exterior BEM gallery renders a connected physical sphere solution first"
