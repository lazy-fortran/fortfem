#!/usr/bin/env bash
set -euo pipefail

repository_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
output_directory="$repository_dir/output/example/tetra_nedelec_p_convergence"

rm -rf -- "$output_directory"
(
    cd "$repository_dir"
    timeout --foreground 10s fpm run --example tetra_nedelec_p_convergence \
        >/tmp/fortfem-tetra-nedelec-gallery.log
)

test -s "$output_directory/nedelec_field_3d.png"
test -s "$output_directory/nedelec_field_slice_2d.png"
test -s "$output_directory/nedelec_field_slice_2d.csv"
test -s "$output_directory/p_convergence_1d.png"
test -s "$output_directory/gallery_sequence.txt"

mapfile -t stages <"$output_directory/gallery_sequence.txt"
test "${#stages[@]}" -eq 2
test "${stages[0]}" = physical_solution
test "${stages[1]}" = diagnostics

primary_name=$(awk '$1 == "tetra_nedelec_p_convergence" {print $2}' \
    "$repository_dir/doc/examples/primary_plots.txt")
test "$primary_name" = nedelec_field_slice_2d.png

# The primary artifact must be a rendered image and must not be the diagnostic
# convergence curve that is emitted last by the example.
test "$(file -b "$output_directory/$primary_name" | cut -d' ' -f1)" = PNG
test "$primary_name" != p_convergence_1d.png

# Check the physical data independently of the rendered PNG.  A scalar-only
# plot can pass the image checks while silently dropping the vector components.
python3 - "$output_directory/nedelec_field_slice_2d.csv" <<'PY'
import csv
import math
import sys

with open(sys.argv[1], newline="") as stream:
    rows = list(csv.DictReader(stream))
if len(rows) < 100:
    raise SystemExit("Nedelec physical slice has too few samples")
for name in ("Ex", "Ey", "Ez", "magnitude"):
    values = [float(row[name]) for row in rows]
    if not all(math.isfinite(value) for value in values):
        raise SystemExit(f"non-finite Nedelec component: {name}")
    if max(abs(value) for value in values) <= 1.0e-12:
        raise SystemExit(f"Nedelec component is visually absent: {name}")
for row in rows:
    field = [float(row[name]) for name in ("Ex", "Ey", "Ez")]
    magnitude = float(row["magnitude"])
    if abs(math.sqrt(sum(value * value for value in field)) - magnitude) > 1.0e-11:
        raise SystemExit("Nedelec magnitude does not match vector components")
PY

echo "tetrahedral Nedelec gallery emits vector solution before diagnostics"
