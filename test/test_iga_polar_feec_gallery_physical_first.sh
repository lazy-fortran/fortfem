#!/usr/bin/env bash
set -euo pipefail

repository_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
output_directory="$repository_dir/output/example/iga_polar_feec"
log_file="${TMPDIR:-/tmp}/fortfem-iga-polar-feec-gallery.log"

rm -rf -- "$output_directory"
(
    cd "$repository_dir"
    timeout --foreground "${EXAMPLE_TIMEOUT:-10}s" \
        fpm run --example iga_polar_feec >"$log_file"
)

test -s "$output_directory/polar_feec_solution_2d.png"
test -s "$output_directory/polar_curvilinear_mesh_2d.png"
test -s "$output_directory/polar_curvilinear_mesh_2d.csv"
test -s "$output_directory/polar_feec_energy_1d.png"
test -s "$output_directory/polar_feec_solution_2d.csv"
test -s "$output_directory/gallery_sequence.txt"
test "$(file -b "$output_directory/polar_feec_solution_2d.png" | \
    cut -d' ' -f1)" = PNG

mapfile -t stages <"$output_directory/gallery_sequence.txt"
test "${#stages[@]}" -eq 2
test "${stages[0]}" = physical_solution
test "${stages[1]}" = diagnostics

primary_name=$(awk '$1 == "iga_polar_feec" {print $2}' \
    "$repository_dir/doc/examples/primary_plots.txt")
test "$primary_name" = polar_feec_solution_2d.png

python3 - "$output_directory" <<'PY'
import csv
import math
import pathlib
import sys

root = pathlib.Path(sys.argv[1])
with (root / "polar_feec_solution_2d.csv").open(
    newline="", encoding="utf-8"
) as stream:
    rows = list(csv.DictReader(stream))

# The oracle checks the physical map rather than reproducing the Fortran
# evaluator: each mapped radial line remains a circle of its parameter radius,
# while the angular samples visit all four quadrants and the magnetic axis.
assert len(rows) == 49 * 128
for row in rows[::31]:
    radius = float(row["radius_parameter"])
    x = float(row["x"])
    y = float(row["y"])
    physical_radius = math.hypot(x, y)
    if radius > 1.0e-12:
        assert abs(physical_radius / radius - 1.0) < 2.0e-3
    else:
        assert physical_radius < 1.0e-12

outer = rows[-128:]
assert abs(float(outer[0]["x"]) - 1.0) < 1.0e-12
quarter = outer[32]
assert abs(float(quarter["x"])) < 2.0e-3
assert abs(float(quarter["y"]) - 1.0) < 2.0e-3
assert max(float(row["solution"]) for row in rows) > 0.1
assert min(float(row["solution"]) for row in rows) < -0.1

# Independent geometry oracle for the line paths sent to FortPlot.  The
# renderer receives one complete 121-point path for every physical ring and
# spoke; checking those paths catches truncated/parameter-space mesh plots
# without depending on renderer-specific PNG pixels.
with (root / "polar_curvilinear_mesh_2d.csv").open(
    newline="", encoding="utf-8"
) as stream:
    mesh_rows = list(csv.DictReader(stream))

expected_lines = 4 + 64
samples_per_line = 121
assert len(mesh_rows) == expected_lines * samples_per_line
for kind, line_count in (("ring", 4), ("spoke", 64)):
    selected = [row for row in mesh_rows if row["line_kind"] == kind]
    assert len(selected) == line_count * samples_per_line
    line_indices = range(2, 6) if kind == "ring" else range(1, 65)
    for line in line_indices:
        path = [row for row in selected if int(row["line_index"]) == line]
        assert len(path) == samples_per_line
        assert [int(row["sample"]) for row in path] == list(
            range(1, samples_per_line + 1)
        )
        assert all(math.isfinite(float(row[name])) for row in path
                   for name in ("parameter", "x", "y"))
        if kind == "ring":
            radius = [math.hypot(float(row["x"]), float(row["y"]))
                      for row in path]
            expected_radius = (line - 1) / 4.0
            assert max(abs(value - expected_radius) for value in radius) < 2.0e-3
            assert abs(float(path[0]["x"]) - float(path[-1]["x"])) < 2.0e-12
            assert abs(float(path[0]["y"]) - float(path[-1]["y"])) < 2.0e-12
        else:
            radii = [math.hypot(float(row["x"]), float(row["y"]))
                     for row in path]
            assert abs(radii[0]) < 2.0e-12
            assert abs(radii[-1] - 1.0) < 2.0e-12
            assert all(radii[index] <= radii[index + 1] + 2.0e-12
                       for index in range(samples_per_line - 1))
PY

echo "polar FEEC gallery renders the physical curvilinear mesh before diagnostics"
