#!/usr/bin/env bash
set -euo pipefail

# Independent output oracle for the physical vector open-boundary gallery.
# It does not reimplement the FEM/BEM or PML algebra; it checks that the
# executable emits finite, nonzero vector data and that physical plots come
# before diagnostics.
repository_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
output_directory="$repository_dir/output/example/maxwell_vector_open_boundary_parity_gallery"
log_file="${TMPDIR:-/tmp}/fortfem-vector-open-boundary-gallery.log"

(
    cd "$repository_dir"
    timeout --foreground "${EXAMPLE_TIMEOUT:-10}s" \
        env MAXWELL_VECTOR_OPEN_BOUNDARY_FAST=1 \
        fo exec --no-build --example maxwell_vector_open_boundary_parity_gallery \
        >"$log_file"
)

test -s "$output_directory/maxwell_vector_solution_2d.png"
test -s "$output_directory/maxwell_vector_pml_solution_2d.png"
test -s "$output_directory/maxwell_vector_torus_geometry_3d.png"
test -s "$output_directory/torus_vector_field.csv"
test -s "$output_directory/pml_vector_field.csv"
test -s "$output_directory/benchmark.txt"
test -s "$output_directory/gallery_sequence.txt"
test -s "$repository_dir/example/maxwell_vector_open_boundary_parity_gallery/provenance.json"
test "$(file -b "$output_directory/maxwell_vector_solution_2d.png" | cut -d' ' -f1)" = PNG
test "$(file -b "$output_directory/maxwell_vector_torus_geometry_3d.png" | cut -d' ' -f1)" = PNG

python3 - "$output_directory" "$repository_dir/example/maxwell_vector_open_boundary_parity_gallery/provenance.json" <<'PY'
import csv
import json
import math
import pathlib
import sys

root = pathlib.Path(sys.argv[1])
provenance = json.loads(pathlib.Path(sys.argv[2]).read_text(encoding="utf-8"))
assert provenance["fixture"] == "maxwell_vector_open_boundary_parity_gallery"
assert provenance["exact_reproduction"] is False
assert len(provenance["references"]) >= 1
for reference in provenance["references"]:
    assert reference["url"].startswith("https://")

with (root / "torus_vector_field.csv").open(newline="", encoding="utf-8") as stream:
    torus_rows = list(csv.DictReader(stream))
assert len(torus_rows) == 32 * 24
torus_norms = []
for row in torus_rows:
    values = [float(row[name]) for name in (
        "x", "z", "Ex_real", "Ey_real", "Ez_real", "E_magnitude")]
    assert all(math.isfinite(value) for value in values)
    expected = math.sqrt(sum(value * value for value in values[2:5]))
    assert abs(expected - values[5]) < 2.0e-12
    torus_norms.append(values[5])
assert max(torus_norms) > 1.0e-4
assert max(torus_norms) - min(torus_norms) > 1.0e-4

with (root / "pml_vector_field.csv").open(newline="", encoding="utf-8") as stream:
    pml_rows = list(csv.DictReader(stream))
assert len(pml_rows) == 24 * 24
pml_norms = []
for row in pml_rows:
    values = [float(row[name]) for name in (
        "x", "y", "Ex_real", "Ey_real", "Ez_real", "Ex_imag", "Ey_imag",
        "Ez_imag", "E_magnitude")]
    assert all(math.isfinite(value) for value in values)
    expected = math.sqrt(sum(value * value for value in values[2:8]))
    assert abs(expected - values[8]) < 2.0e-12
    pml_norms.append(values[8])
assert max(pml_norms) > 1.0e-4

assert (root / "gallery_sequence.txt").read_text().splitlines() == [
    "physical_solution", "diagnostics"]
benchmark = (root / "benchmark.txt").read_text(encoding="utf-8")
assert "dtn_response_norm" in benchmark
assert "larger_domain_difference" in benchmark
assert "fast_gallery,T" in benchmark
print("Vector Maxwell FEM--BEM/DtN/PML gallery emits finite nonzero physical fields")
PY
