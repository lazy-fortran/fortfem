#!/usr/bin/env bash
set -euo pipefail

# Independent fast oracle for the physical scattering gallery.  The reference
# three-by-three torus run is intentionally more expensive; this gate executes
# the same CFIE/DtN/field-reconstruction path with reduced quadrature and keeps
# the ten-second test budget useful for docs builds.
repository_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
output_directory="$repository_dir/output/example/maxwell_torus_curved_scattering"
log_file="${TMPDIR:-/tmp}/fortfem-maxwell-torus-vector-gallery.log"

(
    cd "$repository_dir"
    timeout --foreground "${EXAMPLE_TIMEOUT:-10}s" env MAXWELL_TORUS_FAST=1 \
        fo exec --no-build --example maxwell_torus_curved_scattering >"$log_file"
)

test -s "$output_directory/maxwell_torus_solution_2d.png"
test -s "$output_directory/maxwell_torus_rcs_3d.png"
test -s "$output_directory/scattered_field.csv"
provenance_file="$repository_dir/example/maxwell_torus_curved_scattering/provenance.json"
test -s "$provenance_file"
test "$(file -b "$output_directory/maxwell_torus_solution_2d.png" | \
    cut -d' ' -f1)" = PNG

python3 - "$output_directory" "$provenance_file" <<'PY'
import csv
import json
import math
import pathlib
import sys

root = pathlib.Path(sys.argv[1])
record = json.loads(pathlib.Path(sys.argv[2]).read_text(encoding="utf-8"))
assert record["fixture"] == "maxwell_torus_curved_scattering"
assert record["exact_reproduction"] is False
assert "RWG" in record["discretization"]

with (root / "scattered_field.csv").open(newline="", encoding="utf-8") as stream:
    rows = list(csv.DictReader(stream))
assert len(rows) == 20 * 20
norms = []
for row in rows:
    values = [float(row[name]) for name in (
        "x", "z", "Hx_real", "Hy_real", "Hz_real", "H_scattered_magnitude")
    ]
    assert all(math.isfinite(value) for value in values)
    hx, hy, hz, magnitude = values[2:]
    real_norm = math.sqrt(hx * hx + hy * hy + hz * hz)
    # The magnitude is computed from the full complex field, so it must bound
    # the instantaneous real-vector norm independently of the Fortran code.
    assert magnitude + 2.0e-13 >= real_norm
    norms.append((real_norm, magnitude))
assert max(value for _, value in norms) > 1.0e-4
assert max(value for value, _ in norms) > 1.0e-4
assert max(value for _, value in norms) - min(value for _, value in norms) > 1.0e-4

benchmark = (root / "benchmark.txt").read_text(encoding="utf-8")
assert "CFIE quadrature degree: 2" in benchmark
assert "Lorentz reciprocity relative error:" in benchmark
print("Maxwell torus gallery emits a nonzero reconstructed vector solution and 3-D radiation surface")
PY
