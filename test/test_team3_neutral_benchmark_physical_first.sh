#!/usr/bin/env bash
set -euo pipefail

# Independent output oracle for the solution-first TEAM-3-shaped gallery.
repository_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
output_directory="$repository_dir/output/example/team3_neutral_benchmark"

(
    cd "$repository_dir"
    timeout --foreground 10s fo exec --no-build team3_neutral_benchmark \
        >/tmp/fortfem-team3-neutral-gallery.log
)

test -s "$output_directory/solution.png"
test -s "$output_directory/solution_3d.png"
test -s "$output_directory/probe.png"
test -s "$output_directory/solution.csv"
test -s "$output_directory/material.csv"
test -s "$output_directory/geometry.csv"
test -s "$output_directory/diagnostics.csv"
test -s "$output_directory/gallery_sequence.txt"
test -s "$repository_dir/example/team3_neutral_benchmark/provenance.json"

python3 - "$output_directory/solution.csv" "$output_directory/diagnostics.csv" \
    "$output_directory/gallery_sequence.txt" \
    "$repository_dir/example/team3_neutral_benchmark/provenance.json" <<'PY'
import math
import json
import pathlib
import sys

solution = pathlib.Path(sys.argv[1])
diagnostics = pathlib.Path(sys.argv[2])
sequence = pathlib.Path(sys.argv[3]).read_text().splitlines()
provenance = json.loads(pathlib.Path(sys.argv[4]).read_text())
assert provenance["fixture"] == "team3_neutral_benchmark"
assert provenance["exact_reproduction"] is False
assert provenance["source_arrays_redistributed"] is False
assert sequence.index("physical_solution") < sequence.index("diagnostics")
rows = []
for line in solution.read_text().splitlines()[1:]:
    values = [float(value) for value in line.split(",")]
    assert len(values) == 7
    assert all(math.isfinite(value) for value in values)
    rows.append(values)
assert len(rows) == 81 * 61
expected_energy = sum(0.5 * (row[3] ** 2 + row[4] ** 2) for row in rows)
expected_energy /= len(rows)
values = {}
for line in diagnostics.read_text().splitlines():
    name, value = line.split(",", 1)
    values[name] = float(value)
assert values["mean_magnetic_energy"] > expected_energy
assert values["mean_current_energy"] > 0.0
assert values["mean_source_energy"] > 0.0
assert values["max_relative_permeability"] > 1.0
assert values["divergence_proxy"] < 2.0e-4
PY

echo "TEAM-3 neutral gallery emits a physical 2-D/3-D solution before diagnostics"
