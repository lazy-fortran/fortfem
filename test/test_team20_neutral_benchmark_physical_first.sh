#!/usr/bin/env bash
set -euo pipefail

# Independent behavioral gate for the solution-first TEAM-20-shaped gallery.
repository_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
output_directory="$repository_dir/output/example/team20_neutral_benchmark"

(
    cd "$repository_dir"
    timeout --foreground 10s fo exec --no-build --example team20_neutral_benchmark \
        >/tmp/fortfem-team20-neutral-gallery.log
)

test -s "$output_directory/solution.png"
test -s "$output_directory/solution_3d.png"
test -s "$output_directory/probe.png"
test -s "$output_directory/solution.csv"
test -s "$output_directory/diagnostics.csv"

python3 - "$output_directory/solution.csv" "$output_directory/diagnostics.csv" <<'PY'
import math
import pathlib
import sys

solution = pathlib.Path(sys.argv[1])
diagnostics = pathlib.Path(sys.argv[2])
rows = []
for line in solution.read_text().splitlines()[1:]:
    fields = [float(value) for value in line.split(",")]
    assert len(fields) == 7
    assert all(math.isfinite(value) for value in fields)
    rows.append(fields)
assert len(rows) > 1000
expected_energy = sum(0.5 * (row[3] ** 2 + row[4] ** 2 + row[5] ** 2) for row in rows)
expected_energy /= len(rows)
assert expected_energy > 0.0
values = {}
for line in diagnostics.read_text().splitlines():
    name, value = line.split(",", 1)
    values[name] = float(value)
assert values["mean_magnetic_energy"] > 0.0
assert abs(values["mean_magnetic_energy"] - expected_energy) < 2.0e-12
assert values["divergence_proxy"] < 2.0e-4
PY

echo "TEAM-20 neutral gallery emits a physical 2-D/3-D solution before diagnostics"
