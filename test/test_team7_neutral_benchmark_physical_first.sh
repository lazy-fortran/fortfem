#!/usr/bin/env bash
set -euo pipefail

repository_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
output_directory="$repository_dir/output/example/team7_neutral_benchmark"

(
    cd "$repository_dir"
    timeout --foreground 10s fpm run --example team7_neutral_benchmark \
        >/tmp/fortfem-team7-neutral-gallery.log
)

test -s "$output_directory/solution.png"
test -s "$output_directory/current.png"
test -s "$output_directory/probe.png"
test -s "$output_directory/solution.csv"
test -s "$output_directory/diagnostics.csv"

python3 - "$output_directory/diagnostics.csv" <<'PY'
import pathlib
import sys

values = {}
for line in pathlib.Path(sys.argv[1]).read_text().splitlines():
    name, value = line.split(",", 1)
    values[name] = float(value)
assert values["mean_magnetic_energy"] > 0.0
assert values["mean_current_energy"] > 0.0
assert values["divergence_proxy"] < 2.0e-5
PY

echo "TEAM-7 neutral gallery emits solution before diagnostics"
