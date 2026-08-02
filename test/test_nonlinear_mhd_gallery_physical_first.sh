#!/usr/bin/env bash
set -euo pipefail

repository_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
output_directory="$repository_dir/output/example/nonlinear_resistive_mhd_gallery"

rm -rf "$output_directory"
timeout --foreground 30s fpm run --example nonlinear_resistive_mhd_gallery \
    >/tmp/fortfem-nonlinear-mhd-gallery.log

test -s "$output_directory/gallery_sequence.txt"
mapfile -t stages <"$output_directory/gallery_sequence.txt"
test "${stages[0]}" = "physical_solution"
test "${stages[1]}" = "branch_diagnostics"
test "${stages[2]}" = "ledger_diagnostics"

primary_name=$(awk '$1 == "nonlinear_resistive_mhd_gallery" {print $2}' \
    "$repository_dir/doc/examples/primary_plots.txt")
test "$primary_name" = "island_wall_solution_1d.png"
test -s "$output_directory/$primary_name"
test -s "$output_directory/island_wall_continuation.csv"
test -s "$output_directory/benchmark.json"

python3 - "$output_directory/benchmark.json" <<'PY'
import json
import math
import sys

with open(sys.argv[1], encoding="utf-8") as stream:
    record = json.load(stream)
assert record["branch_multiplicity"] == 2
assert record["hysteresis_detected"] is True
assert record["final_residual_norm"] <= 1.0e-13
assert record["final_dissipation"] >= 0.0
assert all(math.isfinite(record[name]) for name in (
    "max_state_hysteresis", "path_metric", "final_input_power",
    "final_dissipation", "composition_seconds_per_call"))
PY

echo "nonlinear island/wall gallery renders physical solution before diagnostics"
