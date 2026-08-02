#!/usr/bin/env bash
set -euo pipefail

repository_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
output_directory="$repository_dir/output/example/singular_layer_matching_gallery"
log_file="${TMPDIR:-/tmp}/fortfem-singular-layer-gallery.log"

rm -rf -- "$output_directory"
(
    cd "$repository_dir"
    timeout --foreground "${EXAMPLE_TIMEOUT:-30}s" \
        fpm run --example singular_layer_matching_gallery >"$log_file"
)

test -s "$output_directory/singular_layer_solution_1d.png"
test -s "$output_directory/singular_layer_diagnostics_1d.png"
test -s "$output_directory/singular_layer.csv"
test -s "$output_directory/provenance.json"
test -s "$output_directory/gallery_sequence.txt"
test "$(file -b "$output_directory/singular_layer_solution_1d.png" | cut -d' ' -f1)" = PNG

mapfile -t stages <"$output_directory/gallery_sequence.txt"
test "${#stages[@]}" -eq 2
test "${stages[0]}" = physical_solution
test "${stages[1]}" = diagnostics

primary_name=$(awk '$1 == "singular_layer_matching_gallery" {print $2}' \
    "$repository_dir/doc/examples/primary_plots.txt")
test "$primary_name" = singular_layer_solution_1d.png
test -s "$output_directory/$primary_name"

python3 - "$output_directory" <<'PY'
import csv
import json
import math
import pathlib
import sys

root = pathlib.Path(sys.argv[1])
with (root / "provenance.json").open(encoding="utf-8") as stream:
    record = json.load(stream)
assert record["schema"] == "fortfem-singular-layer-gallery-v1"
assert record["contract"] == "singular_layer_matching"
assert record["closure"] == "neutral-caller-owned-traces"
assert record["primary_plot"] == "singular_layer_solution_1d.png"
assert record["sample_count"] == 161
assert record["matched_residual_linf"] < 1.0e-13
assert record["mismatch_residual_linf"] > 1.0e-3
assert record["jvp_central_error"] < 2.0e-9

with (root / "singular_layer.csv").open(newline="", encoding="utf-8") as stream:
    rows = list(csv.DictReader(stream))
assert len(rows) == record["sample_count"]
for row in rows[::17]:
    s = float(row["coordinate"])
    outer_operator = 1.0 + 0.15 * s
    inner_operator = 0.9 - 0.10 * s
    outer = outer_operator * complex(0.80, 0.35)
    inner = inner_operator * complex(0.55, -0.20)
    jump = outer - inner
    mismatch = 0.025 * math.exp(-4.0 * s * s) * complex(1.0, -0.35)
    weight = 1.0 + 0.20 * math.cos(math.pi * s)
    assert abs(complex(float(row["outer_real"]), float(row["outer_imag"])) - outer) < 2.0e-12
    assert abs(complex(float(row["inner_real"]), float(row["inner_imag"])) - inner) < 2.0e-12
    assert abs(complex(float(row["jump_real"]), float(row["jump_imag"])) - jump) < 2.0e-12
    assert abs(complex(float(row["mismatch_real"]), float(row["mismatch_imag"])) - mismatch) < 2.0e-12
    assert abs(float(row["weight"]) - weight) < 2.0e-12
    residual = complex(float(row["residual_real"]), float(row["residual_imag"]))
    assert abs(residual) < 2.0e-13
    mismatch_residual = complex(float(row["mismatch_residual_real"]), float(row["mismatch_residual_imag"]))
    assert abs(mismatch_residual + weight * mismatch) < 2.0e-12
    weight_dot = 0.04 * math.sin(math.pi * s)
    mismatch_dot = 0.011 * math.sin(math.pi * s) * complex(1.0, -0.35)
    residual_jvp = complex(float(row["residual_jvp_real"]), float(row["residual_jvp_imag"]))
    assert abs(residual_jvp + weight_dot * mismatch + weight * mismatch_dot) < 2.0e-12
    assert float(row["jvp_error"]) < 2.0e-9

print("singular-layer gallery renders physical matching before diagnostics")
PY
