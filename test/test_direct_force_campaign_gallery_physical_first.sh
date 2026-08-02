#!/usr/bin/env bash
set -euo pipefail

repo_root=$(cd "$(dirname "$0")/.." && pwd)
cd "$repo_root"
fo exec direct_force_campaign_gallery >/tmp/fortfem_direct_force_gallery.log

output_dir=output/example/direct_force_campaign_gallery
test -s "$output_dir/direct_force_torus_solution_3d.png"
test -s "$output_dir/direct_force_state_2d.png"
test -s "$output_dir/direct_force_objective_diagnostics_1d.png"
test -s "$output_dir/direct_force_surface.csv"
test -s "$output_dir/direct_force_diagnostics.csv"
test -s "$output_dir/benchmark.json"
test "$(head -n 1 "$output_dir/gallery_sequence.txt")" = physical_solution

python3 - "$output_dir" <<'PY'
import csv
import json
import math
import pathlib
import sys

out = pathlib.Path(sys.argv[1])
rows = list(csv.DictReader((out / "direct_force_surface.csv").open()))
assert len(rows) == 25 * 33
for row in rows:
    x, y, z = (float(row[k]) for k in ("x", "y", "z"))
    theta = float(row["theta"])
    radial = math.hypot(x, y)
    expected = 2.4 + 0.7 * math.cos(theta)
    assert abs(radial - expected) < 3e-12
    assert all(math.isfinite(float(row[k])) for k in row if k not in ("theta", "phi"))
report = json.loads((out / "benchmark.json").read_text())
assert report["samples"] == 25 * 33
assert report["fd_error"] < 2e-9
assert report["geometry_error"] < 2e-12
objective = 0.0
for row in rows:
    residual = float(row["residual"])
    weight = float(row["weight"])
    assert weight > 0.0
    objective += 0.5 * weight * residual * residual
assert abs(objective - report["objective"]) <= 2e-12
assert report["provenance"] == "analytic-torus-direct-force-v1"
assert report["primary_plot"] == "direct_force_torus_solution_3d.png"
assert report["closure"] == "neutral-caller-owned-residual"
print("direct-force gallery renders physical torus before diagnostics")
PY
