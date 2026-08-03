#!/usr/bin/env bash
set -euo pipefail

# Independent output oracle for the solution-first Poisson/Ampere examples.
repository_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
fixture_dir=$(mktemp -d /mnt/storage/fortfem-poisson-ampere.XXXXXX)
cleanup() { find "$fixture_dir" -depth -delete; }
trap cleanup EXIT

(
    cd "$repository_dir"
    FO_JOBS=1 timeout --foreground 30s fo build >"$fixture_dir/build.log"
    timeout --foreground 30s fo exec --no-build mixed_poisson \
        >"$fixture_dir/mixed_poisson.log"
    timeout --foreground 30s fo exec --no-build sheet_current_surface_gallery \
        >"$fixture_dir/sheet_current_surface_gallery.log"
)

python3 - "$repository_dir/output/example/mixed_poisson" \
    "$repository_dir/output/example/sheet_current_surface_gallery" <<'PY'
import csv
import json
import math
import pathlib
import sys

poisson_dir, ampere_dir = map(pathlib.Path, sys.argv[1:])

poisson_stages = (poisson_dir / "gallery_sequence.txt").read_text().splitlines()
assert poisson_stages[0] == "physical_solution"
assert poisson_stages[-1] == "diagnostics"
poisson = json.loads((poisson_dir / "benchmark.json").read_text())
assert poisson["schema"] == "fortfem-mixed-poisson-benchmark-v1"
assert poisson["physical_solution_first"] is True
assert poisson["mesh_levels"] == 3
assert poisson["finest_vertex_count"] == 17
for key in ("solve_seconds", "total_seconds"):
    assert math.isfinite(poisson[key]) and poisson[key] > 0.0
assert poisson["total_seconds"] >= poisson["solve_seconds"]
for key in ("finest_flux_l2_error", "finest_pressure_l2_error"):
    assert math.isfinite(poisson[key]) and poisson[key] > 0.0

ampere_stages = (ampere_dir / "gallery_sequence.txt").read_text().splitlines()
assert ampere_stages[:4] == [
    "physical_solution_slab",
    "physical_solution_cylinder",
    "physical_solution_sphere",
    "physical_solution_torus",
]
assert ampere_stages[-1] == "diagnostics"
ampere = json.loads((ampere_dir / "benchmark.json").read_text())
assert ampere["schema"] == "fortfem-sheet-current-surface-gallery-v2"
assert ampere["physical_solution_first"] is True
assert ampere["total_surface_points"] > 300
assert math.isfinite(ampere["elapsed_seconds"])
assert ampere["elapsed_seconds"] > 0.0

ledger = list(csv.DictReader((ampere_dir / "geometry_ledger.csv").open(newline="")))
assert [row["geometry"] for row in ledger] == ["slab", "cylinder", "sphere", "torus"]
errors = []
for row in ledger:
    fitted = [float(row[f"fitted_{axis}"]) for axis in "xyz"]
    resolved = [float(row[f"resolved_{axis}"]) for axis in "xyz"]
    scale = max(1.0, math.sqrt(sum(value * value for value in fitted)))
    error = math.sqrt(sum((resolved[i] - fitted[i]) ** 2 for i in range(3))) / scale
    recorded = float(row["relative_error"])
    assert math.isfinite(error) and error >= 0.0
    assert math.isclose(recorded, error, rel_tol=2.0e-13, abs_tol=2.0e-15)
    errors.append(error)
maximum_error = max(errors)
assert maximum_error < 4.0e-12
assert math.isclose(
    ampere["max_integrated_ampere_relative_error"],
    maximum_error,
    rel_tol=2.0e-13,
    abs_tol=2.0e-15,
)
PY

echo "Poisson and Ampere galleries emit solution-first timing/error metadata"
