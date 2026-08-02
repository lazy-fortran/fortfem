#!/usr/bin/env bash
set -euo pipefail

repository_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
output_directory="$repository_dir/output/example/iga_jorek_flux"

(
    cd "$repository_dir"
    timeout --foreground 10s fpm run --example iga_jorek_flux \
        >/tmp/fortfem-iga-jorek-flux-gallery.log
)

test -s "$output_directory/jorek_flux_initial_2d.png"
test -s "$output_directory/jorek_flux_initial_field.csv"
test -s "$output_directory/jorek_flux_final_2d.png"
test -s "$output_directory/jorek_flux_invariant_1d.png"
test -s "$output_directory/benchmark.txt"

primary_name=$(awk '$1 == "iga_jorek_flux" {print $2}' \
    "$repository_dir/doc/examples/primary_plots.txt")
test "$primary_name" = jorek_flux_initial_2d.png
test "$(file -b "$output_directory/$primary_name" | cut -d' ' -f1)" = PNG
test "$primary_name" != jorek_flux_invariant_1d.png

python3 - "$output_directory/jorek_flux_initial_field.csv" \
    "$output_directory/benchmark.txt" <<'PY'
import csv
import math
import re
import sys

with open(sys.argv[1], newline="") as stream:
    rows = list(csv.DictReader(stream))
if len(rows) != 1600:
    raise SystemExit(f"unexpected JOREK field sample count: {len(rows)}")
for row in rows:
    values = [float(row[name]) for name in
              ("R", "Z", "psi", "Br", "Bz", "B_magnitude")]
    if not all(math.isfinite(value) for value in values):
        raise SystemExit("non-finite JOREK field sample")
    br, bz, magnitude = values[3:]
    if abs(math.hypot(br, bz) - magnitude) > 1.0e-11:
        raise SystemExit("JOREK field magnitude is not vector-derived")
psi = [float(row["psi"]) for row in rows]
br = [float(row["Br"]) for row in rows]
bz = [float(row["Bz"]) for row in rows]
if min(psi) >= -1.0e-2 or max(psi) <= 1.0e-2:
    raise SystemExit("JOREK flux solution has no signed structure")
if max(abs(value) for value in br) <= 1.0e-6:
    raise SystemExit("JOREK radial field component is absent")
if max(abs(value) for value in bz) <= 1.0e-6:
    raise SystemExit("JOREK vertical field component is absent")

bench = open(sys.argv[2], encoding="utf-8").read()

def number(label):
    match = re.search(re.escape(label) + r"\s*:\s*([0-9.Ee+-]+)", bench)
    if match is None:
        raise SystemExit(f"missing benchmark field: {label}")
    value = float(match.group(1))
    if not math.isfinite(value):
        raise SystemExit(f"non-finite benchmark field: {label}")
    return value

if number("maximum relative mass-norm drift") > 5.0e-11:
    raise SystemExit("JOREK midpoint invariant drift is too large")
if number("forward-backward coefficient error") > 5.0e-11:
    raise SystemExit("JOREK midpoint reversibility error is too large")
if number("forward propagation seconds") > 10.0:
    raise SystemExit("JOREK gallery exceeded ten-second runtime")
PY

echo "JOREK gallery emits a signed flux solution and reconstructed vector field first"
