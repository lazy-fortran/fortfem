#!/usr/bin/env bash
set -euo pipefail

repository_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
output_directory="$repository_dir/output/example/maxwell_open_boundary_comparison"

(
    cd "$repository_dir"
    timeout --foreground 10s fo exec --no-build maxwell_open_boundary_comparison \
        >/tmp/fortfem-maxwell-open-boundary-gallery.log
)

test -s "$output_directory/maxwell_pml_field_slice_2d.png"
test -s "$output_directory/maxwell_pml_field_slice.csv"
test -s "$output_directory/maxwell_domain_comparison_1d.png"
test -s "$output_directory/maxwell_dtn_modes_1d.png"
test -s "$output_directory/maxwell_reflection_1d.png"
test -s "$output_directory/benchmark.txt"

primary_name=$(awk '$1 == "maxwell_open_boundary_comparison" {print $2}' \
    "$repository_dir/doc/examples/primary_plots.txt")
test "$primary_name" = maxwell_pml_field_slice_2d.png
test "$(file -b "$output_directory/$primary_name" | cut -d' ' -f1)" = PNG
test "$primary_name" != maxwell_reflection_1d.png

python3 - "$output_directory/maxwell_pml_field_slice.csv" \
    "$output_directory/benchmark.txt" <<'PY'
import csv
import math
import re
import sys
from collections import defaultdict

with open(sys.argv[1], newline="") as stream:
    rows = list(csv.DictReader(stream))
if len(rows) != 1600:
    raise SystemExit(f"unexpected Maxwell field sample count: {len(rows)}")
by_x = defaultdict(list)
components = defaultdict(list)
for row in rows:
    values = [float(row[name]) for name in
              ("x", "y", "Ex_real", "Ey_real", "Ez_real", "Ex_imag",
               "Ey_imag", "Ez_imag", "E_magnitude")]
    if not all(math.isfinite(value) for value in values):
        raise SystemExit("non-finite Maxwell PML field sample")
    ex, ey, ez, exi, eyi, ezi, magnitude = values[2:]
    if abs(math.sqrt(ex * ex + ey * ey + ez * ez + exi * exi + eyi * eyi +
                     ezi * ezi) - magnitude) > 1.0e-11:
        raise SystemExit("Maxwell magnitude is not vector-derived")
    by_x[values[0]].append(magnitude)
    for name, value in zip(("Ex", "Ey", "Ez"), (ex, ey, ez)):
        components[name].append(value)
if max(abs(value) for value in components["Ey"]) <= 1.0e-2:
    raise SystemExit("Maxwell tangential vector component is absent")
if max(abs(value) for value in components["Ex"]) <= 1.0e-4:
    raise SystemExit("Maxwell transverse vector component is absent")
first_x = min(by_x)
last_x = max(by_x)
first_mean = sum(by_x[first_x]) / len(by_x[first_x])
last_mean = sum(by_x[last_x]) / len(by_x[last_x])
if first_mean <= last_mean + 0.1:
    raise SystemExit("PML solution does not show expected attenuation")

bench = open(sys.argv[2], encoding="utf-8").read()

def number(label):
    match = re.search(re.escape(label) + r"\s*:\s*([0-9.Ee+-]+)", bench)
    if match is None:
        raise SystemExit(f"missing benchmark field: {label}")
    value = float(match.group(1))
    if not math.isfinite(value):
        raise SystemExit(f"non-finite benchmark field: {label}")
    return value

if number("maximum TE/TM eigenvalue error") > 1.0e-10:
    raise SystemExit("Maxwell DtN modal error is too large")
if number("Nedelec PML relative edge error") > 2.0e-2:
    raise SystemExit("Maxwell PML edge error is too large")
if number("larger-domain PML field difference") > 1.0e-1:
    raise SystemExit("larger-domain PML comparison is inconsistent")
if number("Nedelec PML solve seconds") > 10.0:
    raise SystemExit("Maxwell PML solve exceeded ten-second runtime")
PY

echo "Maxwell gallery emits an attenuating vector PML solution before diagnostics"
