#!/usr/bin/env bash
set -euo pipefail

repository_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
output_directory="$repository_dir/output/example/field_aligned_anisotropic_3d_gallery"
rm -rf -- "$output_directory"
(
    cd "$repository_dir"
    timeout --foreground 10s fo exec --no-build field_aligned_anisotropic_3d_gallery \
        >/tmp/fortfem-field-aligned-anisotropic-3d.log
)

test -s "$output_directory/solution.png"
test -s "$output_directory/solution_3d.png"
test -s "$output_directory/solution_slice.png"
test -s "$output_directory/solution.csv"
test -s "$output_directory/diagnostics.csv"
test -s "$output_directory/gallery_sequence.txt"

mapfile -t stages <"$output_directory/gallery_sequence.txt"
test "${#stages[@]}" -eq 2
test "${stages[0]}" = physical_solution
test "${stages[1]}" = diagnostics

python3 - "$output_directory/diagnostics.csv" <<'PY'
import pathlib
import sys

values = {}
for line in pathlib.Path(sys.argv[1]).read_text().splitlines()[1:]:
    name, value = line.split(",", 1)
    values[name] = float(value)
assert values["anisotropy_ratio"] >= 1000.0
assert values["maximum_nodal_error"] < 0.15
assert values["linear_residual_inf"] < 5.0e-10
assert values["discrete_energy"] > 0.0
PY

echo "field-aligned anisotropic 3-D gallery emits solution before diagnostics"
