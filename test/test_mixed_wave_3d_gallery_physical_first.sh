#!/usr/bin/env bash
set -euo pipefail

repository_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
output_directory="$repository_dir/output/example/mixed_wave_3d_structure"
log_file="${TMPDIR:-/tmp}/fortfem-mixed-wave-3d-gallery.log"

rm -rf -- "$output_directory"
(
    cd "$repository_dir"
    timeout --foreground "${EXAMPLE_TIMEOUT:-30}s" \
        fpm run --example mixed_wave_3d_structure >"$log_file"
)

test -s "$output_directory/mixed_wave_3d_trajectory_3d.png"
test -s "$output_directory/mixed_wave_3d_components_1d.png"
test -s "$output_directory/mixed_wave_3d_energy_1d.png"
test -s "$output_directory/benchmark.csv"
test -s "$output_directory/gallery_sequence.txt"
test "$(file -b \
    "$output_directory/mixed_wave_3d_trajectory_3d.png" | cut -d' ' -f1)" = PNG

mapfile -t stages <"$output_directory/gallery_sequence.txt"
test "${#stages[@]}" -eq 2
test "${stages[0]}" = physical_solution
test "${stages[1]}" = diagnostics

primary_name=$(awk '$1 == "mixed_wave_3d_structure" {print $2}' \
    "$repository_dir/doc/examples/primary_plots.txt")
test "$primary_name" = mixed_wave_3d_trajectory_3d.png

grep -Fq "mixed_wave_3d_trajectory_3d.png" "$output_directory/benchmark.csv"
echo "3D mixed-wave gallery emits the physical trajectory before diagnostics"
