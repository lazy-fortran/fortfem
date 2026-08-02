#!/usr/bin/env bash
set -euo pipefail

repository_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
output_directory="$repository_dir/output/example/tetra_nedelec_p_convergence"

rm -rf -- "$output_directory"
(
    cd "$repository_dir"
    timeout --foreground 60s fpm run --example tetra_nedelec_p_convergence \
        >/tmp/fortfem-tetra-nedelec-gallery.log
)

test -s "$output_directory/nedelec_field_3d.png"
test -s "$output_directory/nedelec_field_slice_2d.png"
test -s "$output_directory/p_convergence_1d.png"
test -s "$output_directory/gallery_sequence.txt"

mapfile -t stages <"$output_directory/gallery_sequence.txt"
test "${#stages[@]}" -eq 2
test "${stages[0]}" = physical_solution
test "${stages[1]}" = diagnostics

primary_name=$(awk '$1 == "tetra_nedelec_p_convergence" {print $2}' \
    "$repository_dir/doc/examples/primary_plots.txt")
test "$primary_name" = nedelec_field_slice_2d.png

# The primary artifact must be a rendered image and must not be the diagnostic
# convergence curve that is emitted last by the example.
test "$(file -b "$output_directory/$primary_name" | cut -d' ' -f1)" = PNG
test "$primary_name" != p_convergence_1d.png

echo "tetrahedral Nedelec gallery emits vector solution before diagnostics"
