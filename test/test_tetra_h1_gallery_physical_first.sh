#!/usr/bin/env bash
set -euo pipefail

repository_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
output_directory="$repository_dir/output/example/tetra_h1_poisson"

rm -rf -- "$output_directory"
(
    cd "$repository_dir"
    timeout --foreground 30s fpm run --example tetra_h1_poisson \
        >/tmp/fortfem-tetra-h1-gallery.log
)

test -s "$output_directory/tetra_h1_poisson_solution_3d.png"
test -s "$output_directory/tetra_h1_poisson_convergence.png"
test -s "$output_directory/convergence.csv"
test -s "$output_directory/gallery_sequence.txt"

mapfile -t stages <"$output_directory/gallery_sequence.txt"
test "${#stages[@]}" -eq 2
test "${stages[0]}" = physical_solution
test "${stages[1]}" = diagnostics

primary_name=$(awk '$1 == "tetra_h1_poisson" {print $2}' \
    "$repository_dir/doc/examples/primary_plots.txt")
test "$primary_name" = tetra_h1_poisson_solution_3d.png

echo "tetrahedral H1 gallery emits physical solution before diagnostics"
