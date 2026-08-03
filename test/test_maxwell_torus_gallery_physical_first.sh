#!/usr/bin/env bash
set -euo pipefail

repository_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
output_directory="$repository_dir/output/example/maxwell_torus_curved_scattering"

rm -rf -- "$output_directory"
(
    cd "$repository_dir"
    timeout --foreground "${EXAMPLE_TIMEOUT:-180}s" \
        fpm run --example maxwell_torus_curved_scattering \
        >/tmp/fortfem-maxwell-torus-gallery.log
)

test -s "$output_directory/maxwell_torus_solution_2d.png"
test -s "$output_directory/maxwell_torus_rcs_3d.png"
test -s "$output_directory/scattered_field.csv"
test -s "$repository_dir/example/maxwell_torus_curved_scattering/provenance.json"
test -s "$output_directory/gallery_sequence.txt"
test "$(file -b "$output_directory/maxwell_torus_rcs_3d.png" | cut -d' ' -f1)" = PNG

mapfile -t stages <"$output_directory/gallery_sequence.txt"
test "${#stages[@]}" -eq 2
test "${stages[0]}" = physical_solution
test "${stages[1]}" = diagnostics

# A connected surface is required for the 3-D radiation pattern; a point
# cloud silently drops the physical topology and is not a gallery solution.
grep -Fq "add_parametric_surface" \
    "$repository_dir/example/maxwell_torus_curved_scattering/maxwell_torus_curved_scattering.f90"
grep -Fq "normalized RCS surface" \
    "$repository_dir/example/maxwell_torus_curved_scattering/maxwell_torus_curved_scattering.f90"
grep -Fq "call quiver" \
    "$repository_dir/example/maxwell_torus_curved_scattering/maxwell_torus_curved_scattering.f90"

echo "toroidal Maxwell gallery emits a vector solution and connected radiation surface first"
