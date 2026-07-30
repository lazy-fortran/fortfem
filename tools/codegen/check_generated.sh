#!/usr/bin/env bash
set -euo pipefail

codegen_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repository_dir=$(cd "$codegen_dir/../.." && pwd)
temporary_dir=$(mktemp -d)
trap 'rm -rf -- "$temporary_dir"' EXIT

(
    cd "$codegen_dir"
    FORTFEM_CODEGEN_OUTPUT_DIR="$temporary_dir" ./generate.sh
)

for order in 1 2 3 4; do
    filename="fortfem_tetra_nedelec_candidates_order_${order}.f90"
    cmp -- "$repository_dir/src/generated/$filename" \
        "$temporary_dir/$filename"
done

for degree in 0 1 2 3 4; do
    filename="fortfem_tetra_rt_candidates_degree_${degree}.f90"
    cmp -- "$repository_dir/src/generated/$filename" \
        "$temporary_dir/$filename"
done

filename="fortfem_tetra_rt_coefficients.f90"
cmp -- "$repository_dir/src/generated/$filename" \
    "$temporary_dir/$filename"

filename="fortfem_tetra_face_moment_transforms.f90"
cmp -- "$repository_dir/src/generated/$filename" \
    "$temporary_dir/$filename"

filename="fortfem_tetra_nedelec_coefficients.f90"
cmp -- "$repository_dir/src/generated/$filename" \
    "$temporary_dir/$filename"

filename="fortfem_tetra_modal_vector_identities.f90"
cmp -- "$repository_dir/src/generated/$filename" \
    "$temporary_dir/$filename"

filename="fortfem_toroidal_coordinates.f90"
cmp -- "$repository_dir/src/generated/$filename" \
    "$temporary_dir/$filename"

filename="fortfem_interoperability_oracles.f90"
cmp -- "$repository_dir/src/generated/$filename" \
    "$temporary_dir/$filename"

filename="fortfem_tetra_h1_oracle.f90"
cmp -- "$repository_dir/src/generated/$filename" \
    "$temporary_dir/$filename"

filename="fortfem_magnetic_curvilinear_coefficients_2d.f90"
cmp -- "$repository_dir/src/generated/$filename" \
    "$temporary_dir/$filename"

filename="fortfem_nurbs_geometry_products.f90"
cmp -- "$repository_dir/src/generated/$filename" \
    "$temporary_dir/$filename"

echo "generated FortFEM kernels match committed sources"
