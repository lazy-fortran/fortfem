#!/usr/bin/env bash
set -euo pipefail

codegen_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repository_dir=$(cd "$codegen_dir/../.." && pwd)
temporary_root=${FORTFEM_CODEGEN_TMP_ROOT:-/mnt/storage}
if [[ -e "$temporary_root" ]]; then
    root_is_usable=false
    if [[ -d "$temporary_root" && -w "$temporary_root" ]]; then
        root_is_usable=true
    fi
else
    root_parent=$(dirname "$temporary_root")
    root_is_usable=false
    if [[ -d "$root_parent" && -w "$root_parent" ]]; then
        root_is_usable=true
    fi
fi
if [[ "$root_is_usable" != true ]]; then
    temporary_root=${TMPDIR:-/tmp}
fi
mkdir -p "$temporary_root"
temporary_dir=$(mktemp -d "$temporary_root/fortfem-codegen-check.XXXXXX")
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

filename="fortfem_tetra_modal_vector_identities_jvp.f90"
cmp -- "$repository_dir/src/generated/$filename" \
    "$temporary_dir/$filename"

filename="fortfem_toroidal_coordinates.f90"
cmp -- "$repository_dir/src/generated/$filename" \
    "$temporary_dir/$filename"

filename="fortfem_torus_curved_panel_products.f90"
cmp -- "$repository_dir/src/generated/$filename" \
    "$temporary_dir/$filename"

filename="fortfem_toroidal_poisson_products.f90"
cmp -- "$repository_dir/src/generated/$filename" \
    "$temporary_dir/$filename"

filename="fortfem_sphere_curved_panel_products.f90"
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

filename="fortfem_bspline_h1_geometry_products.f90"
cmp -- "$repository_dir/src/generated/$filename" \
    "$temporary_dir/$filename"

filename="fortfem_cartesian_pml_products.f90"
cmp -- "$repository_dir/src/generated/$filename" \
    "$temporary_dir/$filename"

filename="fortfem_planar_helmholtz_dtn_products.f90"
cmp -- "$repository_dir/src/generated/$filename" \
    "$temporary_dir/$filename"

filename="fortfem_fci_parallel_products.f90"
cmp -- "$repository_dir/src/generated/$filename" \
    "$temporary_dir/$filename"

filename="fortfem_fci_power_flux_products.f90"
cmp -- "$repository_dir/src/generated/$filename" \
    "$temporary_dir/$filename"

filename="fortfem_fci_quartic_lagrange.f90"
cmp -- "$repository_dir/src/generated/$filename" \
    "$temporary_dir/$filename"

filename="fortfem_fci_quintic_lagrange.f90"
cmp -- "$repository_dir/src/generated/$filename" \
    "$temporary_dir/$filename"

filename="fortfem_fci_quadrilateral_area_products.f90"
cmp -- "$repository_dir/src/generated/$filename" \
    "$temporary_dir/$filename"

filename="fortfem_fci_curved_quadrilateral_area_products.f90"
cmp -- "$repository_dir/src/generated/$filename" \
    "$temporary_dir/$filename"

filename="fortfem_fci_polygon_edge_area_products.f90"
cmp -- "$repository_dir/src/generated/$filename" \
    "$temporary_dir/$filename"

filename="fortfem_fci_quadratic_bezier_edge_area_products.f90"
cmp -- "$repository_dir/src/generated/$filename" \
    "$temporary_dir/$filename"

filename="fortfem_fci_cubic_bezier_edge_area_products.f90"
cmp -- "$repository_dir/src/generated/$filename" \
    "$temporary_dir/$filename"

filename="fortfem_fci_quartic_bezier_edge_area_products.f90"
cmp -- "$repository_dir/src/generated/$filename" \
    "$temporary_dir/$filename"

filename="fortfem_fci_quintic_bezier_edge_area_products.f90"
cmp -- "$repository_dir/src/generated/$filename" \
    "$temporary_dir/$filename"

filename="fortfem_fci_sextic_bezier_edge_area_products.f90"
cmp -- "$repository_dir/src/generated/$filename" \
    "$temporary_dir/$filename"

filename="fortfem_fci_septic_bezier_edge_area_products.f90"
cmp -- "$repository_dir/src/generated/$filename" \
    "$temporary_dir/$filename"

filename="fortfem_fci_octic_bezier_edge_area_products.f90"
cmp -- "$repository_dir/src/generated/$filename" \
    "$temporary_dir/$filename"

filename="fortfem_fci_nonic_bezier_edge_area_products.f90"
cmp -- "$repository_dir/src/generated/$filename" \
    "$temporary_dir/$filename"

filename="fortfem_surface_integral_products.f90"
cmp -- "$repository_dir/src/generated/$filename" \
    "$temporary_dir/$filename"

filename="fortfem_surface_shape_objective_products.f90"
cmp -- "$repository_dir/src/generated/$filename" \
    "$temporary_dir/$filename"

filename="fortfem_tensor_power_split_products.f90"
cmp -- "$repository_dir/src/generated/$filename" \
    "$temporary_dir/$filename"

filename="fortfem_cgl_pressure_tensor_products.f90"
cmp -- "$repository_dir/src/generated/$filename" \
    "$temporary_dir/$filename"

filename="fortfem_cgl_pressure_divergence_products.f90"
cmp -- "$repository_dir/src/generated/$filename" \
    "$temporary_dir/$filename"

filename="fortfem_fci_support_volume_products.f90"
cmp -- "$repository_dir/src/generated/$filename" \
    "$temporary_dir/$filename"

filename="fortfem_field_aligned_flux_products.f90"
cmp -- "$repository_dir/src/generated/$filename" \
    "$temporary_dir/$filename"

filename="fortfem_field_aligned_hall_products.f90"
cmp -- "$repository_dir/src/generated/$filename" \
    "$temporary_dir/$filename"

filename="fortfem_block_graph_products.f90"
cmp -- "$repository_dir/src/generated/$filename" \
    "$temporary_dir/$filename"

filename="fortfem_force_balance_products.f90"
cmp -- "$repository_dir/src/generated/$filename" \
    "$temporary_dir/$filename"

filename="fortfem_regularized_surface_current_products.f90"
cmp -- "$repository_dir/src/generated/$filename" \
    "$temporary_dir/$filename"

echo "generated FortFEM kernels match committed sources"
