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

mkdir -p "$temporary_dir/lab"
(
    cd "$codegen_dir"
    FORTFEM_CODEGEN_OUTPUT_DIR="$temporary_dir" \
        FORTFEM_LAB_CODEGEN_OUTPUT_DIR="$temporary_dir/lab" ./generate.sh
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

filename="fortfem_verified_rt0_products.f90"
cmp -- "$repository_dir/src/generated/$filename" \
    "$temporary_dir/$filename"

filename="fortfem_verified_bernstein_products.f90"
cmp -- "$repository_dir/src/generated/$filename" \
    "$temporary_dir/$filename"

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

filename="fortfem_tetra_modal_component_curls_jvp.f90"
cmp -- "$repository_dir/src/generated/$filename" \
    "$temporary_dir/$filename"

filename="fortfem_tetra_modal_radial_products.f90"
cmp -- "$repository_dir/src/generated/$filename" \
    "$temporary_dir/$filename"

filename="fortfem_tetra_modal_radial_products_jvp.f90"
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

filename="fortfem_meridian_hcurl_equilibration_products.f90"
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

filename="fortfem_fci_decic_bezier_edge_area_products.f90"
cmp -- "$repository_dir/src/generated/$filename" \
    "$temporary_dir/$filename"

filename="fortfem_fci_hendecic_bezier_edge_area_products.f90"
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

filename="fortfem_reference_scalar_intervals.f90"
cmp -- "$repository_dir/src/generated/$filename" "$temporary_dir/$filename"

filename="fortfem_reference_basis_products.f90"
cmp -- "$repository_dir/src/generated/$filename" "$temporary_dir/$filename"

filename="fortfem_triangle_level_geometry.f90"
cmp -- "$repository_dir/src/generated/$filename" "$temporary_dir/$filename"

filename="fortfem_lagrange_product_jets.f90"
cmp -- "$repository_dir/src/generated/$filename" \
    "$temporary_dir/$filename"

filename="fortfem_mapped_oscillation_products.f90"
cmp -- "$repository_dir/src/generated/$filename" \
    "$temporary_dir/$filename"

filename="fortfem_equilibrated_flux_products.f90"
cmp -- "$repository_dir/src/generated/$filename" \
    "$temporary_dir/$filename"

filename="fortfem_equilibrated_nonconforming_products.f90"
cmp -- "$repository_dir/src/generated/$filename" \
    "$temporary_dir/$filename"

filename="fortfem_equilibrated_flux_3d_products.f90"
cmp -- "$repository_dir/src/generated/$filename" \
    "$temporary_dir/$filename"

filename="fortfem_hcurl_bound_products.f90"
cmp -- "$repository_dir/src/generated/$filename" \
    "$temporary_dir/$filename"

echo "generated FortFEM kernels match committed sources"

filename="fortfem_tetra_whitney_products.f90"
cmp -- "$repository_dir/src/generated/$filename" "$temporary_dir/$filename"

filename="fortfem_triangle_piola_products.f90"
cmp -- "$repository_dir/src/generated/$filename" "$temporary_dir/$filename"

filename="fortfem_triangle_covariant_primal.inc"
cmp -- "$repository_dir/src/generated/$filename" "$temporary_dir/$filename"

filename="fortfem_triangle_contravariant_primal.inc"
cmp -- "$repository_dir/src/generated/$filename" "$temporary_dir/$filename"

filename="fortfem_triangle_covariant_jvp.inc"
cmp -- "$repository_dir/src/generated/$filename" "$temporary_dir/$filename"

filename="fortfem_triangle_contravariant_jvp.inc"
cmp -- "$repository_dir/src/generated/$filename" "$temporary_dir/$filename"

filename="fortfem_affine_map_products.f90"
cmp -- "$repository_dir/src/generated/$filename" "$temporary_dir/$filename"

filename="fortfem_triangle_geometry.inc"
cmp -- "$repository_dir/src/generated/$filename" "$temporary_dir/$filename"

filename="fortfem_triangle_geometry_jvp.inc"
cmp -- "$repository_dir/src/generated/$filename" "$temporary_dir/$filename"

filename="fortfem_triangle_geometry_vjp.inc"
cmp -- "$repository_dir/src/generated/$filename" "$temporary_dir/$filename"

filename="fortfem_tetra_geometry.inc"
cmp -- "$repository_dir/src/generated/$filename" "$temporary_dir/$filename"

filename="fortfem_tetra_geometry_jvp.inc"
cmp -- "$repository_dir/src/generated/$filename" "$temporary_dir/$filename"

filename="fortfem_tetra_geometry_vjp.inc"
cmp -- "$repository_dir/src/generated/$filename" "$temporary_dir/$filename"

filename="fortfem_tetra_covariant.inc"
cmp -- "$repository_dir/src/generated/$filename" "$temporary_dir/$filename"

filename="fortfem_tetra_covariant_jvp.inc"
cmp -- "$repository_dir/src/generated/$filename" "$temporary_dir/$filename"

filename="fortfem_tetra_covariant_vjp.inc"
cmp -- "$repository_dir/src/generated/$filename" "$temporary_dir/$filename"

filename="fortfem_tetra_contravariant.inc"
cmp -- "$repository_dir/src/generated/$filename" "$temporary_dir/$filename"

filename="fortfem_tetra_contravariant_jvp.inc"
cmp -- "$repository_dir/src/generated/$filename" "$temporary_dir/$filename"

filename="fortfem_tetra_contravariant_vjp.inc"
cmp -- "$repository_dir/src/generated/$filename" "$temporary_dir/$filename"


filename="fortfem_polynomial_candidate_jets.f90"
cmp -- "$repository_dir/src/generated/$filename" "$temporary_dir/$filename"

filename="fortfem_power_multiply.inc"
cmp -- "$repository_dir/src/generated/$filename" "$temporary_dir/$filename"

filename="fortfem_monomial_jet3_order2.inc"
cmp -- "$repository_dir/src/generated/$filename" "$temporary_dir/$filename"

filename="fortfem_monomial_jet3_order1.inc"
cmp -- "$repository_dir/src/generated/$filename" "$temporary_dir/$filename"

filename="fortfem_monomial_jet2_order2.inc"
cmp -- "$repository_dir/src/generated/$filename" "$temporary_dir/$filename"

filename="fortfem_monomial_jet2_order1.inc"
cmp -- "$repository_dir/src/generated/$filename" "$temporary_dir/$filename"

filename="fortfem_cached_power_jet.inc"
cmp -- "$repository_dir/src/generated/$filename" "$temporary_dir/$filename"

filename="fortfem_triangle_candidate3_order2.inc"
cmp -- "$repository_dir/src/generated/$filename" "$temporary_dir/$filename"

filename="fortfem_triangle_candidate3_order1.inc"
cmp -- "$repository_dir/src/generated/$filename" "$temporary_dir/$filename"

filename="fortfem_triangle_candidate2_order2.inc"
cmp -- "$repository_dir/src/generated/$filename" "$temporary_dir/$filename"

filename="fortfem_triangle_candidate2_order1.inc"
cmp -- "$repository_dir/src/generated/$filename" "$temporary_dir/$filename"

filename="fortfem_triangle_candidate1_order2.inc"
cmp -- "$repository_dir/src/generated/$filename" "$temporary_dir/$filename"

filename="fortfem_triangle_candidate1_order1.inc"
cmp -- "$repository_dir/src/generated/$filename" "$temporary_dir/$filename"


for component in 1 2 3; do
    for order in 1 2; do
        filename="fortfem_tetra_component${component}_order${order}.inc"
        cmp -- "$repository_dir/src/generated/$filename" "$temporary_dir/$filename"
    done
done

filename="fortfem_meridian_fourier_equilibration_products.f90"
cmp -- "$repository_dir/src/generated/$filename" "$temporary_dir/$filename"

for extension in f90 cuh; do
    filename="fortfem_element_precision_products.$extension"
    cmp -- "$repository_dir/src/generated/$filename" "$temporary_dir/$filename"
done

filename="fortfem_bem_high_order_products.f90"
cmp -- "$repository_dir/src/generated/$filename" "$temporary_dir/$filename"

filename="fortfem_parametric_bem_products.f90"
cmp -- "$repository_dir/src/generated/$filename" \
    "$temporary_dir/$filename"

filename="fortfem_tensor_fem_bem_products.f90"
cmp -- "$repository_dir/src/generated/$filename" \
    "$temporary_dir/$filename"

filename="fortfem_harmonic_field_products.f90"
cmp -- "$repository_dir/src/generated/$filename" \
    "$temporary_dir/$filename"

filename="fortfem_exterior_harmonic_products.f90"
cmp -- "$repository_dir/src/generated/$filename" "$temporary_dir/$filename"

filename="fortfem_tensor_metric_products.f90"
cmp -- "$repository_dir/src/generated/$filename" "$temporary_dir/$filename"

filename="fortfem_lab_gs_source_products.f90"
cmp -- "$repository_dir/benchmark/lab/neq1/generated/$filename" \
    "$temporary_dir/lab/$filename"

for product in scalar geometry transport guenter_kernel vector square_disk; do
    filename="fortfem_lab_${product}.f90"
    cmp -- "$repository_dir/src/generated/$filename" "$temporary_dir/$filename"
done
