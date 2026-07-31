#!/usr/bin/env bash
set -euo pipefail

codegen_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repository_dir=$(cd "$codegen_dir/../.." && pwd)
generated_dir=${FORTFEM_CODEGEN_OUTPUT_DIR:-"$repository_dir/src/generated"}
export FORTFEM_CODEGEN_OUTPUT_DIR="$generated_dir"

# The generator builds FortSym and many independent Fortran applications.
# Host-wide parallelism can exhaust memory on CI runners and make the
# provenance check fail with a compiler segfault.  Keep an explicit caller
# setting, but use a bounded default for reproducible generation.
export FO_JOBS="${FO_JOBS:-2}"

# Fo is the normal code-generation runner.  Some hosted runners have
# exhibited a front-end SIGSEGV while starting `fo build` for this dependency
# graph; FPM provides an equivalent deterministic fallback for CI without
# changing the generated kernels.
codegen_runner=${FORTFEM_CODEGEN_RUNNER:-fo}
if [[ "$codegen_runner" != "fo" && "$codegen_runner" != "fpm" ]]; then
    echo "unsupported FORTFEM_CODEGEN_RUNNER=$codegen_runner" >&2
    exit 2
fi

run_codegen() {
    local target=$1
    if [[ "$codegen_runner" == "fpm" ]]; then
        fpm run --target "$target" --profile release
    else
        fo exec --no-build "$target"
    fi
}

cd "$codegen_dir"
./check_fortsym_revision.sh
if [[ "$codegen_runner" == "fpm" ]]; then
    fpm build --profile release
else
    fo build
fi
run_codegen gen_tetra_nedelec_candidates
run_codegen gen_tetra_modal_vector_identities
run_codegen gen_tetra_rt_candidates
run_codegen gen_tetra_rt_coefficients
run_codegen gen_tetra_face_moment_transforms
run_codegen gen_tetra_nedelec_coefficients
run_codegen gen_toroidal_coordinates
run_codegen gen_torus_curved_panel_products
run_codegen gen_toroidal_poisson_products
run_codegen gen_sphere_curved_panel_products
run_codegen gen_interoperability_oracles
run_codegen gen_tetra_h1_oracle
run_codegen gen_magnetic_curvilinear_coefficients
run_codegen gen_nurbs_geometry_products
run_codegen gen_bspline_h1_geometry_products
run_codegen gen_cartesian_pml_products
run_codegen gen_planar_helmholtz_dtn_products
run_codegen gen_surface_triangle_geometry_products
run_codegen gen_laplace_bem_kernel_products
run_codegen gen_helmholtz_bem_kernel_products
run_codegen gen_laplace_singular_edge_products
run_codegen gen_helmholtz_bem_smooth_products
run_codegen gen_fci_parallel_products
run_codegen gen_fci_quartic_lagrange
run_codegen gen_cgl_pressure_tensor_products
run_codegen gen_cgl_pressure_divergence_products
run_codegen gen_fci_support_volume_products
run_codegen gen_field_aligned_flux_products

cd "$repository_dir"
fo fmt "$generated_dir/fortfem_tetra_face_moment_transforms.f90"
fo fmt "$generated_dir/fortfem_tetra_modal_vector_identities.f90"
fo fmt "$generated_dir/fortfem_tetra_nedelec_coefficients.f90"
fo fmt "$generated_dir"/fortfem_tetra_rt_candidates_degree_*.f90
fo fmt "$generated_dir/fortfem_tetra_rt_coefficients.f90"
fo fmt "$generated_dir/fortfem_toroidal_coordinates.f90"
fo fmt "$generated_dir/fortfem_torus_curved_panel_products.f90"
fo fmt "$generated_dir/fortfem_toroidal_poisson_products.f90"
fo fmt "$generated_dir/fortfem_sphere_curved_panel_products.f90"
fo fmt "$generated_dir/fortfem_interoperability_oracles.f90"
fo fmt "$generated_dir/fortfem_tetra_h1_oracle.f90"
fo fmt "$generated_dir/fortfem_magnetic_curvilinear_coefficients_2d.f90"
fo fmt "$generated_dir/fortfem_nurbs_geometry_products.f90"
fo fmt "$generated_dir/fortfem_bspline_h1_geometry_products.f90"
fo fmt "$generated_dir/fortfem_cartesian_pml_products.f90"
fo fmt "$generated_dir/fortfem_planar_helmholtz_dtn_products.f90"
fo fmt "$generated_dir/fortfem_surface_triangle_geometry_products.f90"
fo fmt "$generated_dir/fortfem_laplace_bem_kernel_products.f90"
fo fmt "$generated_dir/fortfem_helmholtz_bem_kernel_products.f90"
fo fmt "$generated_dir/fortfem_laplace_singular_edge_products.f90"
fo fmt "$generated_dir/fortfem_helmholtz_bem_smooth_products.f90"
fo fmt "$generated_dir/fortfem_fci_parallel_products.f90"
fo fmt "$generated_dir/fortfem_fci_quartic_lagrange.f90"
fo fmt "$generated_dir/fortfem_cgl_pressure_tensor_products.f90"
fo fmt "$generated_dir/fortfem_cgl_pressure_divergence_products.f90"
fo fmt "$generated_dir/fortfem_fci_support_volume_products.f90"
fo fmt "$generated_dir/fortfem_field_aligned_flux_products.f90"
