#!/usr/bin/env bash
set -euo pipefail

codegen_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repository_dir=$(cd "$codegen_dir/../.." && pwd)
generated_dir=${FORTFEM_CODEGEN_OUTPUT_DIR:-"$repository_dir/src/generated"}
export FORTFEM_CODEGEN_OUTPUT_DIR="$generated_dir"

cd "$codegen_dir"
./check_fortsym_revision.sh
fo build
fo exec --no-build gen_tetra_nedelec_candidates
fo exec --no-build gen_tetra_modal_vector_identities
fo exec --no-build gen_tetra_rt_candidates
fo exec --no-build gen_tetra_rt_coefficients
fo exec --no-build gen_tetra_face_moment_transforms
fo exec --no-build gen_tetra_nedelec_coefficients
fo exec --no-build gen_toroidal_coordinates
fo exec --no-build gen_torus_curved_panel_products
fo exec --no-build gen_toroidal_poisson_products
fo exec --no-build gen_sphere_curved_panel_products
fo exec --no-build gen_interoperability_oracles
fo exec --no-build gen_tetra_h1_oracle
fo exec --no-build gen_magnetic_curvilinear_coefficients
fo exec --no-build gen_nurbs_geometry_products
fo exec --no-build gen_bspline_h1_geometry_products
fo exec --no-build gen_cartesian_pml_products
fo exec --no-build gen_planar_helmholtz_dtn_products
fo exec --no-build gen_surface_triangle_geometry_products
fo exec --no-build gen_laplace_bem_kernel_products
fo exec --no-build gen_helmholtz_bem_kernel_products
fo exec --no-build gen_laplace_singular_edge_products
fo exec --no-build gen_helmholtz_bem_smooth_products
fo exec --no-build gen_fci_parallel_products
fo exec --no-build gen_cgl_pressure_tensor_products

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
fo fmt "$generated_dir/fortfem_cgl_pressure_tensor_products.f90"
