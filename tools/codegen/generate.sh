#!/usr/bin/env bash
set -euo pipefail

codegen_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repository_dir=$(cd "$codegen_dir/../.." && pwd)
generated_dir=${FORTFEM_CODEGEN_OUTPUT_DIR:-"$repository_dir/src/generated"}

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
fo exec --no-build gen_interoperability_oracles
fo exec --no-build gen_tetra_h1_oracle
fo exec --no-build gen_magnetic_curvilinear_coefficients
fo exec --no-build gen_nurbs_geometry_products

cd "$repository_dir"
fo fmt "$generated_dir/fortfem_tetra_face_moment_transforms.f90"
fo fmt "$generated_dir/fortfem_tetra_modal_vector_identities.f90"
fo fmt "$generated_dir/fortfem_tetra_nedelec_coefficients.f90"
fo fmt "$generated_dir"/fortfem_tetra_rt_candidates_degree_*.f90
fo fmt "$generated_dir/fortfem_tetra_rt_coefficients.f90"
fo fmt "$generated_dir/fortfem_toroidal_coordinates.f90"
fo fmt "$generated_dir/fortfem_interoperability_oracles.f90"
fo fmt "$generated_dir/fortfem_tetra_h1_oracle.f90"
fo fmt "$generated_dir/fortfem_magnetic_curvilinear_coefficients_2d.f90"
fo fmt "$generated_dir/fortfem_nurbs_geometry_products.f90"
