#!/usr/bin/env bash
set -euo pipefail

codegen_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repository_dir=$(cd "$codegen_dir/../.." && pwd)
generated_dir=${FORTFEM_CODEGEN_OUTPUT_DIR:-"$repository_dir/src/generated"}

cd "$codegen_dir"
./check_fortsym_revision.sh
fo build
fo exec --no-build gen_tetra_nedelec_candidates
fo exec --no-build gen_tetra_face_moment_transforms
fo exec --no-build gen_tetra_nedelec_coefficients

cd "$repository_dir"
fo fmt "$generated_dir/fortfem_tetra_face_moment_transforms.f90"
fo fmt "$generated_dir/fortfem_tetra_nedelec_coefficients.f90"
