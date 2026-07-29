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

echo "generated FortFEM kernels match committed sources"
