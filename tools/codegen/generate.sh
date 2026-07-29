#!/usr/bin/env bash
set -euo pipefail

codegen_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

cd "$codegen_dir"
./check_fortsym_revision.sh
fo build
fo exec --no-build gen_tetra_nedelec_candidates
