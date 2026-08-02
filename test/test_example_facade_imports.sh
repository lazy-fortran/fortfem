#!/usr/bin/env bash
set -euo pipefail

repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "$repo_root"
python3 tools/check_example_facade_imports.py --root "$repo_root"

# Keep the independent behavioral clients on the same canonical surface as
# their gallery programs.  These are intentionally explicit: a new client
# must be assigned a facade in the migration map rather than silently using
# the compatibility umbrella.
for source in \
    test/test_advanced_facade_clients.f90 \
    test/test_bspline_polar_axis.f90 \
    test/test_bem_sphere_surface_solution.f90 \
    test/test_maxwell_torus_fem_bem_solution.f90 \
    test/test_sphere_surface_mesh.f90 \
    test/test_solid_torus_tetra_mesh.f90 \
    test/test_torus_surface_mesh.f90; do
    if rg -n '^\s*use\s+fortfem_api\b' "$source"; then
        echo "ERROR: migrated advanced client still imports fortfem_api: $source" >&2
        exit 1
    fi
done
echo "advanced facade clients: PASS (7 no-umbrella tests)"
