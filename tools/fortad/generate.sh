#!/usr/bin/env bash
# Regenerate fortad derivative kernels for fortfem.
#
# fortad is a build-time source generator, not a runtime dependency: the
# generated .f90 files are committed and compile with any conforming Fortran
# compiler, so fortfem gains no new link-time or plugin dependency.
#
# The products match what the fortsym generator emits, argument for argument
# and contract for contract: tangent-only forward, cotangent-only reverse, and
# only the four Jacobian entries active. A signature difference would make the
# comparison a port of the caller rather than of the derivative.
set -euo pipefail

here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
root=$(cd "$here/../.." && pwd)
fortad_repo=${FORTAD_REPO:-"$root/../fortad"}

fortad_bin=$(find "$fortad_repo/build" -name fortad -type f -perm -u+x 2>/dev/null | head -1)
if [ -z "$fortad_bin" ]; then
    ( cd "$fortad_repo" && fpm build >/dev/null )
    fortad_bin=$(find "$fortad_repo/build" -name fortad -type f -perm -u+x | head -1)
fi

out="$root/src/generated"
mkdir -p "$out"
active=jacobian_11,jacobian_21,jacobian_12,jacobian_22

"$fortad_bin" --indep "$active" --no-primal \
    --name fortfem_h1_geometry_jvp_fortad \
    --module fortfem_fortad_h1_geometry_jvp \
    -o "$out/fortfem_h1_geometry_jvp_fortad.f90" \
    "$here/kernels/h1_geometry.f90"

"$fortad_bin" --mode reverse --indep "$active" --no-primal \
    --name fortfem_h1_geometry_vjp_fortad \
    --module fortfem_fortad_h1_geometry_vjp \
    -o "$out/fortfem_h1_geometry_vjp_fortad.f90" \
    "$here/kernels/h1_geometry.f90"

echo "regenerated fortad kernels in src/generated"
