#!/usr/bin/env bash
# Regenerate fortad derivative kernels for fortfem.
#
# fortad is a build-time source generator, not a runtime dependency: the
# generated .f90 files are committed and compile with any conforming Fortran
# compiler, so fortfem gains no new link-time or plugin dependency.
#
# The kernel list is not written here. tools/fortad/extract.py reads
# src/generated, takes every primal that ships with both a _jvp and a _vjp, and
# records which arguments each product treats as active - read off the committed
# signature rather than guessed. Keeping the list derived means a new operator
# in src/generated is picked up rather than forgotten.
set -euo pipefail

here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
root=$(cd "$here/../.." && pwd)
cd "$root"
fortad_repo=${FORTAD_REPO:-"$root/../fortad"}

fortad_bin=$(find "$fortad_repo/build" -name fortad -type f -perm -u+x 2>/dev/null | head -1)
if [ -z "$fortad_bin" ]; then
    ( cd "$fortad_repo" && fpm build >/dev/null )
    fortad_bin=$(find "$fortad_repo/build" -name fortad -type f -perm -u+x | head -1)
fi

python3 "$here/extract.py"
FORTAD_BIN="$fortad_bin" python3 "$here/emit.py"

echo "regenerated fortad kernels in src/generated/fortad"
