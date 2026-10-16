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

# The committed kernels are reproducible only from the locked toolchain.
while read -r name revision; do
    case "$name" in ''|'#'*) continue ;; esac
    if [ "$name" = fortad ]; then
        checkout=$fortad_repo
    else
        checkout="$fortad_repo/../$name"
    fi
    actual=$(git -C "$checkout" rev-parse HEAD)
    if [ "$actual" != "$revision" ]; then
        echo "$name checkout $checkout is at $actual, fortad.lock wants $revision" >&2
        exit 1
    fi
    if [ -n "$(git -C "$checkout" status --porcelain --untracked-files=no)" ]; then
        echo "$name checkout $checkout has uncommitted changes" >&2
        exit 1
    fi
done < "$here/fortad.lock"

# fo is the build driver, and using it here is not only convention: fpm leaves
# an app binary unrelinked often enough that a stale fortad silently regenerates
# the previous kernels, which then differ from the source they claim to come
# from and pass the equivalence test anyway.
# FORTAD_BIN may name a binary already built from the locked checkout.
if [ -n "${FORTAD_BIN:-}" ]; then
    fortad_bin=$FORTAD_BIN
else
    if [ ! -x "$fortad_repo/build/fo/bin/fortad" ]; then
        ( cd "$fortad_repo" && fo build >/dev/null )
    fi
    fortad_bin="$fortad_repo/build/fo/bin/fortad"
fi

python3 "$here/extract.py"
FORTAD_BIN="$fortad_bin" python3 "$here/emit.py"

echo "regenerated fortad kernels in src/generated/fortad"
