#!/usr/bin/env bash
set -euo pipefail

repository_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
stager="$repository_dir/scripts/stage_gallery_media.sh"
test -x "$stager"

fixture_dir=$(mktemp -d "${TMPDIR:-/tmp}/fortfem-gallery-stage.XXXXXX")
trap 'rm -rf "$fixture_dir"' EXIT
source_dir="$fixture_dir/source"
destination_dir="$fixture_dir/destination"

mkdir -p "$source_dir/nested"
printf 'png fixture\n' >"$source_dir/solution.png"
printf '<svg/>\n' >"$source_dir/nested/mesh.svg"
printf 'x,y\n0,1\n' >"$source_dir/nested/values.csv"

"$stager" "$source_dir" "$destination_dir"
cmp "$source_dir/solution.png" "$destination_dir/solution.png"
cmp "$source_dir/nested/mesh.svg" "$destination_dir/nested/mesh.svg"
cmp "$source_dir/nested/values.csv" "$destination_dir/nested/values.csv"

# Staging is intentionally idempotent so a rerun cannot leave stale gallery
# media behind.
printf 'stale image\n' >"$destination_dir/stale.png"
mkdir -p "$destination_dir/obsolete-example"
printf 'stale nested image\n' >"$destination_dir/obsolete-example/old.png"
"$stager" "$source_dir" "$destination_dir"
test ! -e "$destination_dir/stale.png"
test ! -e "$destination_dir/obsolete-example/old.png"
test ! -e "$destination_dir/obsolete-example"

: >"$source_dir/empty.png"
if "$stager" "$source_dir" "$destination_dir"; then
    echo "empty image was accepted by gallery staging" >&2
    exit 1
fi

echo "gallery media staging preserves files and rejects empty images"
