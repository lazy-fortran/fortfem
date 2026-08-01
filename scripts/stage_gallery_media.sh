#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
    echo "usage: $0 SOURCE_DIR DESTINATION_DIR" >&2
    exit 2
fi

source_dir=$1
destination_dir=$2
if [[ ! -d "$source_dir" ]]; then
    echo "gallery media source is not a directory: $source_dir" >&2
    exit 1
fi

mkdir -p "$destination_dir"
# Copy the complete tree in one operation.  This preserves nested example
# names and avoids a shell read loop silently losing the final artifact.
cp -a -- "$source_dir"/. "$destination_dir"/

if ! diff -qr -- "$source_dir" "$destination_dir"; then
    echo "gallery media tree differs after staging" >&2
    exit 1
fi

empty_media=$(find "$destination_dir" -type f \
    \( -name '*.png' -o -name '*.svg' \) -size 0 -print -quit)
if [[ -n "$empty_media" ]]; then
    echo "gallery media is empty: $empty_media" >&2
    exit 1
fi
