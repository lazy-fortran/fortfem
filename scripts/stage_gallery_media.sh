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
if ! cp -a -- "$source_dir"/. "$destination_dir"/; then
    echo "failed to copy gallery media into: $destination_dir" >&2
    exit 1
fi

mapfile -d '' source_files < <(find "$source_dir" -type f -print0)
for source in "${source_files[@]}"; do
    relative=${source#"$source_dir"/}
    target="$destination_dir/$relative"
    if [[ ! -f "$target" ]]; then
        echo "gallery media was not staged: $relative" >&2
        exit 1
    fi
    if ! cmp -s -- "$source" "$target"; then
        echo "gallery media changed during staging: $relative" >&2
        exit 1
    fi
done

empty_media=$(find "$destination_dir" -type f \
    \( -name '*.png' -o -name '*.svg' \) -size 0 -print -quit)
if [[ -n "$empty_media" ]]; then
    echo "gallery media is empty: $empty_media" >&2
    exit 1
fi
