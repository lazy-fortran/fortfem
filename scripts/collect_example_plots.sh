#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
    echo "usage: $0 SOURCE_DIR DESTINATION_DIR" >&2
    exit 2
fi

source_dir=$1
destination_dir=$2
if [[ ! -d "$source_dir" ]]; then
    echo "example plot source is not a directory: $source_dir" >&2
    exit 1
fi

mkdir -p "$destination_dir"
cp -a -- "$source_dir"/. "$destination_dir"/

plot_count=$(find "$destination_dir" -type f -iname '*.png' | wc -l)
if [[ "$plot_count" -eq 0 ]]; then
    echo "no example PNG plots found under: $source_dir" >&2
    exit 1
fi

empty_plot=$(find "$destination_dir" -type f \
    \( -iname '*.png' -o -iname '*.svg' \) -size 0 -print -quit)
if [[ -n "$empty_plot" ]]; then
    echo "example plot is empty: $empty_plot" >&2
    exit 1
fi

echo "collected $plot_count example PNG plots"
