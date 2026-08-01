#!/usr/bin/env bash
set -euo pipefail

repository_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
collector="$repository_dir/scripts/collect_example_plots.sh"
test -x "$collector"

fixture_dir=$(mktemp -d "${TMPDIR:-/tmp}/fortfem-plot-collection.XXXXXX")
trap 'rm -rf "$fixture_dir"' EXIT
source_dir="$fixture_dir/source"
destination_dir="$fixture_dir/destination"

mkdir -p "$source_dir/simple_poisson" "$source_dir/vector_case"
printf 'solution image\n' >"$source_dir/simple_poisson/primary.png"
printf 'vector image\n' >"$source_dir/vector_case/field_slice.png"
printf 'metadata\n' >"$source_dir/vector_case/field.csv"

"$collector" "$source_dir" "$destination_dir"
cmp "$source_dir/simple_poisson/primary.png" \
    "$destination_dir/simple_poisson/primary.png"
cmp "$source_dir/vector_case/field_slice.png" \
    "$destination_dir/vector_case/field_slice.png"
cmp "$source_dir/vector_case/field.csv" \
    "$destination_dir/vector_case/field.csv"
test "$(find "$destination_dir" -type f -iname '*.png' | wc -l)" -eq 2

: >"$source_dir/vector_case/empty.svg"
if "$collector" "$source_dir" "$destination_dir"; then
    echo "empty SVG was accepted by plot collection" >&2
    exit 1
fi

empty_source="$fixture_dir/empty-source"
mkdir -p "$empty_source"
if "$collector" "$empty_source" "$fixture_dir/empty-destination"; then
    echo "plot collection accepted a source with no PNGs" >&2
    exit 1
fi

echo "example plot collection preserves nested media and rejects invalid output"
