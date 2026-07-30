#!/usr/bin/env bash
set -euo pipefail

repository_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
index_file="$repository_dir/doc/examples/index.md"
generated_dir="$repository_dir/doc/examples/generated"

mapfile -t expected < <(
    find "$repository_dir/example" -mindepth 1 -maxdepth 2 \
        -type f -name '*.f90' -printf '%f\n' |
        sed 's/\.f90$//' |
        sort -u
)
mapfile -t documented < <(
    sed -n 's/^- \[\([^]]*\)\](generated\/[^)]*\.html).*$/\1/p' \
        "$index_file" |
        sort -u
)

if ! diff -u \
        <(printf '%s\n' "${expected[@]}") \
        <(printf '%s\n' "${documented[@]}"); then
    echo "example documentation index does not cover every Fortran example" >&2
    exit 1
fi

for example_name in "${expected[@]}"; do
    page="$generated_dir/$example_name.md"
    test -s "$page"
    grep -Fq "fpm run --example $example_name" "$page"
    grep -Fq '```fortran' "$page"
    grep -Fq "program $example_name" "$page"
done

if grep -n '[[:blank:]]$' "$index_file" "$generated_dir/index.md"; then
    echo "example documentation indexes contain trailing whitespace" >&2
    exit 1
fi

artifacts_dir="$repository_dir/artifacts/plots/examples"
if [[ -d "$artifacts_dir" ]]; then
    while IFS= read -r plot; do
        example_name=$(basename "$(dirname "$plot")")
        plot_name=$(basename "$plot")
        grep -Fq \
            "../../media/examples/$example_name/$plot_name" \
            "$generated_dir/$example_name.md"
    done < <(find "$artifacts_dir" -mindepth 2 -maxdepth 2 \
        -type f -name '*.png' | sort)
fi

echo "example documentation covers ${#expected[@]} executable examples"
