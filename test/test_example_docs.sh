#!/usr/bin/env bash
set -euo pipefail

repository_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
index_file="$repository_dir/doc/examples/index.md"
generated_dir="$repository_dir/doc/examples/generated"
readme_file="$repository_dir/README.md"
order_file="$repository_dir/doc/examples/gallery_order.txt"
primary_file="$repository_dir/doc/examples/primary_plots.txt"
artifacts_dir="$repository_dir/artifacts/plots/examples"

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
mapfile -t gallery_cards < <(
    sed -n 's/^<article class="example-card" data-example="\([^"]*\)">$/\1/p' \
        "$index_file"
)
mapfile -t ordered < <(
    sed -e '/^[[:space:]]*#/d' -e '/^[[:space:]]*$/d' "$order_file"
)
mapfile -t primary_entries < <(
    sed -e '/^[[:space:]]*#/d' -e '/^[[:space:]]*$/d' "$primary_file"
)

if ! diff -u \
        <(printf '%s\n' "${expected[@]}") \
        <(printf '%s\n' "${documented[@]}"); then
    echo "example documentation index does not cover every Fortran example" >&2
    exit 1
fi

if ! diff -u \
        <(printf '%s\n' "${ordered[@]}") \
        <(printf '%s\n' "${gallery_cards[@]}"); then
    echo "example gallery is not in the documented learning order" >&2
    exit 1
fi

if ! diff -u \
        <(printf '%s\n' "${expected[@]}") \
        <(printf '%s\n' "${ordered[@]}" | sort -u); then
    echo "example gallery order does not cover every Fortran example" >&2
    exit 1
fi

mapfile -t primary_names < <(printf '%s\n' "${primary_entries[@]}" | awk '{print $1}' | sort -u)
if ! diff -u \
        <(printf '%s\n' "${expected[@]}") \
        <(printf '%s\n' "${primary_names[@]}"); then
    echo "primary plot map does not cover every Fortran example" >&2
    exit 1
fi

test "${gallery_cards[0]}" = simple_poisson
test "$(grep -c '<img class="example-card-image"' "$index_file")" \
    -eq "${#expected[@]}"
grep -Fq "classic Poisson equation" "$index_file"

for example_name in "${expected[@]}"; do
    page="$generated_dir/$example_name.md"
    test -s "$page"
    grep -Fq "fpm run --example $example_name" "$page"
    grep -Fq '```fortran' "$page"
    grep -Fq "program $example_name" "$page"
    if grep -Fq '*No plot artifact is produced by this example.*' "$page"; then
        echo "example page has no generated plot: $example_name" >&2
        exit 1
    fi
    if ! grep -Eq '\.(png|svg)\)' "$page"; then
        echo "example page has no generated plot or cover link: $example_name" >&2
        exit 1
    fi
    primary_name=$(awk -v name="$example_name" '$1 == name {print $2}' "$primary_file")
    test -n "$primary_name"
    case "$primary_name" in
        *convergence*|*error*|*timing*|*residual*|*accuracy*|*dofs*)
            echo "diagnostic plot cannot be the first gallery view: $example_name/$primary_name" >&2
            exit 1
            ;;
    esac
    # A generated cover is not a substitute for the physical first view.  A
    # mapped source artifact must be present and non-empty for every example;
    # otherwise the gallery can silently fall back to an alphabetical
    # convergence or timing plot.
    test -s "$artifacts_dir/$example_name/$primary_name"
    test -s "$artifacts_dir/$example_name/primary.png"
    first_plot=$(sed -n 's/^!\[[^]]*\](\([^)]*\)).*$/\1/p' "$page" | sed -n '1p')
    test "$first_plot" = "../../media/examples/$example_name/primary.png"
done

if grep -n '[[:blank:]]$' "$index_file" "$generated_dir/index.md"; then
    echo "example documentation indexes contain trailing whitespace" >&2
    exit 1
fi

if grep -Fq 'itpplasma.github.io/fortfem' "$readme_file"; then
    echo "README links to the retired GitHub Pages organization" >&2
    exit 1
fi
grep -Fq \
    'https://lazy-fortran.github.io/fortfem/page/examples/index.html' \
    "$readme_file"

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
