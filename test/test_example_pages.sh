#!/usr/bin/env bash
set -euo pipefail

repository_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
page_root="$repository_dir/build/doc/page"
examples_root="$page_root/examples"

test -s "$examples_root/index.html"
test -s "$examples_root/generated/index.html"

mapfile -t expected < <(
    find "$repository_dir/example" -mindepth 1 -maxdepth 2 \
        -type f -name '*.f90' -printf '%f\n' |
        sed 's/\.f90$//' |
        sort -u
)

for example_name in "${expected[@]}"; do
    page="$examples_root/generated/$example_name.html"
    test -s "$page"
    grep -Fq "generated/$example_name.html" "$examples_root/index.html"
    grep -Fq "fpm" "$page"
    grep -Fq "$example_name" "$page"
done

while IFS= read -r relative_plot; do
    relative_plot=${relative_plot#src=\"}
    relative_plot=${relative_plot%\"}
    resolved=$(realpath -m "$examples_root/generated/$relative_plot")
    case "$resolved" in
        "$page_root"/media/examples/*) ;;
        *)
            echo "example image escapes the gallery media root: $relative_plot" >&2
            exit 1
            ;;
    esac
    test -s "$resolved"
done < <(
    grep -rho 'src="\.\./\.\./media/examples/[^"]*\.png"' \
        "$examples_root/generated"
)

while IFS= read -r relative_plot; do
    relative_plot=${relative_plot#src=\"}
    relative_plot=${relative_plot%\"}
    resolved=$(realpath -m "$examples_root/$relative_plot")
    case "$resolved" in
        "$page_root"/media/examples/*) ;;
        *)
            echo "gallery preview escapes the gallery media root: $relative_plot" >&2
            exit 1
            ;;
    esac
    test -s "$resolved"
done < <(
    grep -o 'src="\.\./media/examples/[^"]*\.png"' \
        "$examples_root/index.html"
)

echo "built GitHub Pages gallery covers ${#expected[@]} executable examples"
