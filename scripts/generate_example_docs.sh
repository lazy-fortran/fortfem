#!/usr/bin/env bash
set -euo pipefail

# The pages depend only on tracked inputs (sources, READMEs, gallery order and
# primary plot map), never on which plot files a run happened to leave behind,
# so the committed pages and the published pages are identical.
#
# usage: generate_example_docs.sh [--check]
#   --check  regenerate into a scratch directory without plot artifacts and
#            fail if the committed pages under doc/examples differ.

check_only=false
if [[ $# -eq 1 && "$1" == --check ]]; then
    check_only=true
elif [[ $# -ne 0 ]]; then
    echo "usage: $0 [--check]" >&2
    exit 2
fi

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repository_dir=$(cd "$script_dir/.." && pwd)
example_dir="$repository_dir/example"
committed_examples_dir="$repository_dir/doc/examples"
artifacts_dir="$repository_dir/artifacts/plots/examples"
order_file="$committed_examples_dir/gallery_order.txt"
primary_file="$committed_examples_dir/primary_plots.txt"

if [[ "$check_only" == true ]]; then
    scratch_dir=$(mktemp -d "${TMPDIR:-/tmp}/fortfem-example-docs.XXXXXX")
    trap 'rm -rf "$scratch_dir"' EXIT
    doc_examples_dir="$scratch_dir/examples"
else
    doc_examples_dir=$committed_examples_dir
fi
generated_dir="$doc_examples_dir/generated"

mkdir -p "$generated_dir"

mapfile -t example_sources < <(
    find "$example_dir" -mindepth 1 -maxdepth 2 \
        -type f -name '*.f90' -print |
        sort
)
if [[ "${#example_sources[@]}" -eq 0 ]]; then
    echo "no Fortran examples found under $example_dir" >&2
    exit 1
fi

declare -A source_for_name
for source in "${example_sources[@]}"; do
    example_name=$(basename "$source" .f90)
    if [[ -n "${source_for_name[$example_name]:-}" ]]; then
        echo "duplicate example target name: $example_name" >&2
        exit 1
    fi
    source_for_name[$example_name]=$source
done

declare -A primary_for_name
while IFS=' ' read -r example_name primary_name extra; do
    [[ -z "$example_name" || "${example_name:0:1}" == "#" ]] && continue
    if [[ -n "$extra" || -z "$primary_name" ]]; then
        echo "invalid primary plot entry: $example_name $primary_name $extra" >&2
        exit 1
    fi
    if [[ -n "${primary_for_name[$example_name]:-}" ]]; then
        echo "primary plot repeats example: $example_name" >&2
        exit 1
    fi
    if [[ -z "${source_for_name[$example_name]:-}" ]]; then
        echo "primary plot names an unknown example: $example_name" >&2
        exit 1
    fi
    primary_for_name[$example_name]=$primary_name
done < "$primary_file"
if ! diff -q \
        <(printf '%s\n' "${!source_for_name[@]}" | sort) \
        <(printf '%s\n' "${!primary_for_name[@]}" | sort) >/dev/null; then
    echo "primary plot map must cover every executable example exactly once" >&2
    exit 1
fi

declare -A group_for_name
example_names=()
current_group=
while IFS= read -r line || [[ -n "$line" ]]; do
    line=${line#"${line%%[![:space:]]*}"}
    line=${line%"${line##*[![:space:]]}"}
    if [[ "$line" == \#* ]]; then
        current_group=${line#\#}
        current_group=${current_group#"${current_group%%[![:space:]]*}"}
    elif [[ -n "$line" ]]; then
        if [[ -z "${source_for_name[$line]:-}" ]]; then
            echo "gallery order names an unknown example: $line" >&2
            exit 1
        fi
        if [[ -n "${group_for_name[$line]:-}" ]]; then
            echo "gallery order repeats example: $line" >&2
            exit 1
        fi
        if [[ -z "$current_group" ]]; then
            echo "gallery example has no section: $line" >&2
            exit 1
        fi
        example_names+=("$line")
        group_for_name[$line]=$current_group
    fi
done < "$order_file"
if ! diff -q \
        <(printf '%s\n' "${!source_for_name[@]}" | sort) \
        <(printf '%s\n' "${example_names[@]}" | sort) >/dev/null; then
    echo "gallery order must cover every executable example exactly once" >&2
    exit 1
fi

if [[ "$check_only" == false ]]; then
    # Every page links its mapped primary plot, so a missing plot stops the
    # run before any page is written.  primary.png is the stable gallery-card
    # alias and is ignored with the other generated artifacts.
    for example_name in "${example_names[@]}"; do
        primary_name=${primary_for_name[$example_name]}
        primary_source="$artifacts_dir/$example_name/$primary_name"
        if [[ ! -s "$primary_source" ]]; then
            echo "missing physical primary plot for $example_name: $primary_source" >&2
            exit 1
        fi
        cp -f "$primary_source" "$artifacts_dir/$example_name/primary.png"
    done
fi

description_for()
{
    local source=$1
    local readme=$2
    local description

    description=
    if [[ -f "$readme" ]]; then
        description=$(awk '
            /^---$/ { in_frontmatter = !in_frontmatter; next }
            in_frontmatter || /^[[:space:]]*#/ || /^[[:space:]]*$/ { next }
            { print; exit }
        ' "$readme")
    fi
    if [[ -z "$description" ]]; then
        description=$(sed -n \
            '/^[[:space:]]*!/ { s/^[[:space:]]*!\+[[:space:]]*//; /^[[:space:]]*$/d; p; q; }' \
            "$source")
    fi
    if [[ -z "$description" ]]; then
        description="Executable FortFEM ${source##*/} example."
    fi
    printf '%s' "$description"
}

html_escape()
{
    printf '%s' "$1" |
        sed -e 's/&/\&amp;/g' -e 's/</\&lt;/g' -e 's/>/\&gt;/g'
}

{
    printf '%s\n' \
        '---' \
        'title: Example Source Code' \
        '---' \
        '' \
        '# Example Source Code and Plots' \
        '' \
        'Every executable example shipped with FortFEM is listed here.' \
        'Pages contain the complete source, usage command, and generated plots.' \
        ''
} > "$generated_dir/index.md"

{
    printf '%s\n' \
        '---' \
        'title: Examples' \
        '---' \
        '' \
        '# FortFEM Examples' \
        '' \
        'This gallery covers every Fortran program under `example/`.' \
        'Examples progress from the first Poisson solve to advanced coupled' \
        'boundary-element and toroidal problems.' \
        'Plots are generated by GitHub Actions and are not' \
        'stored in the source tree.' \
        '' \
        '## Running Examples' \
        '' \
        '```bash' \
        'fpm run --example <example_name>' \
        '```' \
        '' \
        '## Gallery' \
        ''
} > "$doc_examples_dir/index.md"

previous_group=
for example_name in "${example_names[@]}"; do
    source=${source_for_name[$example_name]}
    source_directory=$(dirname "$source")
    readme="$source_directory/README.md"
    if [[ "$source_directory" == "$example_dir" ]]; then
        readme="$example_dir/${example_name}_README.md"
    fi
    description=$(description_for "$source" "$readme")
    escaped_description=$(html_escape "$description")
    page="$generated_dir/$example_name.md"

    {
        printf '%s\n' \
            '---' \
            "title: $example_name Example" \
            '---' \
            '' \
            "# $example_name Example" \
            ''
        if [[ -f "$readme" ]]; then
            cat "$readme"
            printf '\n'
        else
            printf '%s\n\n' "$description"
        fi
        printf '%s\n' \
            '## Usage' \
            '' \
            '```bash' \
            "fpm run --example $example_name" \
            '```' \
            '' \
            '## Source Code' \
            '' \
            '```fortran'
        cat "$source"
        printf '%s\n' \
            '```' \
            '' \
            '## Generated Plots' \
            ''

        for plot_name in primary.png "${primary_for_name[$example_name]}"; do
            printf '%s\n' \
                "### $plot_name" \
                '' \
                "![$plot_name](../../media/examples/$example_name/$plot_name)" \
                ''
        done
        printf '%s\n' \
            '---' \
            '' \
            '[← Back to all examples](../index.html)'
    } > "$page"

    current_group=${group_for_name[$example_name]}
    if [[ "$current_group" != "$previous_group" ]]; then
        if [[ -n "$previous_group" ]]; then
            printf '%s\n' '</div>' '' >> "$doc_examples_dir/index.md"
        fi
        printf '%s\n' \
            "### $current_group" \
            '' \
            '<div class="example-gallery">' \
            >> "$doc_examples_dir/index.md"
        previous_group=$current_group
    fi
    printf '%s\n' \
        "<article class=\"example-card\" data-example=\"$example_name\">" \
        "<a class=\"example-card-preview\" href=\"generated/$example_name.html\">" \
        >> "$doc_examples_dir/index.md"
    printf '%s\n' \
        "<img class=\"example-card-image\"" \
        " src=\"../media/examples/$example_name/primary.png\"" \
        " alt=\"Plot preview for $example_name\" loading=\"lazy\">" \
        >> "$doc_examples_dir/index.md"
    printf '%s\n' \
        '</a>' \
        '<div class="example-card-body">' \
        "<h3><a href=\"generated/$example_name.html\">$example_name</a></h3>" \
        "<p>$escaped_description</p>" \
        '</div>' \
        '</article>' \
        >> "$doc_examples_dir/index.md"
    printf -- '- [%s](%s.html) - %s\n' \
        "$example_name" "$example_name" "$description" \
        >> "$generated_dir/index.md"
done

printf '%s\n' \
    '</div>' \
    '' \
    '## Complete Example Index' \
    '' \
    >> "$doc_examples_dir/index.md"
for example_name in "${example_names[@]}"; do
    source=${source_for_name[$example_name]}
    source_directory=$(dirname "$source")
    readme="$source_directory/README.md"
    if [[ "$source_directory" == "$example_dir" ]]; then
        readme="$example_dir/${example_name}_README.md"
    fi
    description=$(description_for "$source" "$readme")
    printf -- '- [%s](generated/%s.html) - %s\n' \
        "$example_name" "$example_name" "$description" \
        >> "$doc_examples_dir/index.md"
done

printf '\n[← Back to FortFEM documentation](../index.html)\n' \
    >> "$doc_examples_dir/index.md"
printf '\n[← Back to examples](../index.html)\n' \
    >> "$generated_dir/index.md"

if [[ "$check_only" == true ]]; then
    if ! diff -ru --exclude=gallery_order.txt --exclude=primary_plots.txt \
            "$committed_examples_dir" "$doc_examples_dir"; then
        echo "committed example pages are stale; run $0" >&2
        exit 1
    fi
    echo "committed documentation is current for ${#example_names[@]} executable examples"
    exit 0
fi

echo "generated documentation for ${#example_names[@]} executable examples"
