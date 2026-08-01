#!/usr/bin/env bash
set -euo pipefail

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repository_dir=$(cd "$script_dir/.." && pwd)
example_dir="$repository_dir/example"
doc_examples_dir="$repository_dir/doc/examples"
generated_dir="$doc_examples_dir/generated"
artifacts_dir="$repository_dir/artifacts/plots/examples"
order_file="$doc_examples_dir/gallery_order.txt"
primary_file="$doc_examples_dir/primary_plots.txt"

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

python3 "$script_dir/generate_example_covers.py"

# Keep the gallery cover independent of filename sorting.  The producing
# example writes the physical field/geometry plot named by primary_plots.txt;
# this stable alias is ignored with the other generated artifacts.
for example_name in "${example_names[@]}"; do
    primary_name=${primary_for_name[$example_name]}
    primary_source="$artifacts_dir/$example_name/$primary_name"
    if [[ -f "$primary_source" ]]; then
        cp -f "$primary_source" "$artifacts_dir/$example_name/primary.png"
    fi
done

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
        'Plots and fallback covers are generated by GitHub Actions and are not' \
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

        plot_count=0
        if [[ -d "$artifacts_dir/$example_name" ]]; then
            primary_plot="$artifacts_dir/$example_name/primary.png"
            if [[ -f "$primary_plot" ]]; then
                plots=("$primary_plot")
                while IFS= read -r plot; do
                    [[ "$plot" == "$primary_plot" ]] && continue
                    plots+=("$plot")
                done < <(
                    find "$artifacts_dir/$example_name" -maxdepth 1 \
                        -type f -name '*.png' -print |
                        sort
                )
            else
                mapfile -t plots < <(
                    find "$artifacts_dir/$example_name" -maxdepth 1 \
                        -type f -name '*.png' -print |
                        sort
                )
            fi
            for plot in "${plots[@]}"; do
                plot_name=$(basename "$plot")
                printf '%s\n' \
                    "### $plot_name" \
                    '' \
                    "![$plot_name](../../media/examples/$example_name/$plot_name)" \
                    ''
                plot_count=$((plot_count + 1))
            done
            # Optional external benchmark runners may have no local numerical
            # artifact. Keep the provenance-labelled SVG cover visible instead
            # of publishing an empty page; numerical examples are required to
            # emit PNGs by the example-doc test.
            if [[ "$plot_count" -eq 0 && \
                -f "$artifacts_dir/$example_name/cover.svg" ]]; then
                printf '%s\n' \
                    '### cover.svg' \
                    '' \
                    "![cover.svg](../../media/examples/$example_name/cover.svg)" \
                    ''
                plot_count=1
            fi
        fi
        if [[ "$plot_count" -eq 0 ]]; then
            printf '%s\n\n' '*No plot artifact is produced by this example.*'
        fi
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
    preview="$artifacts_dir/$example_name/primary.png"
    if [[ ! -f "$preview" ]]; then
        preview=$(find "$artifacts_dir/$example_name" -maxdepth 1 \
            -type f -name '*.png' -print | sort | sed -n '1p')
    fi
    if [[ -z "$preview" ]]; then
        preview=$(find "$artifacts_dir/$example_name" -maxdepth 1 \
            -type f -name '*.svg' -print | sort | sed -n '1p')
    fi
    if [[ -n "$preview" ]]; then
        preview_name=$(basename "$preview")
        printf '%s\n' \
            "<img class=\"example-card-image\"" \
            " src=\"../media/examples/$example_name/$preview_name\"" \
            " alt=\"Plot preview for $example_name\" loading=\"lazy\">" \
            >> "$doc_examples_dir/index.md"
    else
        echo "gallery cover generation failed for $example_name" >&2
        exit 1
    fi
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

echo "generated documentation for ${#example_names[@]} executable examples"
