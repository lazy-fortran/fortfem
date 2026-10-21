#!/usr/bin/env bash
set -euo pipefail

repository_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "$repository_dir"
primary_file=doc/examples/primary_plots.txt
scratch_dir=$(mktemp -d "${TMPDIR:-/tmp}/fortfem-primary-media.XXXXXX")
trap 'rm -rf "$scratch_dir"' EXIT
unset DISPLAY
export MPLBACKEND=Agg
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-2}
export OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS:-1}

line_number=0
while read -r example_name primary_name extra; do
    line_number=$((line_number + 1))
    [[ -z "$example_name" || "$example_name" == \#* ]] && continue
    if [[ -z "$primary_name" || -n "$extra" ]]; then
        echo "$primary_file:$line_number: invalid primary plot entry" >&2
        exit 1
    fi
    output_dir="output/example/$example_name"
    # Both historical root output and per-example directories are supported.
    # Remove the mapped source before execution so a stale PNG cannot pass.
    rm -f -- "$primary_name" "$output_dir/$primary_name"
    arguments=()
    if [[ "$example_name" == interoperability_benchmarks ]]; then
        arguments=(-- --primary-only "$output_dir")
    fi
    echo "Generating physical primary plot: $example_name/$primary_name"
    # FPM applies the runner after compilation, so a clean library build does
    # not consume the numerical example's execution budget.
    if ! fpm run --example "$example_name" --profile release \
            --runner "timeout --kill-after=10s ${EXAMPLE_TIMEOUT:-180}s" \
            "${arguments[@]}"; then
        echo "$primary_file:$line_number: example failed or timed out: $example_name" >&2
        exit 1
    fi
    source_plot="$output_dir/$primary_name"
    if [[ ! -s "$source_plot" ]]; then
        source_plot=$primary_name
    fi
    if [[ ! -s "$source_plot" ]]; then
        echo "$primary_file:$line_number: missing physical primary plot: $example_name/$primary_name" >&2
        exit 1
    fi
    mkdir -p "$scratch_dir/$example_name"
    cp -- "$source_plot" "$scratch_dir/$example_name/$primary_name"
done < "$primary_file"

# Publish only the mapped physical PNGs, discarding stale and supplemental media.
bash scripts/stage_gallery_media.sh "$scratch_dir" artifacts/plots/examples
