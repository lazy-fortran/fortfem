#!/usr/bin/env bash
set -euo pipefail

repository_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
fixture_dir=$(mktemp -d "${TMPDIR:-/tmp}/fortfem-primary-plots.XXXXXX")
trap 'rm -rf "$fixture_dir"' EXIT
mkdir -p "$fixture_dir"/{scripts,doc/examples,bin}
cp "$repository_dir/scripts/generate_primary_plots.sh" \
    "$repository_dir/scripts/stage_gallery_media.sh" "$fixture_dir/scripts/"
cat > "$fixture_dir/doc/examples/primary_plots.txt" <<'EOF'
# Synthetic examples exercise both output conventions.
root_case field.png
nested_case geometry.png
interoperability_benchmarks manufactured_poisson_solution_2d.png
EOF
cat > "$fixture_dir/bin/fpm" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
test -z "${DISPLAY:-}"
test "$MPLBACKEND" = Agg
test "$1" = run
test "$2" = --example
test "$6" = --runner
name=$3
case "$name" in
    root_case) plot=field.png; directory=. ;;
    nested_case) plot=geometry.png; directory=output/example/$name ;;
    interoperability_benchmarks)
        test "$8" = --
        test "$9" = --primary-only
        directory=${10}
        plot=manufactured_poisson_solution_2d.png
        ;;
    *) exit 1 ;;
esac
mkdir -p "$directory"
if [[ "${HANG_TARGET:-}" == "$name" ]]; then
    read -r -a runner <<< "$7"
    "${runner[@]}" sleep 10
fi
if [[ "${FAIL_TARGET:-}" == "$name" ]]; then
    exit 23
fi
if [[ "${MISSING_TARGET:-}" != "$name" ]]; then
    printf 'synthetic image for %s\n' "$name" > "$directory/$plot"
fi
printf 'diagnostic image\n' > "$directory/diagnostic.png"
EOF
chmod +x "$fixture_dir/bin/fpm"
export PATH="$fixture_dir/bin:$PATH"
generator="$fixture_dir/scripts/generate_primary_plots.sh"

bash "$generator"
destination="$fixture_dir/artifacts/plots/examples"
test "$(find "$destination" -type f | wc -l)" -eq 3
grep -Fq 'root_case' "$destination/root_case/field.png"
grep -Fq 'nested_case' "$destination/nested_case/geometry.png"
test -s "$destination/interoperability_benchmarks/manufactured_poisson_solution_2d.png"

# Reruns discard obsolete published media and cannot accept a stale source.
printf 'obsolete\n' > "$destination/obsolete.png"
bash "$generator"
test ! -e "$destination/obsolete.png"
if MISSING_TARGET=nested_case bash "$generator" > "$fixture_dir/missing.log" 2>&1; then
    echo "missing primary plot was accepted" >&2
    exit 1
fi
grep -Fq 'primary_plots.txt:3' "$fixture_dir/missing.log"
grep -Fq 'missing physical primary plot' "$fixture_dir/missing.log"
if FAIL_TARGET=root_case bash "$generator" > "$fixture_dir/failure.log" 2>&1; then
    echo "failed example was accepted" >&2
    exit 1
fi
grep -Fq 'primary_plots.txt:2' "$fixture_dir/failure.log"
grep -Fq 'example failed or timed out' "$fixture_dir/failure.log"
if EXAMPLE_TIMEOUT=0.1 HANG_TARGET=root_case bash "$generator" > "$fixture_dir/timeout.log" 2>&1; then
    echo "hung example was accepted" >&2
    exit 1
fi
grep -Fq 'example failed or timed out' "$fixture_dir/timeout.log"

echo "primary plot generation stages only fresh mapped plots and rejects failures"
