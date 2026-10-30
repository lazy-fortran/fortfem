#!/usr/bin/env bash
set -euo pipefail

repository_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
generator="$repository_dir/scripts/generate_example_docs.sh"

# The committed pages must equal a fresh generation from tracked inputs alone;
# no plot artifacts are needed for this.
bash "$generator" --check

fixture_dir=$(mktemp -d "${TMPDIR:-/tmp}/fortfem-example-docs-current.XXXXXX")
trap 'rm -rf "$fixture_dir"' EXIT
mkdir -p "$fixture_dir/scripts" "$fixture_dir/doc"
cp -- "$generator" "$fixture_dir/scripts/"
cp -a -- "$repository_dir/example" "$fixture_dir/"
cp -a -- "$repository_dir/doc/examples" "$fixture_dir/doc/"

# An edited example source without regenerated pages is drift.
printf '! drift fixture\n' >>"$fixture_dir/example/simple_poisson.f90"
if bash "$fixture_dir/scripts/generate_example_docs.sh" --check \
        >/dev/null 2>&1; then
    echo "stale example page was accepted by the drift check" >&2
    exit 1
fi

# Without plot artifacts, a full generation must stop before rewriting pages.
cp -a -- "$fixture_dir/doc/examples" "$fixture_dir/pages-before"
if bash "$fixture_dir/scripts/generate_example_docs.sh" >/dev/null 2>&1; then
    echo "generation succeeded without the mapped primary plots" >&2
    exit 1
fi
diff -r "$fixture_dir/pages-before" "$fixture_dir/doc/examples"

echo "committed example pages are current and drift is rejected"
