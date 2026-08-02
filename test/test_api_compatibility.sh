#!/usr/bin/env bash
set -euo pipefail

repository_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
checker="$repository_dir/tools/check_api_compatibility.py"
fixture="$repository_dir/test/fixtures/api_compatibility_duplicate"

python3 "$checker" --root "$repository_dir"

output=$(python3 "$checker" --root "$fixture" 2>&1 || true)
if ! grep -Fq "duplicate canonical export" <<<"$output"; then
    echo "duplicate facade fixture did not report duplicate canonical export" >&2
    printf '%s\n' "$output" >&2
    exit 1
fi

echo "API compatibility gate passed"
