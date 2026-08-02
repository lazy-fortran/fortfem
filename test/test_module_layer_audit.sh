#!/usr/bin/env bash
set -euo pipefail

repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
checker="$repo_root/scripts/check_module_layers.py"

python3 "$checker" --root "$repo_root/src"

fixtures="$repo_root/test/fixtures/module_layer_audit"
run_negative() {
    local name=$1
    local expected=$2
    local output
    if output=$(python3 "$checker" --root "$fixtures/$name" 2>&1); then
        echo "negative fixture unexpectedly passed: $name" >&2
        return 1
    fi
    grep -Fq "$expected" <<<"$output" || {
        echo "negative fixture $name did not report $expected" >&2
        printf '%s\n' "$output" >&2
        return 1
    }
}

run_negative bad_umbrella "forbidden umbrella import"
run_negative bad_generated "forbidden generated-layer import"
run_negative bad_core_application "layer violation"
run_negative bad_cycle "module dependency cycle"
echo "module-layer negative fixtures passed"
