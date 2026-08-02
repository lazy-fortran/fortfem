#!/usr/bin/env bash
set -euo pipefail

repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
checker="$repo_root/scripts/check_generated_visibility.py"

python3 "$checker" --root "$repo_root"

fixtures="$repo_root/test/fixtures/generated_visibility"
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

run_negative bad_facade "public facade imports generated module"
run_negative bad_client "public client imports generated module"
run_negative bad_wrapper "must declare PRIVATE"
run_negative bad_export "leaks generated names"
echo "generated visibility negative fixtures passed"
