#!/usr/bin/env bash
set -euo pipefail

repository_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
release_checker="$repository_dir/scripts/check_api_release_gate.py"
layer_gate="$repository_dir/test/test_module_layer_audit.sh"
generated_gate="$repository_dir/test/test_generated_visibility.sh"
fixtures="$repository_dir/tools/fixtures/api_release_gate"

python3 "$release_checker" --root "$repository_dir"
bash "$layer_gate"
bash "$generated_gate"

run_negative() {
    local label=$1
    local expected=$2
    shift 2
    local output
    if output=$("$@" 2>&1); then
        echo "negative fixture unexpectedly passed: $label" >&2
        return 1
    fi
    grep -Fq "$expected" <<<"$output" || {
        echo "negative fixture $label did not report $expected" >&2
        printf '%s\n' "$output" >&2
        return 1
    }
}

run_negative stale-name \
    "stale internal API spelling" \
    python3 "$release_checker" \
    --root "$fixtures/stale_name" --stale-only
run_negative stale-inventory \
    "public API inventory is stale" \
    python3 "$release_checker" \
    --root "$repository_dir" --inventory "$fixtures/stale_inventory.md" \
    --inventory-only

# Keep the release path intentionally small: these two analytical tests cover
# value/JVP/VJP parity and the larger-domain value/JVP closure.  Full galleries,
# external-code comparisons and cross-compiler jobs run asynchronously.
timeout --foreground "${API_RELEASE_TEST_TIMEOUT:-10}s" \
    fo test test_boundary_operator_parity test_larger_domain_parity

echo "API release gate passed (fast inventory/layer/generated/derivative path)"
