#!/usr/bin/env bash
set -euo pipefail

codegen_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
locked_revision=$(tr -d '[:space:]' < "$codegen_dir/fortsym.lock")
fortsym_dir=$(cd "$codegen_dir/../../../fortsym" && pwd)
actual_revision=$(git -C "$fortsym_dir" rev-parse HEAD)

if [[ "$actual_revision" != "$locked_revision" ]]; then
    echo "fortsym checkout does not match fortsym.lock" >&2
    echo "locked: $locked_revision" >&2
    echo "actual: $actual_revision" >&2
    exit 1
fi
