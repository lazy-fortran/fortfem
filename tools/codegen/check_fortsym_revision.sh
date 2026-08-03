#!/usr/bin/env bash
set -euo pipefail

codegen_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
locked_revision=$(tr -d '[:space:]' < "$codegen_dir/fortsym.lock")
fortsym_dir=${FORTSYM_DIR:-"$codegen_dir/../../../fortsym"}
if [[ ! -d "$fortsym_dir" ]]; then
    echo "FortSym checkout does not exist: $fortsym_dir" >&2
    exit 1
fi
fortsym_dir=$(cd "$fortsym_dir" && pwd)
actual_revision=$(git -C "$fortsym_dir" rev-parse HEAD)

if [[ "$actual_revision" != "$locked_revision" ]]; then
    echo "fortsym checkout does not match fortsym.lock" >&2
    echo "locked: $locked_revision" >&2
    echo "actual: $actual_revision" >&2
    exit 1
fi

if [[ -n "$(git -C "$fortsym_dir" status --porcelain)" ]]; then
    echo "FortSym checkout has uncommitted changes: $fortsym_dir" >&2
    exit 1
fi
