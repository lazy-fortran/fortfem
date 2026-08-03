#!/usr/bin/env bash
set -euo pipefail

repository_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
fixture_dir="$repository_dir/tools/fixtures/module_layer_audit"

# These sources are checker inputs, not package test sources.  Keeping the
# assertion explicit prevents a future fixture from being moved back below
# FPM's recursive test discovery root by accident.
if [[ "$fixture_dir" == "$repository_dir/test"/* ]]; then
    echo "module-layer fixtures must stay outside FPM's test discovery root" >&2
    exit 1
fi
test -f "$fixture_dir/bad_core_application/core/core.f90"
test -f "$fixture_dir/bad_core_application/plot/plot.f90"

# FPM's source graph is the behavioral oracle: before the relocation this
# command failed while resolving bad_plot_application.  Listing targets is
# enough to exercise discovery without compiling or running the full suite.
if ! command -v fpm >/dev/null 2>&1; then
    echo "fpm fixture-discovery oracle skipped: fpm is unavailable"
    exit 0
fi

timeout --foreground "${FPM_FIXTURE_DISCOVERY_TIMEOUT:-10}s" \
    fpm test --list >/dev/null 2>&1
echo "fpm test discovery excludes module-layer negative fixtures"
