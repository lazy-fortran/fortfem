#!/usr/bin/env bash
set -euo pipefail

repository_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "$repository_dir"

# Keep this downstream-style API smoke bounded independently of the full suite.
timeout --foreground "${API_CONSUMER_TEST_TIMEOUT:-10}s" \
    fo test test_api07_free_boundary_consumer

echo "API-07 free-boundary consumer smoke passed"
