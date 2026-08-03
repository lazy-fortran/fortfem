#!/usr/bin/env bash
set -euo pipefail

repository_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "$repository_dir"

timeout --foreground "${API_CONSUMER_TEST_TIMEOUT:-10}s" \
    fo test test_api07_tensor_consumer

echo "API-07 tensor consumer smoke passed"
