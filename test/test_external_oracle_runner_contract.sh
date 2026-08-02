#!/usr/bin/env bash
set -euo pipefail

repository_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
manifest="$repository_dir/benchmark/external_oracles/runner_manifest.json"
validator="$repository_dir/tools/validate_external_oracle_runner.py"

python3 "$validator" "$manifest"

# The validator is also an independent behavioral oracle: a skipped runner
# without an explanation must not be publishable, even when every other field
# is present.  Keep this mutation in build/ so no fixture or external data is
# written to the repository.
invalid_manifest="$repository_dir/build/invalid_external_oracle_runner.json"
python3 - "$manifest" "$invalid_manifest" <<'PY'
import json
import pathlib
import sys

source = pathlib.Path(sys.argv[1])
target = pathlib.Path(sys.argv[2])
data = json.loads(source.read_text())
data["runners"][0]["skip_reason"] = ""
target.write_text(json.dumps(data, indent=2) + "\n")
PY

if python3 "$validator" "$invalid_manifest" >/dev/null 2>&1; then
    echo "invalid skipped runner unexpectedly validated" >&2
    exit 1
fi

echo "external oracle runner contract validates and rejects missing skip reasons"
