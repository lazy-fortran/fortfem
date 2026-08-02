#!/usr/bin/env bash
set -euo pipefail

repository_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
contract="$repository_dir/benchmark/external_oracles/biro_paper_manifest.json"
validator="$repository_dir/tools/validate_biro_external_manifest.py"
runner="$repository_dir/benchmark/external_oracles/run_biro_paper_adapter.py"

python3 "$validator" "$contract"
skip_output=$(python3 "$runner")
grep -Fq "SKIP: exact Biro paper data is absent" <<<"$skip_output"

temporary_dir=$(mktemp -d)
trap 'rm -rf "$temporary_dir"' EXIT
invalid_manifest="$temporary_dir/invalid.json"
python3 - "$contract" "$invalid_manifest" <<'PY'
import json
import pathlib
import sys

source = pathlib.Path(sys.argv[1])
target = pathlib.Path(sys.argv[2])
document = json.loads(source.read_text(encoding="utf-8"))
document["paper"]["provenance_uri"] = "https://example.invalid/not-the-paper"
target.write_text(json.dumps(document, indent=2) + "\n", encoding="utf-8")
PY
if python3 "$validator" "$invalid_manifest" >/dev/null 2>&1; then
    echo "invalid Biro provenance manifest unexpectedly validated" >&2
    exit 1
fi

# Exercise the READY path with a synthetic payload only. It is deliberately
# not geometry/material/source data from the paper; its sole purpose is to keep
# the checksum and path-consumption behavior independently covered.
python3 - "$contract" "$temporary_dir" <<'PY'
import hashlib
import json
import pathlib
import sys

contract_path = pathlib.Path(sys.argv[1])
root = pathlib.Path(sys.argv[2])
payload = b"synthetic adapter checksum fixture; not paper data\n"
artifact_name = "synthetic-reviewed-artifact.bin"
(root / artifact_name).write_bytes(payload)
digest = hashlib.sha256(payload).hexdigest()
contract = json.loads(contract_path.read_text(encoding="utf-8"))
contract["external_data"].update({
    "manifest_uri": "https://github.com/lazy-fortran/fortfem-benchmark-data/blob/fixture/biro-paper-data.json",
    "archive_sha256": digest,
    "availability": "available",
})
contract["adapter"].update({"status": "ready", "skip_reason": None})
(root / "ready-contract.json").write_text(
    json.dumps(contract, indent=2) + "\n", encoding="utf-8")
data_manifest = {
    "schema": "fortfem-biro-paper-data-1",
    "case_id": "biro-1996-tree-cotree-magnetostatic-v1",
    "provenance_uri": "https://doi.org/10.1109/20.497322",
    "manifest_uri": "https://github.com/lazy-fortran/fortfem-benchmark-data/blob/fixture/biro-paper-data.json",
    "artifact_uri": "https://github.com/lazy-fortran/fortfem-benchmark-data/blob/fixture/synthetic-reviewed-artifact.bin",
    "archive_sha256": digest,
    "artifact_file": artifact_name,
    "members": {"geometry": "geometry", "material": "material", "source": "source"},
}
(root / "data-manifest.json").write_text(
    json.dumps(data_manifest, indent=2) + "\n", encoding="utf-8")
PY
ready_output=$(python3 "$runner" \
    --contract "$temporary_dir/ready-contract.json" \
    --data-manifest "$temporary_dir/data-manifest.json" \
    --data-root "$temporary_dir")
grep -Fq "READY: verified exact Biro artifact" <<<"$ready_output"

echo "Biro external adapter validates, rejects wrong provenance, skips absent data, and verifies a pinned synthetic checksum"
