#!/usr/bin/env bash
set -euo pipefail

repository_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
contract="$repository_dir/benchmark/external_oracles/team_manifest.json"
adapter="$repository_dir/benchmark/external_oracles/run_team_adapter.py"
gallery="$repository_dir/benchmark/external_oracles/run_team_gallery.py"

python3 "$adapter" --contract "$contract" >/dev/null
skip_output=$(python3 "$gallery" --contract "$contract")
grep -Fq "SKIP: exact TEAM payloads are absent" <<<"$skip_output"

temporary_dir=$(mktemp -d)
trap 'rm -rf "$temporary_dir"' EXIT
artifact_name="synthetic-team-payload.json"

# The contract must reject a provenance edit, independently of repository state.
python3 - "$contract" "$temporary_dir/bad-contract.json" <<'PY'
import json
import pathlib
import sys

source = pathlib.Path(sys.argv[1])
target = pathlib.Path(sys.argv[2])
document = json.loads(source.read_text(encoding="utf-8"))
document["paper"]["provenance_uri"] = "https://example.invalid/not-team"
target.write_text(json.dumps(document, indent=2) + "\n", encoding="utf-8")
PY
if python3 "$adapter" --contract "$temporary_dir/bad-contract.json" >/dev/null 2>&1; then
    echo "wrong TEAM provenance unexpectedly validated" >&2
    exit 1
fi

# This is a deliberately tiny synthetic payload.  It exercises the checksum
# and renderer paths without copying any TEAM geometry, curves, or arrays.
python3 - "$contract" "$temporary_dir" <<'PY'
import hashlib
import json
import pathlib
import sys

contract_path = pathlib.Path(sys.argv[1])
root = pathlib.Path(sys.argv[2])
contract = json.loads(contract_path.read_text(encoding="utf-8"))
cases = {}
for case_id in ("TEAM-3", "TEAM-7", "TEAM-13", "TEAM-20"):
    cases[case_id] = {
        "provenance_uri": "https://www.osti.gov/biblio/7179128",
        "geometry": {
            "nodes": [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]],
            "elements": [[0, 1, 2]],
        },
        "solution": {"values": [0.0, 1.0, 2.0]},
    }
artifact_name = "synthetic-team-payload.json"
artifact_path = root / artifact_name
artifact_path.write_text(json.dumps({
    "schema": "fortfem-team-payload-1",
    "cases": cases,
}, indent=2) + "\n", encoding="utf-8")
digest = hashlib.sha256(artifact_path.read_bytes()).hexdigest()
manifest_uri = "https://github.com/lazy-fortran/fortfem-benchmark-data/blob/fixture/team.json"
artifact_uri = "https://github.com/lazy-fortran/fortfem-benchmark-data/blob/fixture/" + artifact_name
contract["external_data"].update({
    "manifest_uri": manifest_uri,
    "archive_sha256": digest,
    "availability": "available",
})
contract["adapter"].update({"status": "ready", "skip_reason": None})
contract["gallery"].update({"status": "ready", "skip_reason": None})
(root / "ready-contract.json").write_text(
    json.dumps(contract, indent=2) + "\n", encoding="utf-8")
(root / "data-manifest.json").write_text(json.dumps({
    "schema": "fortfem-team-data-1",
    "case_ids": ["TEAM-3", "TEAM-7", "TEAM-13", "TEAM-20"],
    "provenance_uri": "https://www.osti.gov/biblio/7179128",
    "manifest_uri": manifest_uri,
    "artifact_uri": artifact_uri,
    "archive_sha256": digest,
    "artifact_file": artifact_name,
    "license": "synthetic test fixture; not TEAM data",
    "payload_schema": "fortfem-team-payload-1",
    "members": {case_id: case_id for case_id in
                 ("TEAM-3", "TEAM-7", "TEAM-13", "TEAM-20")},
}, indent=2) + "\n", encoding="utf-8")
PY

ready_output=$(python3 "$adapter" \
    --contract "$temporary_dir/ready-contract.json" \
    --data-manifest "$temporary_dir/data-manifest.json" \
    --data-root "$temporary_dir")
grep -Fq "READY: verified exact TEAM artifact" <<<"$ready_output"

gallery_output=$(python3 "$gallery" \
    --contract "$temporary_dir/ready-contract.json" \
    --data-manifest "$temporary_dir/data-manifest.json" \
    --data-root "$temporary_dir" \
    --case-id TEAM-3 --output-dir "$temporary_dir/gallery")
grep -Fq "SOLUTION_FIRST: wrote verified TEAM gallery" <<<"$gallery_output"
test -s "$temporary_dir/gallery/TEAM-3/solution.svg"
test -s "$temporary_dir/gallery/TEAM-3/solution.csv"
grep -Fq "TEAM-3 exact-data solution" "$temporary_dir/gallery/TEAM-3/solution.svg"

# A checksum mismatch must stop consumption before the payload is rendered.
printf '\n' >> "$temporary_dir/$artifact_name"
if python3 "$adapter" \
    --contract "$temporary_dir/ready-contract.json" \
    --data-manifest "$temporary_dir/data-manifest.json" \
    --data-root "$temporary_dir" >/dev/null 2>&1; then
    echo "tampered TEAM artifact unexpectedly validated" >&2
    exit 1
fi

echo "TEAM external adapter validates provenance, schema, checksum, synthetic solution rendering, and tamper rejection"
