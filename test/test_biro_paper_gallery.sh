#!/usr/bin/env bash
set -euo pipefail

repository_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
contract="$repository_dir/benchmark/external_oracles/biro_paper_manifest.json"
validator="$repository_dir/tools/validate_biro_external_manifest.py"
gallery="$repository_dir/benchmark/external_oracles/run_biro_paper_gallery.py"

python3 "$validator" "$contract"
skip_output=$(python3 "$gallery")
grep -Fq "SKIP: exact Biro paper payload is absent" <<<"$skip_output"

temporary_dir=$(mktemp -d)
trap 'rm -rf "$temporary_dir"' EXIT

python3 - "$contract" "$temporary_dir" <<'PY'
import hashlib
import json
import pathlib
import sys

contract_path = pathlib.Path(sys.argv[1])
root = pathlib.Path(sys.argv[2])
payload = {
    "schema": "fortfem-biro-paper-payload-1",
    "case_id": "biro-1996-tree-cotree-magnetostatic-v1",
    "provenance_uri": "https://doi.org/10.1109/20.497322",
    "geometry": {
        "nodes": [[0.0, 0.0, 0.0], [1.0, 0.0, 0.2],
                   [1.0, 1.0, 0.8], [0.0, 1.0, 0.4]],
        "elements": [[0, 1, 2], [0, 2, 3]],
    },
    "material": {"member": "material"},
    "source": {"member": "source"},
    "solution": {"name": "vector-potential magnitude", "values": [
        [0.0, 0.0, 0.0], [1.0, 0.0, 0.0],
        [0.6, 0.8, 0.0], [0.0, 0.3, 0.4]
    ]},
}
artifact_name = "synthetic-reviewed-biro-payload.json"
artifact = root / artifact_name
artifact.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
digest = hashlib.sha256(artifact.read_bytes()).hexdigest()
contract = json.loads(contract_path.read_text(encoding="utf-8"))
manifest_uri = "https://github.com/lazy-fortran/fortfem-benchmark-data/blob/fixture/biro-gallery.json"
contract["external_data"].update({
    "manifest_uri": manifest_uri,
    "archive_sha256": digest,
    "availability": "available",
})
contract["adapter"].update({"status": "ready", "skip_reason": None})
contract["gallery"].update({"status": "ready", "skip_reason": None})
(root / "ready-contract.json").write_text(
    json.dumps(contract, indent=2) + "\n", encoding="utf-8")
data_manifest = {
    "schema": "fortfem-biro-paper-data-1",
    "case_id": "biro-1996-tree-cotree-magnetostatic-v1",
    "provenance_uri": "https://doi.org/10.1109/20.497322",
    "manifest_uri": manifest_uri,
    "artifact_uri": "https://github.com/lazy-fortran/fortfem-benchmark-data/blob/fixture/" + artifact_name,
    "archive_sha256": digest,
    "artifact_file": artifact_name,
    "members": {"geometry": "geometry", "material": "material", "source": "source"},
}
(root / "data-manifest.json").write_text(
    json.dumps(data_manifest, indent=2) + "\n", encoding="utf-8")
PY

output_dir="$temporary_dir/gallery"
ready_output=$(python3 "$gallery" \
    --contract "$temporary_dir/ready-contract.json" \
    --data-manifest "$temporary_dir/data-manifest.json" \
    --data-root "$temporary_dir" \
    --output-dir "$output_dir")
grep -Fq "SOLUTION_FIRST: wrote verified exact Biro gallery" <<<"$ready_output"
test -s "$output_dir/solution.svg"
test -s "$output_dir/solution.csv"
grep -Fq "solution-first" "$output_dir/solution.svg"
grep -Fq "<polygon" "$output_dir/solution.svg"
python3 - "$output_dir/provenance.json" <<'PY'
import json
import pathlib
import sys

metadata = json.loads(pathlib.Path(sys.argv[1]).read_text(encoding="utf-8"))
assert metadata["status"] == "ready"
assert metadata["exact_data"] is True
assert metadata["manufactured"] is False
assert metadata["comparison_mode"] == "external-benchmark-payload"
assert metadata["analytical_reference"].startswith("provenance-only")
assert metadata["solution_first"] is True
PY

python3 - "$output_dir/solution.csv" "$output_dir/solution.svg" <<'PY'
import csv
import pathlib
import sys

rows = list(csv.DictReader(pathlib.Path(sys.argv[1]).open(encoding="utf-8")))
assert [float(row["solution"]) for row in rows] == [0.0, 1.0, 1.0, 0.5]
svg = pathlib.Path(sys.argv[2]).read_text(encoding="utf-8")
assert svg.count("<polygon ") == 2
assert svg.count("<circle ") == 4
assert "solution-first" in svg
PY

# A provenance-valid payload still must describe drawable elements.  Reject a
# duplicate-node element before it can produce a misleading solution plot.
cp "$temporary_dir/synthetic-reviewed-biro-payload.json" \
   "$temporary_dir/degenerate-biro-payload.json"
python3 - "$temporary_dir/degenerate-biro-payload.json" <<'PY'
import json
import pathlib
import sys

path = pathlib.Path(sys.argv[1])
document = json.loads(path.read_text(encoding="utf-8"))
document["geometry"]["elements"] = [[0, 0, 1]]
path.write_text(json.dumps(document, indent=2) + "\n", encoding="utf-8")
PY
python3 - "$temporary_dir/data-manifest.json" \
    "$temporary_dir/degenerate-data-manifest.json" \
    "$temporary_dir/degenerate-biro-payload.json" <<'PY'
import hashlib
import json
import pathlib
import sys

source = pathlib.Path(sys.argv[1])
target = pathlib.Path(sys.argv[2])
artifact = pathlib.Path(sys.argv[3])
document = json.loads(source.read_text(encoding="utf-8"))
document["artifact_file"] = artifact.name
document["artifact_uri"] = document["artifact_uri"].replace(
    "synthetic-reviewed-biro-payload.json", artifact.name)
document["archive_sha256"] = hashlib.sha256(artifact.read_bytes()).hexdigest()
target.write_text(json.dumps(document, indent=2) + "\n", encoding="utf-8")
PY
if python3 "$gallery" \
    --contract "$temporary_dir/ready-contract.json" \
    --data-manifest "$temporary_dir/degenerate-data-manifest.json" \
    --data-root "$temporary_dir" \
    --output-dir "$temporary_dir/degenerate" >/dev/null 2>&1; then
    echo "degenerate Biro element unexpectedly rendered" >&2
    exit 1
fi

# A changed payload must be rejected before any plot is emitted.
printf 'tampered\n' >> "$temporary_dir/synthetic-reviewed-biro-payload.json"
if python3 "$gallery" \
    --contract "$temporary_dir/ready-contract.json" \
    --data-manifest "$temporary_dir/data-manifest.json" \
    --data-root "$temporary_dir" \
    --output-dir "$temporary_dir/tampered" >/dev/null 2>&1; then
    echo "tampered Biro payload unexpectedly rendered" >&2
    exit 1
fi

echo "Biro gallery skips absent data, renders verified solution-first SVG, and rejects tampering"
