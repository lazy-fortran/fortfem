#!/usr/bin/env bash
set -euo pipefail

repository_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
python3 - "$repository_dir" <<'PY'
import json
import pathlib
import sys

root = pathlib.Path(sys.argv[1])
expected = {
    "team3_neutral_benchmark": "TEAM-3",
    "team7_neutral_benchmark": "TEAM-7",
    "team13_neutral_benchmark": "TEAM-13",
    "team20_neutral_benchmark": "TEAM-20",
}
for fixture, label in expected.items():
    path = root / "example" / fixture / "provenance.json"
    document = json.loads(path.read_text(encoding="utf-8"))
    assert document["fixture"] == fixture
    assert document["benchmark_family"] == "TEAM"
    assert document["exact_reproduction"] is False
    assert document["source_arrays_redistributed"] is False
    assert document["references"]
    assert all(reference["url"].startswith("https://")
               for reference in document["references"])
    assert label in document["shape_target"]
    assert "external" in document["license_note"]
print("TEAM neutral fixtures have explicit provenance and non-exact status")
PY
