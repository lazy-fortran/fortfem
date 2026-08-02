#!/usr/bin/env python3
"""Gate an optional exact Biro-data run without bundling paper arrays.

With no ``--data-manifest`` (the normal repository path), this command exits
successfully with an explicit ``SKIP``. A future sister-repository checkout
can provide a reviewed data manifest and archive; this adapter then verifies
the pinned URI, case, provenance, and SHA-256 before a downstream solver is
allowed to consume the artifact. It does not download, execute, or copy an
external solver.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

# Keep the adapter runnable directly from the repository root without turning
# ``tools`` into a Python package or adding a runtime dependency.
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tools"))

from validate_biro_external_manifest import (
    CASE_ID,
    ContractError,
    canonical_file_sha256,
    require,
    validate_external_data_manifest,
    validate_manifest,
)


def load_json(path: Path) -> dict:
    document = json.loads(path.read_text(encoding="utf-8"))
    require(isinstance(document, dict), f"{path} must contain a JSON object")
    return document


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument(
        "--contract",
        type=Path,
        default=Path(__file__).with_name("biro_paper_manifest.json"),
        help="repository adapter contract",
    )
    result.add_argument(
        "--data-manifest",
        type=Path,
        help="local manifest from the reviewed sister-repository artifact",
    )
    result.add_argument(
        "--data-root",
        type=Path,
        help="checkout root containing the artifact_file named by --data-manifest",
    )
    return result


def main() -> int:
    args = parser().parse_args()
    try:
        contract = load_json(args.contract)
        validate_manifest(contract)
        if args.data_manifest is None:
            status = contract["adapter"]["status"]
            if status == "skipped":
                print("SKIP: exact Biro paper data is absent; no external solver was run")
                return 0
            raise ContractError("ready contract requires --data-manifest")

        if not args.data_manifest.is_file():
            if contract["adapter"]["status"] == "skipped":
                print(f"SKIP: external data manifest is absent: {args.data_manifest}")
                return 0
            raise ContractError(f"external data manifest is absent: {args.data_manifest}")
        require(args.data_root is not None,
                "--data-root is required when --data-manifest is supplied")
        data_manifest = load_json(args.data_manifest)
        validate_external_data_manifest(data_manifest)
        require(contract["external_data"]["availability"] == "available",
                "contract has no pinned external artifact; exact run remains skipped")
        require(data_manifest["manifest_uri"] ==
                contract["external_data"]["manifest_uri"],
                "external data manifest URI does not match the contract")
        repository_uri = contract["external_data"]["repository_uri"].rstrip("/")
        require(data_manifest["artifact_uri"].startswith(repository_uri + "/"),
                "external artifact must come from the sister repository")
        require(data_manifest["archive_sha256"] ==
                contract["external_data"]["archive_sha256"],
                "external artifact checksum does not match the contract")
        artifact = (args.data_root / data_manifest["artifact_file"]).resolve()
        root = args.data_root.resolve()
        require(artifact == root or root in artifact.parents,
                "external artifact escapes --data-root")
        require(artifact.is_file(), f"external artifact is absent: {artifact}")
        actual = canonical_file_sha256(artifact)
        require(actual == data_manifest["archive_sha256"],
                "external artifact bytes do not match the manifest checksum")
        require(data_manifest["case_id"] == CASE_ID, "unexpected Biro case")
        print(f"READY: verified exact Biro artifact {data_manifest['artifact_uri']}")
        print(f"      sha256={actual}; downstream solver may consume the reviewed archive")
        return 0
    except (OSError, json.JSONDecodeError, ContractError) as error:
        print(f"Biro adapter: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
