#!/usr/bin/env python3
"""Validate optional exact TEAM benchmark data without bundling arrays.

The checked-in contract is metadata only.  Exact TEAM meshes, material curves,
excitations, and observations must be supplied by a separately licensed sister
repository.  The adapter verifies the pinned manifest, provenance, license
declaration, and archive bytes; it never downloads or invokes an external
solver.  A missing artifact is an explicit successful ``SKIP``.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from pathlib import Path
from typing import Any


SCHEMA = "fortfem-team-external-adapter-1"
DATA_SCHEMA = "fortfem-team-data-1"
PAYLOAD_SCHEMA = "fortfem-team-payload-1"
CONTRACT_CASES = ("TEAM-3", "TEAM-7", "TEAM-13", "TEAM-20")
PAPER_URI = "https://www.osti.gov/biblio/7179128"
SHA256 = re.compile(r"^[0-9a-f]{64}$")


class ContractError(ValueError):
    """A human-readable provenance-contract violation."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ContractError(message)


def https_uri(value: Any, label: str) -> None:
    require(isinstance(value, str) and value.startswith("https://") and
            value.strip() == value and len(value) > len("https://"),
            f"{label} must be an HTTPS URI")


def safe_relative_path(value: Any, label: str) -> None:
    require(isinstance(value, str) and value and
            not Path(value).is_absolute() and ".." not in Path(value).parts,
            f"{label} must be a relative safe path")


def canonical_file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    document = json.loads(path.read_text(encoding="utf-8"))
    require(isinstance(document, dict), f"{path} must contain a JSON object")
    return document


def case_map(document: dict[str, Any]) -> dict[str, dict[str, Any]]:
    cases = document.get("cases")
    require(isinstance(cases, list) and
            tuple(case.get("id") for case in cases if isinstance(case, dict)) ==
            CONTRACT_CASES, "cases must list TEAM-3, TEAM-7, TEAM-13, TEAM-20")
    result: dict[str, dict[str, Any]] = {}
    for case in cases:
        require(isinstance(case.get("id"), str), "case id is required")
        require(case["id"] not in result, f"duplicate case {case['id']}")
        require(isinstance(case.get("title"), str) and case["title"].strip(),
                f"{case['id']} title is required")
        links = case.get("provenance_links")
        require(isinstance(links, list) and links,
                f"{case['id']} provenance_links are required")
        for link in links:
            https_uri(link, f"{case['id']} provenance link")
        require(case.get("arrays") ==
                ["geometry", "material", "source", "observables"],
                f"{case['id']} arrays must remain metadata-only")
        result[case["id"]] = case
    return result


def validate_manifest(document: dict[str, Any]) -> None:
    require(document.get("schema") == SCHEMA, f"schema must be {SCHEMA}")
    require(document.get("benchmark_family") == "TEAM",
            "benchmark_family must be TEAM")
    paper = document.get("paper")
    require(isinstance(paper, dict), "paper metadata is required")
    require(isinstance(paper.get("citation"), str) and paper["citation"].strip(),
            "paper citation is required")
    require(paper.get("provenance_uri") == PAPER_URI,
            f"paper provenance URI must be {PAPER_URI}")
    links = paper.get("provenance_links")
    require(isinstance(links, list) and links, "paper provenance_links are required")
    for link in links:
        https_uri(link, "paper provenance link")
    require(isinstance(paper.get("content_policy"), str) and
            paper["content_policy"].strip(), "paper content policy is required")
    case_map(document)

    external = document.get("external_data")
    require(isinstance(external, dict), "external_data metadata is required")
    https_uri(external.get("repository_uri"), "external_data.repository_uri")
    require(external.get("payload_schema") == PAYLOAD_SCHEMA,
            f"external_data.payload_schema must be {PAYLOAD_SCHEMA}")
    require(isinstance(external.get("license"), str) and
            external["license"].strip(), "external_data.license is required")
    availability = external.get("availability")
    require(availability in ("absent", "available"),
            "external_data.availability must be absent or available")
    if availability == "absent":
        require(external.get("manifest_uri") is None,
                "absent external data cannot publish a manifest URI")
        require(external.get("archive_sha256") is None,
                "absent external data cannot publish an archive checksum")
    else:
        https_uri(external.get("manifest_uri"), "external_data.manifest_uri")
        digest = external.get("archive_sha256")
        require(isinstance(digest, str) and SHA256.fullmatch(digest) is not None,
                "available external data must pin a lower-case SHA-256 digest")

    adapter = document.get("adapter")
    require(isinstance(adapter, dict), "adapter metadata is required")
    require(adapter.get("runner") ==
            "benchmark/external_oracles/run_team_adapter.py",
            "adapter.runner must identify the repository adapter")
    status = adapter.get("status")
    require(status in ("skipped", "ready"),
            "adapter.status must be skipped or ready")
    if status == "skipped":
        require(external["availability"] == "absent",
                "skipped adapter requires absent external data")
        require(isinstance(adapter.get("skip_reason"), str) and
                adapter["skip_reason"].strip(),
                "skipped adapter requires an explicit skip_reason")
    else:
        require(external["availability"] == "available",
                "ready adapter requires available external data")
        require(adapter.get("skip_reason") is None,
                "ready adapter must not retain skip_reason")

    gallery = document.get("gallery")
    require(isinstance(gallery, dict), "gallery metadata is required")
    require(gallery.get("runner") ==
            "benchmark/external_oracles/run_team_gallery.py",
            "gallery.runner must identify the solution-first gallery")
    require(gallery.get("payload_schema") == PAYLOAD_SCHEMA,
            f"gallery.payload_schema must be {PAYLOAD_SCHEMA}")
    require(gallery.get("plot") == "solution-first-svg",
            "gallery.plot must be solution-first-svg")
    require(gallery.get("status") in ("skipped", "ready"),
            "gallery.status must be skipped or ready")
    if gallery["status"] == "skipped":
        require(isinstance(gallery.get("skip_reason"), str) and
                gallery["skip_reason"].strip(),
                "skipped gallery requires an explicit skip_reason")
    else:
        require(gallery.get("skip_reason") is None,
                "ready gallery must not retain skip_reason")


def validate_external_data_manifest(document: dict[str, Any]) -> None:
    require(document.get("schema") == DATA_SCHEMA,
            f"data schema must be {DATA_SCHEMA}")
    case_ids = document.get("case_ids")
    require(case_ids == list(CONTRACT_CASES),
            "data case_ids must list all TEAM cases in contract order")
    require(document.get("provenance_uri") == PAPER_URI,
            f"data provenance URI must be {PAPER_URI}")
    https_uri(document.get("manifest_uri"), "data manifest_uri")
    https_uri(document.get("artifact_uri"), "data artifact_uri")
    digest = document.get("archive_sha256")
    require(isinstance(digest, str) and SHA256.fullmatch(digest) is not None,
            "data archive_sha256 must be a lower-case SHA-256 digest")
    safe_relative_path(document.get("artifact_file"), "data artifact_file")
    require(isinstance(document.get("license"), str) and
            document["license"].strip(), "data license is required")
    require(document.get("payload_schema") == PAYLOAD_SCHEMA,
            f"data payload_schema must be {PAYLOAD_SCHEMA}")
    members = document.get("members")
    require(isinstance(members, dict) and tuple(members) == CONTRACT_CASES,
            "data members must list all TEAM cases")
    for case_id in CONTRACT_CASES:
        require(isinstance(members[case_id], str) and members[case_id].strip(),
                f"data members.{case_id} must be a non-empty selector")


def locate_verified_artifact(contract: dict[str, Any], data_manifest_path: Path | None,
                             data_root: Path | None) -> tuple[Path | None, dict[str, Any] | None]:
    if data_manifest_path is None:
        require(contract["adapter"]["status"] == "skipped",
                "ready contract requires --data-manifest")
        return None, None
    if not data_manifest_path.is_file():
        require(contract["adapter"]["status"] == "skipped",
                f"external data manifest is absent: {data_manifest_path}")
        return None, None
    require(data_root is not None, "--data-root is required with --data-manifest")
    data_manifest = load_json(data_manifest_path)
    validate_external_data_manifest(data_manifest)
    external = contract["external_data"]
    require(external["availability"] == "available",
            "contract has no pinned external artifact; exact run remains skipped")
    require(data_manifest["manifest_uri"] == external["manifest_uri"],
            "external data manifest URI does not match contract")
    require(data_manifest["archive_sha256"] == external["archive_sha256"],
            "external artifact checksum does not match contract")
    repository_uri = external["repository_uri"].rstrip("/")
    require(data_manifest["artifact_uri"].startswith(repository_uri + "/"),
            "external artifact must come from the sister repository")
    artifact = (data_root / data_manifest["artifact_file"]).resolve()
    root = data_root.resolve()
    require(root in artifact.parents, "external artifact escapes --data-root")
    require(artifact.is_file(), f"external artifact is absent: {artifact}")
    require(canonical_file_sha256(artifact) == data_manifest["archive_sha256"],
            "external artifact bytes do not match the manifest checksum")
    return artifact, data_manifest


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--contract", type=Path,
                        default=Path(__file__).with_name("team_manifest.json"))
    result.add_argument("--data-manifest", type=Path)
    result.add_argument("--data-root", type=Path)
    return result


def main() -> int:
    args = parser().parse_args()
    try:
        contract = load_json(args.contract)
        validate_manifest(contract)
        artifact, _ = locate_verified_artifact(contract, args.data_manifest, args.data_root)
        if artifact is None:
            print("SKIP: exact TEAM data is absent; no external solver was run")
            return 0
        print(f"READY: verified exact TEAM artifact {artifact.name}")
        print(f"       sha256={canonical_file_sha256(artifact)}")
        return 0
    except (OSError, json.JSONDecodeError, ContractError) as error:
        print(f"TEAM adapter: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
