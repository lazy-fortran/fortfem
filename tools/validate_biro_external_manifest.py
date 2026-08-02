#!/usr/bin/env python3
"""Validate the provenance boundary for the exact Biro paper adapter.

The repository deliberately contains no geometry, material, or source arrays
from Biro, Preis, and Richter. This validator only checks a metadata contract
which points at an independently licensed sister-repository artifact. A future
adapter may turn the ``absent`` record into an ``available`` record after the
artifact has been reviewed and its immutable checksum is known.
"""

from __future__ import annotations

import hashlib
import json
import re
import sys
from pathlib import Path
from typing import Any


SCHEMA = "fortfem-biro-paper-adapter-1"
DATA_SCHEMA = "fortfem-biro-paper-data-1"
PAYLOAD_SCHEMA = "fortfem-biro-paper-payload-1"
CASE_ID = "biro-1996-tree-cotree-magnetostatic-v1"
PAPER_URI = "https://doi.org/10.1109/20.497322"
ARRAY_KEYS = ("geometry", "material", "source")
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


def nullable_sha256(value: Any, label: str) -> None:
    require(value is None or
            (isinstance(value, str) and SHA256.fullmatch(value) is not None),
            f"{label} must be null or a lower-case SHA-256 digest")


def validate_manifest(document: dict[str, Any]) -> None:
    """Validate the checked-in adapter contract.

    ``arrays`` is intentionally metadata-only. An entry has a semantic name
    and an artifact member name, but it cannot contain numerical values. This
    prevents accidentally committing the paper's source arrays while still
    making the expected external payload explicit.
    """

    require(document.get("schema") == SCHEMA, f"schema must be {SCHEMA}")

    paper = document.get("paper")
    require(isinstance(paper, dict), "paper metadata is required")
    require(isinstance(paper.get("citation"), str) and paper["citation"].strip(),
            "paper citation is required")
    require(paper.get("provenance_uri") == PAPER_URI,
            f"paper provenance URI must be {PAPER_URI}")
    require(isinstance(paper.get("content_policy"), str) and
            paper["content_policy"].strip(),
            "paper content policy is required")

    case = document.get("case")
    require(isinstance(case, dict), "case metadata is required")
    require(case.get("id") == CASE_ID, f"case id must be {CASE_ID}")
    require(isinstance(case.get("method"), str) and case["method"].strip(),
            "case method is required")
    arrays = case.get("arrays")
    require(isinstance(arrays, dict), "case arrays metadata is required")
    require(tuple(arrays) == ARRAY_KEYS,
            "arrays must list geometry, material, then source")
    for name in ARRAY_KEYS:
        entry = arrays[name]
        require(isinstance(entry, dict), f"arrays.{name} must be an object")
        require(entry.get("state") == "external",
                f"arrays.{name} must remain external")
        require(entry.get("embedded") is False,
                f"arrays.{name} must not be embedded")
        require(entry.get("artifact_member") == name,
                f"arrays.{name}.artifact_member must be {name!r}")
        require(set(entry).issubset({"state", "embedded", "artifact_member",
                                     "description"}),
                f"arrays.{name} contains unsupported data fields")
        require(isinstance(entry.get("description"), str) and
                entry["description"].strip(),
                f"arrays.{name}.description is required")

    external = document.get("external_data")
    require(isinstance(external, dict), "external_data metadata is required")
    https_uri(external.get("repository_uri"),
              "external_data.repository_uri")
    availability = external.get("availability")
    require(availability in ("absent", "available"),
            "external_data.availability must be absent or available")
    manifest_uri = external.get("manifest_uri")
    archive_sha256 = external.get("archive_sha256")
    if availability == "absent":
        require(manifest_uri is None,
                "absent external data cannot publish a manifest URI")
        require(archive_sha256 is None,
                "absent external data cannot publish an archive checksum")
    else:
        https_uri(manifest_uri, "external_data.manifest_uri")
        nullable_sha256(archive_sha256, "external_data.archive_sha256")
        require(isinstance(archive_sha256, str),
                "available external data must pin an archive checksum")

    adapter = document.get("adapter")
    require(isinstance(adapter, dict), "adapter metadata is required")
    status = adapter.get("status")
    require(status in ("skipped", "ready"),
            "adapter.status must be skipped or ready")
    require(adapter.get("runner") ==
            "benchmark/external_oracles/run_biro_paper_adapter.py",
            "adapter.runner must identify the repository adapter")
    reason = adapter.get("skip_reason")
    if status == "skipped":
        require(isinstance(reason, str) and reason.strip(),
                "skipped adapter requires an explicit skip_reason")
        require(availability == "absent",
                "skipped adapter requires absent external data")
    else:
        require(reason is None, "ready adapter must not retain skip_reason")
        require(availability == "available",
                "ready adapter requires available external data")

    validate_gallery_metadata(document)


def validate_external_data_manifest(document: dict[str, Any]) -> None:
    """Validate a sister-repository artifact manifest without reading arrays."""

    require(document.get("schema") == DATA_SCHEMA,
            f"data schema must be {DATA_SCHEMA}")
    require(document.get("case_id") == CASE_ID,
            f"data case_id must be {CASE_ID}")
    require(document.get("provenance_uri") == PAPER_URI,
            f"data provenance URI must be {PAPER_URI}")
    https_uri(document.get("manifest_uri"), "data manifest_uri")
    https_uri(document.get("artifact_uri"), "data artifact_uri")
    digest = document.get("archive_sha256")
    require(isinstance(digest, str) and SHA256.fullmatch(digest) is not None,
            "data archive_sha256 must be a lower-case SHA-256 digest")
    member = document.get("artifact_file")
    require(isinstance(member, str) and member and
            not Path(member).is_absolute() and ".." not in Path(member).parts,
            "data artifact_file must be a relative safe path")
    members = document.get("members")
    require(isinstance(members, dict) and tuple(members) == ARRAY_KEYS,
            "data members must list geometry, material, then source")
    for name in ARRAY_KEYS:
        require(isinstance(members[name], str) and members[name].strip(),
                f"data members.{name} must be a non-empty artifact selector")


def validate_gallery_metadata(document: dict[str, Any]) -> None:
    """Validate the optional solution-first gallery contract."""

    gallery = document.get("gallery")
    require(isinstance(gallery, dict), "gallery metadata is required")
    require(gallery.get("runner") ==
            "benchmark/external_oracles/run_biro_paper_gallery.py",
            "gallery.runner must identify the solution-first gallery runner")
    require(gallery.get("payload_schema") == PAYLOAD_SCHEMA,
            f"gallery.payload_schema must be {PAYLOAD_SCHEMA}")
    require(gallery.get("plot") == "solution-first-svg",
            "gallery.plot must be solution-first-svg")
    status = gallery.get("status")
    require(status in ("skipped", "ready"),
            "gallery.status must be skipped or ready")
    if status == "skipped":
        require(isinstance(gallery.get("skip_reason"), str) and
                gallery["skip_reason"].strip(),
                "skipped gallery requires an explicit skip_reason")
    else:
        require(gallery.get("skip_reason") is None,
                "ready gallery must not retain skip_reason")


def canonical_file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> int:
    if len(sys.argv) != 2:
        print(f"usage: {Path(sys.argv[0]).name} MANIFEST.json", file=sys.stderr)
        return 2
    path = Path(sys.argv[1])
    try:
        document = json.loads(path.read_text(encoding="utf-8"))
        require(isinstance(document, dict), "manifest root must be an object")
        validate_manifest(document)
    except (OSError, json.JSONDecodeError, ContractError) as error:
        print(f"Biro adapter contract: {error}", file=sys.stderr)
        return 1
    print(f"Biro adapter contract valid: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
