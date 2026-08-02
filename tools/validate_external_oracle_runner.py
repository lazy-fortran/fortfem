#!/usr/bin/env python3
"""Validate the optional external-code oracle runner contract.

The contract deliberately uses only the Python standard library.  It validates
metadata and a nine-point analytical Poisson reference; it never imports,
executes, or downloads FreeFEM, MFEM, or FEniCSx.  Real runner records can be
published in the separate benchmark-data repository after this check passes.
"""

from __future__ import annotations

import hashlib
import json
import math
import re
import sys
from pathlib import Path
from typing import Any


SCHEMA = "fortfem-external-oracle-runner-1"
TARGETS = ("freefem", "mfem", "fenicsx")
SHA256 = re.compile(r"^[0-9a-f]{64}$")


class ContractError(ValueError):
    """A human-readable contract violation."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ContractError(message)


def finite_nonnegative(value: Any, label: str) -> None:
    require(isinstance(value, (int, float)) and not isinstance(value, bool),
            f"{label} must be numeric")
    require(math.isfinite(float(value)) and float(value) >= 0.0,
            f"{label} must be finite and non-negative")


def canonical_digest(value: Any) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(encoded.encode("utf-8")).hexdigest()


def validate_case(case: dict[str, Any]) -> None:
    require(case.get("id") == "unit-square-poisson-mms-v1",
            "case id must identify the shared Poisson MMS")
    require(case.get("equation") == "-laplacian(u) = f",
            "case equation is not the analytical Poisson contract")
    require(case.get("exact_solution") == "sin(pi*x)*sin(pi*y)",
            "case exact solution is not the analytical Poisson contract")
    require(case.get("source") == "2*pi^2*sin(pi*x)*sin(pi*y)",
            "case source is not the analytical Poisson contract")
    sampling = case.get("sampling")
    require(isinstance(sampling, dict), "case sampling must be an object")
    require(sampling.get("coordinate_system") == "physical-cartesian",
            "case sampling must use physical Cartesian coordinates")
    points = sampling.get("points")
    values = sampling.get("reference_values")
    weights = sampling.get("weights")
    require(isinstance(points, list) and isinstance(values, list),
            "analytical sampling points and values are required")
    require(isinstance(weights, list), "analytical sampling weights are required")
    require(len(points) == len(values) == len(weights) == 9,
            "the Poisson metadata fixture must contain nine samples")
    require(sampling.get("sample_count") == len(points),
            "sample_count must agree with the analytical sample list")
    require(all(isinstance(point, list) and len(point) == 2 for point in points),
            "every analytical point must be a two-dimensional coordinate")
    require(all(0.0 <= float(x) <= 1.0 and 0.0 <= float(y) <= 1.0
                for x, y in points),
            "analytical points must lie in the unit square")
    require(all(isinstance(value, (int, float)) for value in values),
            "analytical reference values must be numeric")
    require(all(isinstance(weight, (int, float)) and float(weight) > 0.0
                for weight in weights),
            "analytical sampling weights must be positive")
    require(math.isclose(sum(float(weight) for weight in weights), 1.0,
                         rel_tol=0.0, abs_tol=1.0e-14),
            "analytical sampling weights must sum to one")
    for point, value in zip(points, values):
        x, y = (float(component) for component in point)
        expected = math.sin(math.pi * x) * math.sin(math.pi * y)
        require(math.isclose(float(value), expected, rel_tol=0.0,
                             abs_tol=1.0e-14),
                f"analytical Poisson sample mismatch at ({x}, {y})")
    reference_artifact = sampling.get("reference_artifact")
    require(isinstance(reference_artifact, dict),
            "reference artifact provenance is required")
    uri = reference_artifact.get("repository_uri")
    require(isinstance(uri, str) and uri.startswith("https://"),
            "reference artifact must point to a sister repository URI")
    coordinate_digest = reference_artifact.get("coordinates_sha256")
    value_digest = reference_artifact.get("values_sha256")
    require(isinstance(coordinate_digest, str) and SHA256.fullmatch(coordinate_digest),
            "reference coordinate checksum must be SHA-256")
    require(isinstance(value_digest, str) and SHA256.fullmatch(value_digest),
            "reference value checksum must be SHA-256")
    require(coordinate_digest == canonical_digest(points),
            "reference coordinate checksum does not match points")
    require(value_digest == canonical_digest(values),
            "reference value checksum does not match values")
    tolerances = case.get("tolerances")
    require(isinstance(tolerances, dict), "case tolerances must be an object")
    for name in ("coordinate", "absolute", "relative", "residual"):
        finite_nonnegative(tolerances.get(name), f"case tolerances.{name}")


def validate_runner(runner: dict[str, Any], case: dict[str, Any]) -> None:
    runner_id = runner.get("id")
    require(runner_id in TARGETS, f"unexpected runner id: {runner_id!r}")
    require(isinstance(runner.get("code"), str) and runner["code"].strip(),
            f"{runner_id}: code name is required")
    require(isinstance(runner.get("license"), str) and runner["license"].strip(),
            f"{runner_id}: license is required")
    require(isinstance(runner.get("source_uri"), str) and
            runner["source_uri"].startswith("https://"),
            f"{runner_id}: source URI is required")
    launch = runner.get("launch")
    require(isinstance(launch, dict), f"{runner_id}: launch metadata is required")
    command = launch.get("command")
    require(isinstance(command, list) and command and
            all(isinstance(item, str) and item for item in command),
            f"{runner_id}: launch command must be a non-empty argument list")
    require(isinstance(launch.get("working_directory"), str) and
            launch["working_directory"],
            f"{runner_id}: working directory is required")
    finite_nonnegative(launch.get("timeout_seconds"),
                       f"{runner_id}: launch timeout_seconds")
    require(float(launch["timeout_seconds"]) > 0.0,
            f"{runner_id}: launch timeout_seconds must be positive")
    revision = runner.get("revision")
    require(revision is None or (isinstance(revision, str) and revision.strip()),
            f"{runner_id}: revision must be null or a non-empty immutable label")
    uv = runner.get("uv")
    require(isinstance(uv, dict) and isinstance(uv.get("supported"), bool),
            f"{runner_id}: uv support declaration is required")
    if uv["supported"]:
        uv_command = uv.get("postprocess_command")
        require(isinstance(uv_command, list) and uv_command and
                uv_command[0] == "uv" and
                all(isinstance(item, str) and item for item in uv_command),
                f"{runner_id}: uv postprocess command must start with uv")
        require(isinstance(uv.get("purpose"), str) and uv["purpose"].strip(),
                f"{runner_id}: uv purpose is required")
    else:
        require(isinstance(uv.get("reason"), str) and uv["reason"].strip(),
                f"{runner_id}: reason is required when uv is unavailable")

    sampling = runner.get("sampling")
    require(isinstance(sampling, dict), f"{runner_id}: sampling metadata is required")
    for name in ("coordinate_system", "rule"):
        require(isinstance(sampling.get(name), str) and sampling[name].strip(),
                f"{runner_id}: sampling.{name} is required")
    require(sampling.get("sample_count") == case["sampling"]["sample_count"],
            f"{runner_id}: sample count must match the shared case")

    tolerances = runner.get("tolerances")
    require(isinstance(tolerances, dict), f"{runner_id}: tolerances are required")
    for name in ("coordinate", "absolute", "relative", "residual"):
        finite_nonnegative(tolerances.get(name), f"{runner_id}: tolerances.{name}")

    timing = runner.get("timing")
    require(isinstance(timing, dict), f"{runner_id}: timing metadata is required")
    status = runner.get("status")
    require(status in ("passed", "skipped"),
            f"{runner_id}: status must be passed or skipped")
    phase_names = ("mesh_seconds", "assembly_seconds", "factorization_seconds",
                   "solve_seconds", "output_seconds")
    if status == "passed":
        require(isinstance(revision, str) and revision.strip(),
                f"{runner_id}: passed runner requires an immutable revision")
        for name in (*phase_names, "total_seconds", "peak_memory_bytes"):
            finite_nonnegative(timing.get(name), f"{runner_id}: timing.{name}")
        require(float(timing["total_seconds"]) >= sum(float(timing[name])
                    for name in phase_names),
                f"{runner_id}: total timing must cover all phases")
        for name in ("warmup_count", "repetition_count"):
            require(isinstance(timing.get(name), int) and timing[name] >= 0,
                    f"{runner_id}: timing.{name} must be a non-negative integer")
        require(timing["repetition_count"] >= 1,
                f"{runner_id}: repetition_count must be positive")
    else:
        require(isinstance(runner.get("skip_reason"), str) and
                runner["skip_reason"].strip(),
                f"{runner_id}: skipped runner requires an explicit skip_reason")
        for name in (*phase_names, "total_seconds", "peak_memory_bytes"):
            require(timing.get(name) is None,
                    f"{runner_id}: skipped timing.{name} must be null")
        require(timing.get("warmup_count") is None and
                timing.get("repetition_count") is None,
                f"{runner_id}: skipped repetition metadata must be null")

    artifact = runner.get("artifact")
    require(isinstance(artifact, dict), f"{runner_id}: artifact metadata is required")
    artifact_uri = artifact.get("repository_uri")
    require(isinstance(artifact_uri, str) and artifact_uri.startswith("https://"),
            f"{runner_id}: artifact repository URI is required")
    for name in ("manifest_uri", "coordinates_sha256", "samples_sha256"):
        require(name in artifact, f"{runner_id}: artifact.{name} key is required")
    if status == "passed":
        require(isinstance(artifact["manifest_uri"], str) and
                artifact["manifest_uri"].startswith("https://"),
                f"{runner_id}: passed artifact manifest URI is required")
        for name in ("coordinates_sha256", "samples_sha256"):
            require(isinstance(artifact[name], str) and SHA256.fullmatch(artifact[name]),
                    f"{runner_id}: passed artifact.{name} must be SHA-256")
    else:
        require(artifact["manifest_uri"] is None and
                artifact["coordinates_sha256"] is None and
                artifact["samples_sha256"] is None,
                f"{runner_id}: skipped artifact checksums must be null")


def validate_manifest(document: dict[str, Any]) -> None:
    require(document.get("schema") == SCHEMA,
            f"schema must be {SCHEMA}")
    require(document.get("data_repository_uri", "").startswith("https://"),
            "data_repository_uri must point to the sister benchmark repository")
    case = document.get("case")
    require(isinstance(case, dict), "case object is required")
    validate_case(case)
    runners = document.get("runners")
    require(isinstance(runners, list) and len(runners) == len(TARGETS),
            "the contract must list exactly the three optional runners")
    ids = [runner.get("id") for runner in runners if isinstance(runner, dict)]
    require(tuple(ids) == TARGETS,
            "runners must be ordered FreeFEM, MFEM, then FEniCSx")
    for runner in runners:
        require(isinstance(runner, dict), "runner entries must be objects")
        validate_runner(runner, case)


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
        print(f"external oracle runner contract: {error}", file=sys.stderr)
        return 1
    print(f"external oracle runner contract valid: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
