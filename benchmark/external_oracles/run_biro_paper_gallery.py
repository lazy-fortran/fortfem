#!/usr/bin/env python3
"""Render a verified, optional Biro-paper solution payload.

The repository intentionally ships no paper-specific arrays.  With the
checked-in contract the command therefore exits successfully with ``SKIP``.
An independently licensed sister-repository checkout may provide a JSON
payload (``fortfem-biro-paper-payload-1``); after the adapter verifies its
provenance and bytes, this runner writes a solution-first SVG and CSV.  The
renderer is standard-library-only so the provenance gate can be exercised in
CI without a plotting or FEM dependency.
"""

from __future__ import annotations

import argparse
import html
import json
import math
import sys
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tools"))

from validate_biro_external_manifest import (  # noqa: E402
    CASE_ID,
    PAPER_URI,
    PAYLOAD_SCHEMA,
    ContractError,
    canonical_file_sha256,
    require,
    validate_external_data_manifest,
    validate_manifest,
)


def load_json(path: Path) -> dict[str, Any]:
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
    result.add_argument(
        "--output-dir",
        type=Path,
        default=Path("output/external_oracle/biro_paper_gallery"),
        help="directory for the solution-first SVG and CSV",
    )
    return result


def verified_artifact(
    contract: dict[str, Any], data_manifest_path: Path | None,
    data_root: Path | None,
) -> tuple[Path | None, dict[str, Any] | None]:
    """Return the checked artifact path and manifest, or the explicit skip."""

    if data_manifest_path is None:
        require(contract["adapter"]["status"] == "skipped",
                "ready contract requires --data-manifest")
        return None, None
    if not data_manifest_path.is_file():
        require(contract["adapter"]["status"] == "skipped",
                f"external data manifest is absent: {data_manifest_path}")
        return None, None

    require(data_root is not None,
            "--data-root is required when --data-manifest is supplied")
    data_manifest = load_json(data_manifest_path)
    validate_external_data_manifest(data_manifest)
    external = contract["external_data"]
    require(external["availability"] == "available",
            "contract has no pinned external artifact; exact gallery remains skipped")
    require(data_manifest["manifest_uri"] == external["manifest_uri"],
            "external data manifest URI does not match the contract")
    repository_uri = external["repository_uri"].rstrip("/")
    require(data_manifest["artifact_uri"].startswith(repository_uri + "/"),
            "external artifact must come from the sister repository")
    require(data_manifest["archive_sha256"] == external["archive_sha256"],
            "external artifact checksum does not match the contract")
    artifact = (data_root / data_manifest["artifact_file"]).resolve()
    root = data_root.resolve()
    require(artifact == root or root in artifact.parents,
            "external artifact escapes --data-root")
    require(artifact.is_file(), f"external artifact is absent: {artifact}")
    actual = canonical_file_sha256(artifact)
    require(actual == data_manifest["archive_sha256"],
            "external artifact bytes do not match the manifest checksum")
    require(data_manifest["case_id"] == CASE_ID, "unexpected Biro case")
    return artifact, data_manifest


def finite_number(value: Any, label: str) -> float:
    require(isinstance(value, (int, float)) and not isinstance(value, bool),
            f"{label} must be numeric")
    number = float(value)
    require(math.isfinite(number), f"{label} must be finite")
    return number


def solution_scalar(value: Any, label: str) -> float:
    if isinstance(value, list):
        require(value, f"{label} vector must not be empty")
        components = [finite_number(entry, f"{label}[{index}]")
                      for index, entry in enumerate(value)]
        return math.sqrt(sum(component * component for component in components))
    return finite_number(value, label)


def load_payload(path: Path) -> tuple[list[list[float]], list[list[int]], list[float]]:
    payload = load_json(path)
    require(payload.get("schema") == PAYLOAD_SCHEMA,
            f"payload schema must be {PAYLOAD_SCHEMA}")
    require(payload.get("case_id") == CASE_ID,
            f"payload case_id must be {CASE_ID}")
    require(payload.get("provenance_uri") == PAPER_URI,
            f"payload provenance URI must be {PAPER_URI}")

    geometry = payload.get("geometry")
    require(isinstance(geometry, dict), "payload.geometry must be an object")
    raw_nodes = geometry.get("nodes")
    raw_elements = geometry.get("elements")
    require(isinstance(raw_nodes, list) and len(raw_nodes) >= 3,
            "payload.geometry.nodes must contain at least three nodes")
    require(isinstance(raw_elements, list) and raw_elements,
            "payload.geometry.elements must not be empty")
    nodes: list[list[float]] = []
    for node_index, node in enumerate(raw_nodes):
        require(isinstance(node, list) and len(node) in (2, 3),
                f"geometry.nodes[{node_index}] must have two or three coordinates")
        coordinates = [finite_number(value, f"geometry.nodes[{node_index}]")
                       for value in node]
        if len(coordinates) == 2:
            coordinates.append(0.0)
        nodes.append(coordinates)

    raw_indices: list[list[int]] = []
    for element_index, element in enumerate(raw_elements):
        require(isinstance(element, list) and len(element) >= 3,
                f"geometry.elements[{element_index}] must have at least three nodes")
        indices: list[int] = []
        for local_index, value in enumerate(element):
            require(isinstance(value, int) and not isinstance(value, bool),
                    f"geometry.elements[{element_index}][{local_index}] must be integer")
            indices.append(value)
        raw_indices.append(indices)
    all_indices = [index for element in raw_indices for index in element]
    if any(index == 0 for index in all_indices):
        offset = 0
    else:
        offset = 1
    elements = [[index - offset for index in element] for element in raw_indices]
    for element_index, element in enumerate(elements):
        require(all(0 <= index < len(nodes) for index in element),
                f"geometry.elements[{element_index}] references an invalid node")

    raw_solution = payload.get("solution")
    if isinstance(raw_solution, dict):
        raw_solution = raw_solution.get("values")
    require(isinstance(raw_solution, list) and len(raw_solution) == len(nodes),
            "payload.solution.values must match the node count")
    solution = [solution_scalar(value, f"solution[{index}]")
                for index, value in enumerate(raw_solution)]
    return nodes, elements, solution


def color(value: float, lower: float, upper: float) -> str:
    span = max(upper - lower, 1.0e-15)
    t = min(1.0, max(0.0, (value - lower) / span))
    red = int(round(255.0 * t))
    green = int(round(100.0 + 120.0 * (1.0 - abs(2.0 * t - 1.0))))
    blue = int(round(255.0 * (1.0 - t)))
    return f"#{red:02x}{green:02x}{blue:02x}"


def write_gallery(
    output_dir: Path, nodes: list[list[float]], elements: list[list[int]],
    solution: list[float], artifact: Path, digest: str,
) -> None:
    projected = [(node[0] - 0.45 * node[1], node[2] + 0.30 * node[1])
                 for node in nodes]
    min_u = min(point[0] for point in projected)
    max_u = max(point[0] for point in projected)
    min_v = min(point[1] for point in projected)
    max_v = max(point[1] for point in projected)
    span = max(max_u - min_u, max_v - min_v, 1.0e-12)

    def screen(point: tuple[float, float]) -> tuple[float, float]:
        return (90.0 + 760.0 * (point[0] - min_u) / span,
                610.0 - 520.0 * (point[1] - min_v) / span)

    lower = min(solution)
    upper = max(solution)
    polygons: list[str] = []
    for element in elements:
        vertices = [screen(projected[index]) for index in element]
        average = sum(solution[index] for index in element) / len(element)
        points = " ".join(f"{x:.3f},{y:.3f}" for x, y in vertices)
        polygons.append(
            f'<polygon points="{points}" fill="{color(average, lower, upper)}" '
            'stroke="#25324a" stroke-width="1.2" fill-opacity="0.88"/>')
    circles = []
    for index, point in enumerate(projected):
        x, y = screen(point)
        circles.append(
            f'<circle cx="{x:.3f}" cy="{y:.3f}" r="3.0" '
            f'fill="{color(solution[index], lower, upper)}" stroke="#101820" '
            f'stroke-width="0.7"/>')
    title = "Bíro 1996 exact-data solution (verified)"
    subtitle = "solution-first • tree/cotree magnetostatic vector-potential payload"
    svg = [
        '<svg xmlns="http://www.w3.org/2000/svg" width="960" height="700" '
        'viewBox="0 0 960 700">',
        '<rect width="960" height="700" fill="#ffffff"/>',
        '<text x="60" y="38" font-family="sans-serif" font-size="22" '
        'font-weight="600" fill="#17233b">' + html.escape(title) + '</text>',
        '<text x="60" y="64" font-family="sans-serif" font-size="13" '
        'fill="#394865">' + html.escape(subtitle) + '</text>',
        *polygons,
        *circles,
        '<line x1="885" y1="130" x2="885" y2="610" stroke="#202b3d" '
        'stroke-width="1"/>',
    ]
    for step in range(11):
        value = lower + (upper - lower) * step / 10.0
        y = 610.0 - 480.0 * step / 10.0
        svg.append(f'<line x1="880" y1="{y:.1f}" x2="890" y2="{y:.1f}" '
                   'stroke="#202b3d" stroke-width="1"/>')
        svg.append(f'<text x="898" y="{y + 4.0:.1f}" font-family="sans-serif" '
                   f'font-size="11" fill="#202b3d">{value:.5g}</text>')
    svg.extend([
        '<text x="850" y="112" font-family="sans-serif" font-size="12" '
        'fill="#202b3d">|A| / solution</text>',
        '<text x="60" y="658" font-family="sans-serif" font-size="11" '
        'fill="#52617a">DOI: 10.1109/20.497322 • SHA-256: ' + digest + '</text>',
        '<text x="60" y="678" font-family="sans-serif" font-size="11" '
        'fill="#52617a">artifact: ' + html.escape(str(artifact.name)) +
        ' • paper arrays remain external</text>',
        '</svg>',
    ])
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "solution.svg").write_text("\n".join(svg) + "\n", encoding="utf-8")
    with (output_dir / "solution.csv").open("w", encoding="utf-8") as stream:
        stream.write("node,x,y,z,solution\n")
        for index, node in enumerate(nodes):
            stream.write(f"{index},{node[0]:.16e},{node[1]:.16e},"
                         f"{node[2]:.16e},{solution[index]:.16e}\n")
    metadata = {
        "status": "ready",
        "case_id": CASE_ID,
        "provenance_uri": PAPER_URI,
        "exact_data": True,
        "manufactured": False,
        "solution_first": True,
        "plot": "solution.svg",
        "artifact_sha256": digest,
        "node_count": len(nodes),
        "element_count": len(elements),
    }
    (output_dir / "provenance.json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8")


def main() -> int:
    args = parser().parse_args()
    try:
        contract = load_json(args.contract)
        validate_manifest(contract)
        artifact, _ = verified_artifact(contract, args.data_manifest, args.data_root)
        if artifact is None:
            print("SKIP: exact Biro paper payload is absent; no solution plot generated")
            return 0
        nodes, elements, solution = load_payload(artifact)
        digest = canonical_file_sha256(artifact)
        write_gallery(args.output_dir, nodes, elements, solution, artifact, digest)
        print(f"SOLUTION_FIRST: wrote verified exact Biro gallery to {args.output_dir}")
        print(f"               plot={args.output_dir / 'solution.svg'}; sha256={digest}")
        return 0
    except (OSError, json.JSONDecodeError, ContractError) as error:
        print(f"Biro solution gallery: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
